///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallTimeIntegration.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <utility>

namespace Nektar
{
namespace
{
using Vector = RANSWallTimeIntegration::Vector;
using Snapshot = RANSWallTimeIntegration::Snapshot;
using Checkpoint = RANSWallTimeIntegration::Checkpoint;
constexpr std::uint64_t magic = 0x4550415746544931ULL;
constexpr std::uint64_t endian = 0x0102030405060708ULL;
constexpr std::uint64_t schema = 1;
static_assert(sizeof(NekDouble) == 8 && sizeof(double) == 8,
              "Wall IMEX checkpoints require 64-bit doubles");

void Require(bool valid, const std::string &message)
{
    if (!valid)
    {
        throw std::runtime_error("Wall IMEX: " + message);
    }
}

void Finite(const Vector &value, const std::string &name)
{
    Require(std::all_of(value.begin(), value.end(),
                        [](double x) { return std::isfinite(x); }),
            "nonfinite " + name);
}

bool SameTime(double a, double b)
{
    return std::isfinite(a) && std::isfinite(b) &&
           std::abs(a - b) <= 1e-12 * std::max({1., std::abs(a), std::abs(b)});
}

void Validate(const Snapshot &saved, std::size_t order, double dt,
              std::size_t expectedSize)
{
    Require(order == 1 || order == 2, "only native IMEX orders 1 and 2 are supported");
    Require(saved.order == order, "checkpoint order mismatch");
    Require(std::isfinite(dt) && dt > 0 && std::isfinite(saved.dt) &&
                std::abs(saved.dt - dt) <= 1e-13 * dt,
            "checkpoint timestep mismatch");
    Require(std::isfinite(saved.time), "nonfinite checkpoint time");
    Require(expectedSize > 0 && saved.state.size() == expectedSize,
            "checkpoint differential shape mismatch");
    Require(saved.history.size() == 2 * order && saved.times.size() == 2 * order,
            "incomplete native IMEX history");
    Finite(saved.state, "checkpoint differential state");
    for (std::size_t slot = 0; slot < saved.history.size(); ++slot)
    {
        Require(saved.history[slot].size() == expectedSize,
                "checkpoint history shape mismatch");
        Finite(saved.history[slot], "native history");
        const double label = slot < order ? saved.time - slot * dt : dt;
        Require(SameTime(saved.times[slot], label), "checkpoint time labels mismatch");
    }
    Require(saved.history[0] == saved.state,
            "checkpoint current state differs from native history");
}

void Validate(const Checkpoint &saved)
{
    Validate(saved.integration, saved.integration.order, saved.integration.dt,
             saved.integration.state.size());
    Require(!saved.signature.empty(), "empty checkpoint signature");
    Require(saved.coupledState.size() >= saved.integration.state.size(),
            "checkpoint coupled shape mismatch");
    Finite(saved.coupledState, "checkpoint coupled state");
    Finite(saved.fixedReference, "fixed reference");
    Finite(saved.fixedDefect, "fixed reference defect");
    Require(saved.eta == 0 || saved.eta == 1,
            "checkpoint eta must be the declared global zero or one");
    Require(saved.fixedReference.empty() ||
                saved.fixedReference.size() == saved.coupledState.size(),
            "checkpoint fixed-reference shape mismatch");
    Require(saved.fixedDefect.empty() ||
                saved.fixedDefect.size() == saved.coupledState.size() ||
                saved.fixedDefect.size() == saved.integration.state.size(),
            "checkpoint fixed-defect shape mismatch");
    Require(saved.eta == 0 ||
                (!saved.fixedReference.empty() && !saved.fixedDefect.empty()),
            "anchored checkpoint is missing its fixed reference or defect");
}

std::uint64_t Hash(const char *data, std::size_t size)
{
    std::uint64_t hash = 14695981039346656037ULL;
    for (std::size_t i = 0; i < size; ++i)
    {
        hash ^= static_cast<unsigned char>(data[i]);
        hash *= 1099511628211ULL;
    }
    return hash;
}

template <class T> void Put(std::string &bytes, const T &value)
{
    bytes.append(reinterpret_cast<const char *>(&value), sizeof(T));
}

void Number(std::string &bytes, std::uint64_t value) { Put(bytes, value); }

void PutVector(std::string &bytes, const Vector &value)
{
    Number(bytes, value.size());
    if (!value.empty())
    {
        bytes.append(reinterpret_cast<const char *>(value.data()),
                     value.size() * sizeof(double));
    }
}

struct Reader
{
    const std::string &bytes;
    std::size_t cursor = 0;

    template <class T> T Get()
    {
        Require(cursor <= bytes.size() && sizeof(T) <= bytes.size() - cursor,
                "truncated checkpoint");
        T result;
        std::memcpy(&result, bytes.data() + cursor, sizeof(T));
        cursor += sizeof(T);
        return result;
    }

    std::size_t Count(std::size_t limit)
    {
        const auto n = Get<std::uint64_t>();
        Require(n <= limit, "invalid checkpoint count");
        return static_cast<std::size_t>(n);
    }

    Vector GetVector(std::size_t limit)
    {
        const std::size_t n = Count(limit);
        Require(cursor <= bytes.size() && n <= (bytes.size() - cursor) / sizeof(double),
                "truncated checkpoint vector");
        Vector result(n);
        if (n)
        {
            std::memcpy(result.data(), bytes.data() + cursor, n * sizeof(double));
        }
        cursor += n * sizeof(double);
        return result;
    }

    std::string Text(std::size_t limit)
    {
        const std::size_t n = Count(limit);
        Require(cursor <= bytes.size() && n <= bytes.size() - cursor,
                "truncated checkpoint signature");
        const std::string result = bytes.substr(cursor, n);
        cursor += n;
        return result;
    }
};
}

RANSWallTimeIntegration::RANSWallTimeIntegration(
    std::size_t order, double dt, ImplicitCallback implicit,
    ExplicitCallback explicitRhs)
    : m_order(order), m_dt(dt), m_implicit(std::move(implicit)),
      m_explicit(std::move(explicitRhs))
{
    Require(order == 1 || order == 2, "only native IMEX orders 1 and 2 are supported");
    Require(std::isfinite(dt) && dt > 0, "constant timestep must be positive and finite");
    Require(static_cast<bool>(m_implicit) && static_cast<bool>(m_explicit),
            "missing coupled implicit or explicit callback");
}

void RANSWallTimeIntegration::BuildScheme(const Vector &initial, double time,
                                          const Snapshot *restored)
{
    m_scheme = std::make_shared<RANSIMEX>(m_order);
    Operators operators;
    operators.DefineImplicitSolve(&RANSWallTimeIntegration::Implicit, this);
    operators.DefineOdeRhs(&RANSWallTimeIntegration::Explicit, this);
    operators.DefineProjection(&RANSWallTimeIntegration::Projection, this);
    if (restored)
    {
        LibUtilities::TripleArray history(restored->history.size());
        Array<OneD, NekDouble> times(restored->times.size());
        for (std::size_t slot = 0; slot < restored->history.size(); ++slot)
        {
            history[slot] = LibUtilities::DoubleArray(1);
            history[slot][0] = Array<OneD, NekDouble>(m_size, 0.);
            std::copy(restored->history[slot].begin(), restored->history[slot].end(),
                      history[slot][0].data());
            times[slot] = restored->times[slot];
        }
        m_scheme->Restore(history, times, restored->steps);
    }
    LibUtilities::DoubleArray native(1);
    native[0] = Array<OneD, NekDouble>(m_size, 0.);
    std::copy(initial.begin(), initial.end(), native[0].data());
    m_scheme->InitializeScheme(m_dt, native, time, operators);
}

void RANSWallTimeIntegration::Initialise(const Vector &state, double time)
{
    Require(!state.empty() && std::isfinite(time), "invalid initial differential state/time");
    Finite(state, "initial differential state");
    m_size = state.size();
    m_steps = m_localSteps = 0;
    m_time = time;
    m_state = state;
    BuildScheme(state, time, nullptr);
}

RANSWallTimeIntegration::Vector RANSWallTimeIntegration::Advance()
{
    Require(static_cast<bool>(m_scheme), "time integration is not initialized");
    Require(m_steps < std::numeric_limits<std::size_t>::max() &&
                std::isfinite(m_time + m_dt) && m_time + m_dt > m_time,
            "time/step range exceeded");
    const auto &result = m_scheme->TimeIntegrate(m_localSteps, m_dt);
    m_state = Input(result);
    ++m_localSteps;
    ++m_steps;
    m_time += m_dt;
    return m_state;
}

void RANSWallTimeIntegration::Restore(const Snapshot &saved)
{
    Validate(saved, m_order, m_dt, saved.state.size());
    m_size = saved.state.size();
    m_steps = saved.steps;
    m_localSteps = 0;
    m_time = saved.time;
    m_state = saved.state;
    BuildScheme(saved.state, saved.time, &saved);
}

RANSWallTimeIntegration::Snapshot RANSWallTimeIntegration::Capture() const
{
    Require(static_cast<bool>(m_scheme), "time integration is not initialized");
    Snapshot result;
    result.order = m_order;
    result.steps = m_steps;
    result.time = m_time;
    result.dt = m_dt;
    result.state = m_state;
    const auto &native = m_scheme->GetSolutionVector();
    result.history.reserve(native.size());
    result.times.reserve(native.size());
    for (std::size_t slot = 0; slot < native.size(); ++slot)
    {
        Require(native[slot].size() == 1 && native[slot][0].size() == m_size,
                "unexpected native history shape");
        result.history.emplace_back(native[slot][0].data(),
                                    native[slot][0].data() + m_size);
        result.times.push_back(slot < m_order ? m_time - slot * m_dt : m_dt);
    }
    Validate(result, m_order, m_dt, m_size);
    return result;
}

RANSWallTimeIntegration::Vector RANSWallTimeIntegration::Input(
    Operators::InArrayType &input) const
{
    Require(input.size() == 1 && input[0].size() == m_size,
            "unexpected native differential shape");
    Vector result(input[0].data(), input[0].data() + m_size);
    Finite(result, "native differential input");
    return result;
}

void RANSWallTimeIntegration::Output(const Vector &value,
                                     Operators::OutArrayType &output) const
{
    Require(value.size() == m_size && output.size() == 1 &&
                output[0].size() == m_size,
            "coupled callback differential shape mismatch");
    Finite(value, "coupled differential output");
    std::copy(value.begin(), value.end(), output[0].data());
}

void RANSWallTimeIntegration::Implicit(Operators::InArrayType &input,
                                       Operators::OutArrayType &output,
                                       double time, double lambda)
{
    Require(std::isfinite(time) && std::isfinite(lambda) && lambda > 0,
            "invalid native implicit stage time/coefficient");
    Output(m_implicit(Input(input), time, lambda), output);
}

void RANSWallTimeIntegration::Explicit(Operators::InArrayType &input,
                                       Operators::OutArrayType &output,
                                       double time)
{
    Require(std::isfinite(time), "invalid native explicit stage time");
    Output(m_explicit(Input(input), time), output);
}

void RANSWallTimeIntegration::Projection(Operators::InArrayType &input,
                                         Operators::OutArrayType &output,
                                         double time)
{
    Require(std::isfinite(time), "invalid native projection time");
    Output(Input(input), output);
}

std::string RANSWallTimeIntegration::WriteCheckpoint(
    const std::string &baseDirectory, const Checkpoint &saved)
{
    Validate(saved);
    Require(!baseDirectory.empty(), "empty checkpoint directory");
    std::filesystem::path directory;
    for (std::size_t generation = 0; generation < 1000000; ++generation)
    {
        const std::filesystem::path candidate =
            baseDirectory + (generation ? "." + std::to_string(generation) : "");
        std::error_code error;
        if (std::filesystem::create_directory(candidate, error))
        {
            directory = candidate;
            break;
        }
        Require(!error || error == std::errc::file_exists,
                "cannot reserve checkpoint directory: " + error.message());
    }
    Require(!directory.empty(), "checkpoint generation range exceeded");
    std::string bytes;
    Number(bytes, magic);
    Number(bytes, endian);
    Number(bytes, schema);
    Number(bytes, saved.signature.size());
    bytes += saved.signature;
    Number(bytes, saved.integration.order);
    Number(bytes, saved.integration.steps);
    Put(bytes, saved.integration.time);
    Put(bytes, saved.integration.dt);
    Put(bytes, saved.eta);
    PutVector(bytes, saved.integration.state);
    Number(bytes, saved.integration.history.size());
    for (const auto &slot : saved.integration.history)
    {
        PutVector(bytes, slot);
    }
    PutVector(bytes, saved.integration.times);
    PutVector(bytes, saved.coupledState);
    PutVector(bytes, saved.fixedReference);
    PutVector(bytes, saved.fixedDefect);
    Number(bytes, Hash(bytes.data(), bytes.size()));
    const auto temporary = directory / "State.bin.tmp";
    const auto published = directory / "State.bin";
    std::ofstream stream(temporary, std::ios::binary | std::ios::trunc);
    Require(static_cast<bool>(stream), "cannot open checkpoint payload");
    stream.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
    stream.close();
    Require(static_cast<bool>(stream), "cannot write checkpoint payload");
    std::filesystem::rename(temporary, published);
    return directory.string();
}

RANSWallTimeIntegration::Checkpoint RANSWallTimeIntegration::ReadCheckpoint(
    const std::string &directory, const std::string &expectedSignature,
    std::size_t expectedOrder, double expectedDt,
    std::size_t expectedCoupledSize, std::size_t expectedDifferentialSize)
{
    Require(expectedOrder == 1 || expectedOrder == 2, "invalid expected IMEX order");
    Require(expectedDifferentialSize > 0 &&
                expectedCoupledSize >= expectedDifferentialSize &&
                expectedCoupledSize <= 100000000 && !expectedSignature.empty(),
            "invalid expected checkpoint shape/signature");
    const auto file = std::filesystem::path(directory) / "State.bin";
    Require(std::filesystem::is_regular_file(file), "missing checkpoint payload");
    const std::uintmax_t maximum =
        4096 + expectedSignature.size() + sizeof(double) *
            (3 * expectedCoupledSize + (1 + 2 * expectedOrder) *
                                          expectedDifferentialSize);
    const auto length = std::filesystem::file_size(file);
    Require(length >= 10 * sizeof(std::uint64_t) && length <= maximum,
            "invalid checkpoint payload size");
    std::ifstream stream(file, std::ios::binary);
    Require(static_cast<bool>(stream), "cannot read checkpoint payload");
    std::string bytes((std::istreambuf_iterator<char>(stream)),
                      std::istreambuf_iterator<char>());
    Require(bytes.size() == length, "incomplete checkpoint payload read");
    std::uint64_t checksum;
    std::memcpy(&checksum, bytes.data() + bytes.size() - sizeof(checksum),
                sizeof(checksum));
    Require(checksum == Hash(bytes.data(), bytes.size() - sizeof(checksum)),
            "checkpoint checksum mismatch");
    bytes.resize(bytes.size() - sizeof(checksum));
    Reader reader{bytes};
    Require(reader.Get<std::uint64_t>() == magic &&
                reader.Get<std::uint64_t>() == endian &&
                reader.Get<std::uint64_t>() == schema,
            "checkpoint format/schema mismatch");
    Checkpoint result;
    result.signature = reader.Text(bytes.size());
    Require(result.signature == expectedSignature, "checkpoint signature mismatch");
    result.integration.order = reader.Count(2);
    result.integration.steps = reader.Count(std::numeric_limits<std::size_t>::max() - 1);
    result.integration.time = reader.Get<double>();
    result.integration.dt = reader.Get<double>();
    result.eta = reader.Get<double>();
    result.integration.state = reader.GetVector(expectedDifferentialSize);
    const std::size_t slots = reader.Count(2 * expectedOrder);
    result.integration.history.reserve(slots);
    for (std::size_t slot = 0; slot < slots; ++slot)
    {
        result.integration.history.push_back(reader.GetVector(expectedDifferentialSize));
    }
    result.integration.times = reader.GetVector(2 * expectedOrder);
    result.coupledState = reader.GetVector(expectedCoupledSize);
    result.fixedReference = reader.GetVector(expectedCoupledSize);
    result.fixedDefect = reader.GetVector(expectedCoupledSize);
    Require(reader.cursor == bytes.size(), "checkpoint has trailing data");
    Require(result.coupledState.size() == expectedCoupledSize,
            "checkpoint coupled shape mismatch");
    Validate(result.integration, expectedOrder, expectedDt, expectedDifferentialSize);
    Validate(result);
    return result;
}
}
