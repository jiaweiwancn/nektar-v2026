///////////////////////////////////////////////////////////////////////////////
// File: RANSHistory.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSHistory.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace Nektar
{
using LibUtilities::TripleArray;
void RANSIMEX::Restore(const TripleArray &values,
    const Array<OneD, NekDouble> &times, size_t steps)
{
    m_restored = values; m_times = times; m_offset = steps;
}

void RANSIMEX::v_InitializeScheme(NekDouble dt, LibUtilities::ConstDoubleArray &initial,
    NekDouble time, const LibUtilities::TimeIntegrationSchemeOperators &op)
{
    IMEXTimeIntegrationScheme::v_InitializeScheme(dt, initial, time, op);
    if (!m_restored.size()) return;
    // This constructor marks all saved values/derivatives as initialized.
    m_solVector = MemoryManager<LibUtilities::TimeIntegrationSolutionGLM>::AllocateSharedPtr(
        m_integration_phases.back().get(), initial.size(), initial[0].size());
    m_solVector->UpdateSolutionVector() = m_restored;
    m_solVector->UpdateTimeVector() = m_times;
    m_restored = TripleArray();
    m_times = Array<OneD, NekDouble>();
}

LibUtilities::ConstDoubleArray &RANSIMEX::v_TimeIntegrate(size_t step, NekDouble dt)
{
    return IMEXTimeIntegrationScheme::v_TimeIntegrate(step + m_offset, dt);
}

namespace
{
constexpr uint64_t magic = 0x52414e5348495331ULL, endian = 0x0102030405060708ULL;
static_assert(sizeof(NekDouble) == 8, "RANS history schema requires 64-bit NekDouble");
uint64_t Hash(const char *data, size_t size)
{
    uint64_t hash = 14695981039346656037ULL;
    for (size_t i = 0; i < size; ++i) { hash ^= static_cast<unsigned char>(data[i]); hash *= 1099511628211ULL; }
    return hash;
}
template<class T> void Put(std::string &bytes, const T *value, size_t count = 1)
{
    bytes.append(reinterpret_cast<const char *>(value), count*sizeof(T));
}
void Number(std::string &bytes, uint64_t value) { Put(bytes, &value); }
struct Reader
{
    const std::string &data; size_t cursor = 0;
    template<class T> void Get(T *value, size_t count = 1)
    {
        if (count > (data.size()-cursor)/sizeof(T)) throw std::runtime_error("truncated payload");
        std::memcpy(value, data.data()+cursor, count*sizeof(T)); cursor += count*sizeof(T);
    }
    uint64_t Number() { uint64_t n; Get(&n); return n; }
};
void Require(const LibUtilities::CommSharedPtr &comm, bool valid, const std::string &message)
{
    int bad = valid ? 0 : 1; comm->AllReduce(bad, LibUtilities::ReduceMax);
    ASSERTL0(!bad, "invalid RANS history: " + message);
}
uint64_t Basis(const MultiRegions::ExpListSharedPtr &field, int e)
{
    std::ostringstream s; auto exp = field->GetExp(e);
    s << exp->DetShapeType();
    for (int d = 0; d < exp->GetNumBases(); ++d)
        s << ':' << exp->GetBasis(d)->GetBasisType() << ':' << exp->GetBasis(d)->GetNumModes()
          << ':' << exp->GetBasis(d)->GetPointsType() << ':' << exp->GetBasis(d)->GetNumPoints();
    const auto value = s.str(); return Hash(value.data(), value.size());
}
std::string Precise(double value) { std::ostringstream s; s << std::setprecision(17) << value; return s.str(); }
}

std::string RANSHistory::Write(const std::string &baseDirectory, const Fields &fields,
    const RANSIMEX &scheme, const std::string &signature,
    size_t steps, NekDouble time, NekDouble dt)
{
    auto comm = fields[0]->GetComm(); bool valid = true;
    int generation = 0;
    auto candidate = [&]() { return baseDirectory + (generation ? "."+std::to_string(generation) : ""); };
    try
    {
        if (comm->GetRank() == 0)
            while (std::filesystem::exists(candidate())) ++generation;
    }
    catch (...) { valid = false; }
    Require(comm, valid, "cannot select checkpoint directory");
    comm->Bcast(generation,0);
    const std::string directory = candidate();
    try
    {
        // Never overwrite payloads referenced by a native backup field file.
        if (comm->GetRank() == 0 && !std::filesystem::create_directory(directory)) valid = false;
    }
    catch (...) { valid = false; }
    Require(comm, valid, "cannot reserve checkpoint directory"); comm->Block();
    const auto &history = scheme.GetSolutionVector();
    const int nf = fields.size()-1, nq = fields[0]->GetTotPoints();
    Array<OneD, NekDouble> x(nq, 0.0), y(nq, 0.0), z(nq, 0.0);
    fields[0]->GetCoords(x, y, z);
    try
    {
        std::string bytes; Number(bytes, magic); Number(bytes, endian);
        Number(bytes, history.size()); Number(bytes, nf); Number(bytes, fields[0]->GetExpSize());
        for (int e = 0; e < fields[0]->GetExpSize(); ++e)
        {
            const int q = fields[0]->GetExp(e)->GetTotPoints(), offset = fields[0]->GetPhys_Offset(e);
            const int nc = fields[0]->GetExp(e)->GetNcoeffs(), pc = fields[nf]->GetExp(e)->GetNcoeffs();
            Number(bytes, fields[0]->GetExp(e)->GetGeom()->GetGlobalID());
            Number(bytes, q); Number(bytes, nc); Number(bytes, pc);
            Number(bytes, Basis(fields[0], e)); Number(bytes, Basis(fields[nf], e));
            for (const auto &v : {x, y, z}) Put(bytes, v.data()+offset, q);
            for (int f = 0; f < fields.size(); ++f)
                Put(bytes, fields[f]->GetCoeffs().data()+fields[f]->GetCoeff_Offset(e), f == nf ? pc : nc);
            for (int slot = 0; slot < history.size(); ++slot)
                for (int f = 0; f < nf; ++f) Put(bytes, history[slot][f].data()+offset, q);
        }
        Number(bytes, Hash(bytes.data(), bytes.size()));
        const auto file = std::filesystem::path(directory)/("rank-"+std::to_string(comm->GetRank())+".bin");
        std::ofstream stream(file.string()+".tmp", std::ios::binary);
        stream.write(bytes.data(), bytes.size()); stream.close();
        if (!stream) throw std::runtime_error("write failed");
        std::filesystem::rename(file.string()+".tmp", file);
    }
    catch (...) { valid = false; }
    Require(comm, valid, "cannot publish rank payload"); comm->Block();
    if (comm->GetRank() == 0)
    {
        try
        {
            TiXmlDocument document; auto node = new TiXmlElement("RANSHISTORY"); document.LinkEndChild(node);
            node->SetAttribute("SCHEMA", 1); node->SetAttribute("RANKS", comm->GetSize());
            node->SetAttribute("STEPS", std::to_string(steps)); node->SetAttribute("ORDER", static_cast<int>(scheme.GetOrder()));
            node->SetAttribute("TIME", Precise(time)); node->SetAttribute("DT", Precise(dt));
            node->SetAttribute("SIGNATURE", signature);
            // Native one-stage IMEX updates/uses only time[0]; its other
            // time-vector entries are unspecified. Store canonical finite
            // labels instead of serializing those unused memory locations.
            for (int slot = 0; slot < history.size(); ++slot)
            {
                const double value = slot < static_cast<int>(scheme.GetOrder()) ? time-slot*dt : dt;
                auto entry = new TiXmlElement("T"); entry->SetAttribute("VALUE", Precise(value)); node->LinkEndChild(entry);
            }
            const auto file = std::filesystem::path(directory)/"Info.xml";
            if (!document.SaveFile((file.string()+".tmp").c_str())) throw std::runtime_error("manifest write failed");
            std::filesystem::rename(file.string()+".tmp", file);
        }
        catch (...) { valid = false; }
    }
    Require(comm, valid, "cannot publish history manifest"); comm->Block();
    return directory;
}

size_t RANSHistory::Read(const std::string &directory, const Fields &fields,
    const LibUtilities::DoubleArray &importedCoefficients,
    RANSIMEX &scheme, const std::string &signature, NekDouble time, NekDouble dt)
{
    const auto comm = fields[0]->GetComm(); const int nf = fields.size()-1;
    const int slots = 2*scheme.GetOrder(), nq = fields[0]->GetTotPoints();
    TripleArray history(slots);
    for (auto &slot : history)
    {
        slot = LibUtilities::DoubleArray(nf);
        for (auto &v : slot) v = Array<OneD, NekDouble>(nq, 0.0);
    }
    Array<OneD, NekDouble> times(slots), x(nq, 0.0), y(nq, 0.0), z(nq, 0.0);
    fields[0]->GetCoords(x, y, z);
    Array<OneD, Array<OneD, NekDouble>> coefficients(fields.size());
    for (int f = 0; f < fields.size(); ++f) coefficients[f] = Array<OneD, NekDouble>(fields[f]->GetNcoeffs());
    std::map<int, int> owned; std::set<int> seen, loaded;
    for (int e = 0; e < fields[0]->GetExpSize(); ++e) owned[fields[0]->GetExp(e)->GetGeom()->GetGlobalID()] = e;
    int total = fields[0]->GetExpSize(); comm->AllReduce(total, LibUtilities::ReduceSum);
    bool valid = true; std::string error; size_t steps = 0;
    try
    {
        TiXmlDocument document((std::filesystem::path(directory)/"Info.xml").string());
        if (!document.LoadFile()) throw std::runtime_error("missing manifest");
        auto root = document.FirstChildElement("RANSHISTORY");
        if (!root) throw std::runtime_error("missing manifest root");
        auto attribute = [&](const char *key) {
            const char *v = root->Attribute(key); if (!v) throw std::runtime_error("missing manifest attribute"); return std::string(v);
        };
        if (attribute("SCHEMA") != "1" || std::stoi(attribute("ORDER")) != static_cast<int>(scheme.GetOrder()) ||
            attribute("SIGNATURE") != signature) throw std::runtime_error("model/configuration/order mismatch");
        if (std::abs(std::stod(attribute("TIME"))-time) > 1e-12*std::max(1.0,std::abs(time)) ||
            std::abs(std::stod(attribute("DT"))-dt) > 1e-13*dt) throw std::runtime_error("time or timestep mismatch");
        steps = std::stoull(attribute("STEPS"));
        if (steps == 0) throw std::runtime_error("uninitialized history");
        int slot = 0;
        for (auto item = root->FirstChildElement("T"); item; item = item->NextSiblingElement("T"))
        {
            if (slot >= slots || !item->Attribute("VALUE")) throw std::runtime_error("bad time vector");
            times[slot++] = std::stod(item->Attribute("VALUE"));
        }
        if (slot != slots) throw std::runtime_error("incomplete time vector");
        for (int j = 0; j < slots; ++j)
        {
            const double expected = j < static_cast<int>(scheme.GetOrder()) ? time-j*dt : dt;
            if (!std::isfinite(times[j]) || std::abs(times[j]-expected) > 1e-12*std::max(1.0,std::abs(expected)))
                throw std::runtime_error("invalid time labels");
        }
        const int ranks = std::stoi(attribute("RANKS"));
        if (ranks <= 0 || ranks > total) throw std::runtime_error("invalid rank count");
        for (int rank = 0; rank < ranks; ++rank)
        {
            const auto file = std::filesystem::path(directory)/("rank-"+std::to_string(rank)+".bin");
            std::ifstream stream(file, std::ios::binary);
            if (!stream) throw std::runtime_error("missing rank payload");
            std::string bytes((std::istreambuf_iterator<char>(stream)), std::istreambuf_iterator<char>());
            if (bytes.size() < 6*sizeof(uint64_t)) throw std::runtime_error("short rank payload");
            uint64_t checksum; std::memcpy(&checksum, bytes.data()+bytes.size()-sizeof(checksum), sizeof(checksum));
            if (checksum != Hash(bytes.data(),bytes.size()-sizeof(checksum))) throw std::runtime_error("payload checksum mismatch");
            Reader reader{bytes};
            if (reader.Number() != magic || reader.Number() != endian || reader.Number() != static_cast<uint64_t>(slots) ||
                reader.Number() != static_cast<uint64_t>(nf)) throw std::runtime_error("payload schema/architecture mismatch");
            const auto count = reader.Number();
            if (count > static_cast<uint64_t>(total)) throw std::runtime_error("invalid element count");
            for (uint64_t n = 0; n < count; ++n)
            {
                const auto id = reader.Number(), q = reader.Number(), nc = reader.Number(), pc = reader.Number();
                const auto basis = reader.Number(), pressureBasis = reader.Number();
                if (!seen.insert(id).second) throw std::runtime_error("duplicate element");
                if (q > bytes.size()/sizeof(double) || nc > bytes.size()/sizeof(double) || pc > bytes.size()/sizeof(double))
                    throw std::runtime_error("invalid payload dimensions");
                const auto owner = owned.find(id); const bool here = owner != owned.end();
                int e = here ? owner->second : 0, offset = here ? fields[0]->GetPhys_Offset(e) : 0;
                if (here && (q != static_cast<uint64_t>(fields[0]->GetExp(e)->GetTotPoints()) ||
                    nc != static_cast<uint64_t>(fields[0]->GetExp(e)->GetNcoeffs()) ||
                    pc != static_cast<uint64_t>(fields[nf]->GetExp(e)->GetNcoeffs()) ||
                    basis != Basis(fields[0],e) || pressureBasis != Basis(fields[nf],e)))
                    throw std::runtime_error("mesh, basis or quadrature mismatch");
                Array<OneD, NekDouble> buffer(std::max({q,nc,pc}));
                for (const auto &coordinate : {x,y,z})
                {
                    reader.Get(buffer.data(),q);
                    if (here) for (int i = 0; i < q; ++i)
                        if (std::abs(buffer[i]-coordinate[offset+i]) > 1e-12*std::max(1.0,std::abs(coordinate[offset+i])))
                            throw std::runtime_error("element coordinates/order mismatch");
                }
                for (int f = 0; f < fields.size(); ++f)
                {
                    const auto size = f == nf ? pc : nc; reader.Get(buffer.data(),size);
                    if (here)
                    {
                        const int start = fields[f]->GetCoeff_Offset(e);
                        for (int i = 0; i < size; ++i)
                        {
                            if (!std::isfinite(buffer[i]) || std::abs(buffer[i]-importedCoefficients[f][start+i]) >
                                1e-8*std::max(1.0,std::abs(buffer[i]))) throw std::runtime_error("field/history state mismatch: field " +
                                    std::to_string(f) + " element " + std::to_string(id));
                            coefficients[f][start+i] = buffer[i];
                        }
                    }
                }
                for (int j = 0; j < slots; ++j)
                    for (int f = 0; f < nf; ++f)
                    {
                        reader.Get(buffer.data(),q);
                        if (here) for (int i = 0; i < q; ++i)
                        {
                            if (!std::isfinite(buffer[i])) throw std::runtime_error("nonfinite history value");
                            history[j][f][offset+i] = buffer[i];
                        }
                    }
                if (here) loaded.insert(id);
            }
            if (reader.cursor + sizeof(uint64_t) != bytes.size()) throw std::runtime_error("unexpected trailing data");
        }
        if (seen.size() != static_cast<size_t>(total) || loaded.size() != owned.size()) throw std::runtime_error("incomplete mesh coverage");
    }
    catch (const std::exception &e) { valid = false; error = e.what(); }
    Require(comm, valid, error);
    for (int f = 0; f < fields.size(); ++f)
    {
        Vmath::Vcopy(coefficients[f].size(), coefficients[f], 1, fields[f]->UpdateCoeffs(), 1);
        if (f < nf) Vmath::Vcopy(nq, history[0][f], 1, fields[f]->UpdatePhys(), 1);
        else fields[f]->BwdTrans(fields[f]->GetCoeffs(), fields[f]->UpdatePhys());
        fields[f]->SetPhysState(true);
    }
    scheme.Restore(history,times,steps);
    return steps;
}
}
