///////////////////////////////////////////////////////////////////////////////
// File: HITSpectrum.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Turbulence/HITSpectrum.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <locale>
#include <regex>
#include <sstream>
#include <stdexcept>

namespace Nektar::SolverUtils
{
namespace
{
void Require(bool condition, const std::string &message)
{
    if (!condition)
    {
        throw std::invalid_argument("HITSpectrum: " + message);
    }
}

/// log(a+b) without overflowing the sum of positive finite numbers.
long double LogSum(NekDouble a, NekDouble b)
{
    NekDouble largest = std::max(a, b);
    NekDouble smallest = std::min(a, b);
    return std::log(static_cast<long double>(largest)) +
           std::log1p(static_cast<long double>(smallest) / largest);
}

/// log(a/b) for a>b>0. log1p retains the difference of adjacent values;
/// separate logarithms handle ratios that would overflow or underflow.
long double LogRatio(NekDouble a, NekDouble b)
{
    NekDouble difference = a - b;
    if (difference <= b)
    {
        return std::log1p(static_cast<long double>(difference) / b);
    }
    return std::log(static_cast<long double>(a)) -
           std::log(static_cast<long double>(b));
}

NekDouble SpectrumExp(long double logarithm)
{
    const long double maxLog = std::log(
        static_cast<long double>(std::numeric_limits<NekDouble>::max()));
    const long double minLog = std::log(
        static_cast<long double>(std::numeric_limits<NekDouble>::denorm_min()));
    if (logarithm < minLog)
    {
        return 0.0;
    }
    if (!std::isfinite(logarithm) || logarithm > maxLog)
    {
        throw std::overflow_error("HITSpectrum: spectrum exceeds numerical range");
    }
    long double value = std::exp(logarithm);
    if (!std::isfinite(value) ||
        value > std::numeric_limits<NekDouble>::max())
    {
        throw std::overflow_error("HITSpectrum: spectrum exceeds numerical range");
    }
    return static_cast<NekDouble>(value);
}

/// For expressions mathematically bounded by an input value, clamp roundoff
/// at that bound. Underflow below the smallest subnormal is permitted.
NekDouble BoundedExp(long double logarithm, NekDouble upper)
{
    const long double minLog = std::log(
        static_cast<long double>(std::numeric_limits<NekDouble>::denorm_min()));
    if (logarithm < minLog)
    {
        return 0.0;
    }
    return static_cast<NekDouble>(
        std::min(std::exp(logarithm), static_cast<long double>(upper)));
}

std::string Trim(const std::string &text)
{
    const std::string whitespace = " \t\r\n\f\v";
    auto first = text.find_first_not_of(whitespace);
    if (first == std::string::npos)
    {
        return "";
    }
    return text.substr(first, text.find_last_not_of(whitespace) - first + 1);
}

NekDouble ParseNumber(const std::string &text, std::size_t line)
{
    static const std::regex decimal(
        R"([+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?)");
    Require(std::regex_match(text, decimal),
            "invalid decimal token on line " + std::to_string(line));
    std::istringstream input(text);
    input.imbue(std::locale::classic());
    NekDouble value = 0.0;
    input >> value;
    Require(!input.fail() && input.eof() && std::isfinite(value),
            "decimal outside finite numerical range on line " +
                std::to_string(line));
    return value;
}
} // namespace

HITSpectrum HITSpectrum::VonKarman(NekDouble length, NekDouble amplitude)
{
    Require(std::isfinite(length) && length > 0.0,
            "von Karman length must be finite and positive");
    Require(std::isfinite(amplitude) && amplitude >= 0.0,
            "von Karman amplitude must be finite and nonnegative");
    HITSpectrum result;
    result.m_length = length;
    result.m_amplitude = amplitude;
    result.m_identicallyZero = amplitude == 0.0;
    return result;
}

HITSpectrum HITSpectrum::Tabulated(
    const std::vector<NekDouble> &wavenumbers,
    const std::vector<NekDouble> &energies, HITInterpolation interpolation,
    HITLowWavenumber low, HITHighWavenumber high)
{
    Require(interpolation == HITInterpolation::Linear ||
                interpolation == HITInterpolation::LogLinear,
            "unknown interpolation policy");
    Require(low == HITLowWavenumber::Zero || low == HITLowWavenumber::K4,
            "unknown low-wavenumber policy");
    Require(high == HITHighWavenumber::Zero ||
                high == HITHighWavenumber::Gaussian,
            "unknown high-wavenumber policy");
    Require(wavenumbers.size() >= 2 && wavenumbers.size() == energies.size(),
            "table requires at least two matching k and E rows");
    for (std::size_t i = 0; i < wavenumbers.size(); ++i)
    {
        Require(std::isfinite(wavenumbers[i]) && wavenumbers[i] > 0.0,
                "table wavenumbers must be finite and positive");
        Require(i == 0 || wavenumbers[i] > wavenumbers[i - 1],
                "table wavenumbers must strictly increase");
        Require(std::isfinite(energies[i]) && energies[i] >= 0.0,
                "table energies must be finite and nonnegative");
        Require(interpolation != HITInterpolation::LogLinear ||
                    energies[i] > 0.0,
                "LogLinear requires strictly positive energies");
    }

    HITSpectrum result;
    result.m_tabulated = true;
    result.m_wavenumbers = wavenumbers;
    result.m_energies = energies;
    result.m_interpolation = interpolation;
    result.m_low = low;
    result.m_high = high;
    result.m_identicallyZero =
        std::all_of(energies.begin(), energies.end(),
                    [](NekDouble energy) { return energy == 0.0; });
    if (high == HITHighWavenumber::Gaussian)
    {
        std::size_t last = energies.size() - 1;
        Require(energies[last - 1] > energies[last] && energies[last] > 0.0,
                "Gaussian tail requires a decreasing positive final pair");
        result.m_gaussianLogBeta = static_cast<NekDouble>(
            std::log(LogRatio(energies[last - 1], energies[last])) -
            std::log(static_cast<long double>(wavenumbers[last] -
                                              wavenumbers[last - 1])) -
            LogSum(wavenumbers[last], wavenumbers[last - 1]));
    }
    return result;
}

HITSpectrum HITSpectrum::ReadTable(std::istream &input,
                                  HITInterpolation interpolation,
                                  HITLowWavenumber low,
                                  HITHighWavenumber high)
{
    const std::string expectedHeader =
        "# NEKTAR_HIT_SPECTRUM_V1 KIND=Shell3D K_UNITS=rad/m E_UNITS=m3/s2";
    bool header = false;
    std::vector<NekDouble> wavenumbers, energies;
    std::string line;
    std::size_t lineNumber = 0;
    while (std::getline(input, line))
    {
        ++lineNumber;
        line = Trim(line);
        if (line.empty())
        {
            continue;
        }
        if (!header)
        {
            Require(line == expectedHeader,
                    "first nonblank line must declare the version-1 SI "
                    "Shell3D spectrum");
            header = true;
            continue;
        }
        if (line[0] == '#')
        {
            continue;
        }
        std::istringstream row(line);
        row.imbue(std::locale::classic());
        std::string k, energy, extra;
        Require(static_cast<bool>(row >> k >> energy) && !(row >> extra),
                "expected exactly two columns on line " +
                    std::to_string(lineNumber));
        wavenumbers.push_back(ParseNumber(k, lineNumber));
        energies.push_back(ParseNumber(energy, lineNumber));
    }
    Require(!input.bad(), "failed while reading spectrum input");
    Require(header, "missing spectrum header");
    return Tabulated(wavenumbers, energies, interpolation, low, high);
}

HITSpectrum HITSpectrum::FromFile(const std::string &filename,
                                 HITInterpolation interpolation,
                                 HITLowWavenumber low,
                                 HITHighWavenumber high)
{
    std::ifstream input(filename);
    Require(input.is_open(), "cannot open spectrum file '" + filename + "'");
    return ReadTable(input, interpolation, low, high);
}

NekDouble HITSpectrum::Evaluate(NekDouble k) const
{
    Require(std::isfinite(k) && k >= 0.0,
            "evaluation wavenumber must be finite and nonnegative");
    if (k == 0.0 || IsIdenticallyZero())
    {
        return 0.0;
    }
    if (!m_tabulated)
    {
        long double logL = std::log(static_cast<long double>(m_length));
        long double logQ = logL + std::log(static_cast<long double>(k));
        long double twiceLogQ = 2.0L * logQ;
        long double logDenominator =
            twiceLogQ > 0.0L
                ? twiceLogQ + std::log1p(std::exp(-twiceLogQ))
                : std::log1p(std::exp(twiceLogQ));
        return SpectrumExp(std::log(static_cast<long double>(m_amplitude)) +
                           (5.0L / 3.0L) * logL + 4.0L * logQ -
                           (17.0L / 6.0L) * logDenominator);
    }

    auto position = std::lower_bound(m_wavenumbers.begin(),
                                     m_wavenumbers.end(), k);
    if (position != m_wavenumbers.end() && *position == k)
    {
        return m_energies[position - m_wavenumbers.begin()];
    }
    if (position == m_wavenumbers.begin())
    {
        if (m_low == HITLowWavenumber::Zero || m_energies.front() == 0.0)
        {
            return 0.0;
        }
        return BoundedExp(
            std::log(static_cast<long double>(m_energies.front())) -
                4.0L * LogRatio(m_wavenumbers.front(), k),
            m_energies.front());
    }
    if (position == m_wavenumbers.end())
    {
        if (m_high == HITHighWavenumber::Zero)
        {
            return 0.0;
        }
        NekDouble lastK = m_wavenumbers.back();
        long double logDecay =
            m_gaussianLogBeta +
            std::log(static_cast<long double>(k - lastK)) + LogSum(k, lastK);
        // exp(10) already exceeds the log-ratio of the largest finite value
        // to the smallest subnormal, so the tail is necessarily zero here.
        if (logDecay > 10.0L)
        {
            return 0.0;
        }
        return BoundedExp(
            std::log(static_cast<long double>(m_energies.back())) -
                std::exp(logDecay),
            m_energies.back());
    }

    std::size_t right = position - m_wavenumbers.begin();
    std::size_t left = right - 1;
    NekDouble eLeft = m_energies[left], eRight = m_energies[right];
    if (eLeft == eRight)
    {
        return eLeft;
    }
    if (m_interpolation == HITInterpolation::Linear)
    {
        NekDouble fraction = (k - m_wavenumbers[left]) /
                             (m_wavenumbers[right] - m_wavenumbers[left]);
        return eLeft <= eRight
                   ? std::fma(fraction, eRight - eLeft, eLeft)
                   : std::fma(1.0 - fraction, eLeft - eRight, eRight);
    }
    long double fraction = LogRatio(k, m_wavenumbers[left]) /
                           LogRatio(m_wavenumbers[right], m_wavenumbers[left]);
    long double logValue =
        (1.0L - fraction) * std::log(static_cast<long double>(eLeft)) +
        fraction * std::log(static_cast<long double>(eRight));
    NekDouble value = BoundedExp(logValue, std::max(eLeft, eRight));
    return std::max(value, std::min(eLeft, eRight));
}

bool HITSpectrum::IsIdenticallyZero() const
{
    return m_identicallyZero;
}

const std::vector<NekDouble> &HITSpectrum::GetWavenumbers() const
{
    return m_wavenumbers;
}

const std::vector<NekDouble> &HITSpectrum::GetEnergies() const
{
    return m_energies;
}
} // namespace Nektar::SolverUtils
