///////////////////////////////////////////////////////////////////////////////
// File: HITInflowConfig.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////

#include "HITInflowConfig.h"

#include <LibUtilities/Interpreter/Interpreter.h>
#include <tinyxml.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <climits>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Nektar
{
namespace
{
using LibUtilities::CommSharedPtr;
using LibUtilities::SessionReaderSharedPtr;
using SolverUtils::HITGenerator;
using SolverUtils::HITGeneratorConfig;
using SolverUtils::HITHighWavenumber;
using SolverUtils::HITInterpolation;
using SolverUtils::HITLowWavenumber;
using SolverUtils::HITNormalization;
using SolverUtils::HITSpectrum;

[[noreturn]] void Invalid(const std::string &message)
{
    throw std::invalid_argument("HITInflow: " + message);
}

std::string Trim(const std::string &value)
{
    auto first = value.find_first_not_of(" \t\r\n");
    return first == std::string::npos
               ? std::string()
               : value.substr(first, value.find_last_not_of(" \t\r\n") -
                                         first + 1);
}

std::string Upper(std::string value)
{
    for (char &character : value)
    {
        character = static_cast<char>(
            std::toupper(static_cast<unsigned char>(character)));
    }
    return value;
}

std::string XML(const TiXmlElement *element)
{
    TiXmlPrinter printer;
    printer.SetStreamPrinting();
    element->Accept(&printer);
    return printer.CStr();
}

void Attributes(const TiXmlElement *element,
                const std::set<std::string> &required,
                const std::set<std::string> &optional = {})
{
    std::set<std::string> found;
    for (auto attribute = element->FirstAttribute(); attribute;
         attribute = attribute->Next())
    {
        std::string name = attribute->Name();
        if ((!required.count(name) && !optional.count(name)) ||
            !found.insert(name).second)
        {
            Invalid("unknown or duplicate attribute " + name + " on " +
                    element->ValueStr());
        }
        if (Trim(attribute->Value()).empty())
        {
            Invalid("empty attribute " + name + " on " +
                    element->ValueStr());
        }
    }
    for (const auto &name : required)
    {
        if (!found.count(name))
        {
            Invalid("missing " + element->ValueStr() + "/@" + name);
        }
    }
}

std::map<std::string, const TiXmlElement *> Children(
    const TiXmlElement *element, const std::set<std::string> &required)
{
    std::map<std::string, const TiXmlElement *> found;
    for (auto child = element->FirstChild(); child; child = child->NextSibling())
    {
        if (child->ToComment())
        {
            continue;
        }
        if (child->ToText() && Trim(child->Value()).empty())
        {
            continue;
        }
        auto childElement = child->ToElement();
        if (!childElement || !required.count(childElement->ValueStr()) ||
            !found.emplace(childElement->ValueStr(), childElement).second)
        {
            Invalid("unknown, duplicate, or non-element child of " +
                    element->ValueStr());
        }
    }
    for (const auto &name : required)
    {
        if (!found.count(name))
        {
            Invalid("missing " + element->ValueStr() + "/" + name);
        }
    }
    return found;
}

int RegionId(const std::string &text)
{
    std::uint64_t value = HITGenerator::ParseSeed(text);
    if (value > static_cast<std::uint64_t>(INT_MAX))
    {
        Invalid("boundary region REF/ID must lie in 0..INT_MAX");
    }
    return static_cast<int>(value);
}

// Comm's vector interface does not resize receiving buffers, and its count is
// an int. Avoid empty-vector pointers and reject oversized input before Bcast.
void BroadcastString(const CommSharedPtr &comm, std::string &value, int root)
{
    int count = 0;
    if (comm->GetRank() == root)
    {
        count = value.size() > static_cast<std::size_t>(INT_MAX)
                    ? -1
                    : static_cast<int>(value.size());
    }
    comm->Bcast(count, root);
    if (count < 0)
    {
        Invalid("metadata exceeds the communicator's int count limit");
    }
    if (count == 0)
    {
        value.clear();
        return;
    }
    std::vector<char> bytes(count);
    if (comm->GetRank() == root)
    {
        std::copy(value.begin(), value.end(), bytes.begin());
    }
    comm->Bcast(bytes, root);
    value.assign(bytes.begin(), bytes.end());
}

// Every caller enters this routine, including ranks with no inlet elements.
// Choosing the first failing rank gives all participants one identical error.
void AgreeError(const CommSharedPtr &comm, std::string error)
{
    int failedRank = error.empty() ? INT_MAX : comm->GetRank();
    comm->AllReduce(failedRank, LibUtilities::ReduceMin);
    if (failedRank != INT_MAX)
    {
        BroadcastString(comm, error, failedRank);
        throw std::invalid_argument(error);
    }
}

std::string CurrentException()
{
    try
    {
        throw;
    }
    catch (const std::exception &error)
    {
        return error.what();
    }
    catch (const std::string &error)
    {
        return "HITInflow: " + error;
    }
    catch (...)
    {
        return "HITInflow: unknown configuration error";
    }
}

// Validate syntax and identifiers before using Interpreter. Its public API does
// not expose dependencies, and unknown parameters can trigger debug assertions
// rather than exceptions. This follows its expression grammar, function arities
// and constants; evaluation itself still uses the loaded session interpreter.
class ScalarExpression
{
public:
    ScalarExpression(const SessionReaderSharedPtr &session,
                     const TiXmlElement *conditions)
        : m_interpreter(session->GetInterpreter())
    {
        auto parameters = conditions->FirstChildElement("PARAMETERS");
        if (parameters)
        {
            for (auto parameter = parameters->FirstChildElement("P"); parameter;
                 parameter = parameter->NextSiblingElement("P"))
            {
                for (auto node = parameter->FirstChild(); node;
                     node = node->NextSibling())
                {
                    if (node->ToText())
                    {
                        std::string declaration = node->Value();
                        auto equal = declaration.find('=');
                        if (equal != std::string::npos)
                        {
                            // SessionReader retains these exact spellings in
                            // Interpreter, although its own lookup folds case.
                            m_parameters.insert(Trim(declaration.substr(0, equal)));
                        }
                        break;
                    }
                }
            }
        }
    }

    NekDouble Evaluate(const std::string &text, const std::string &setting)
    {
        m_text = Trim(text);
        m_setting = setting;
        m_position = 0;
        if (m_text.empty())
        {
            Fail("empty expression");
        }
        LogicalOr();
        Space();
        if (m_position != m_text.size())
        {
            Fail("invalid expression syntax");
        }
        NekDouble result;
        try
        {
            // Empty vlist disallows all spatial/time variables. The identifier
            // check above also prevents cached x/y/z/t-dependent expressions.
            result = m_interpreter->Evaluate(
                m_interpreter->DefineFunction("", m_text));
        }
        catch (...)
        {
            Invalid(setting + " could not be evaluated: " + CurrentException());
        }
        if (!std::isfinite(result))
        {
            Fail("expression must evaluate to a finite value");
        }
        return result;
    }

    NekDouble Attribute(const TiXmlElement *element, const std::string &name)
    {
        return Evaluate(element->Attribute(name),
                        element->ValueStr() + "/@" + name);
    }

private:
    [[noreturn]] void Fail(const std::string &reason) const
    {
        Invalid(m_setting + ": " + reason + " [" + m_text + "]");
    }

    void Space()
    {
        while (m_position < m_text.size() &&
               std::isspace(static_cast<unsigned char>(m_text[m_position])))
        {
            ++m_position;
        }
    }

    bool Take(const std::string &token)
    {
        Space();
        if (m_text.compare(m_position, token.size(), token) != 0)
        {
            return false;
        }
        m_position += token.size();
        return true;
    }

    void LogicalOr()
    {
        LogicalAnd();
        while (Take("||"))
        {
            LogicalAnd();
        }
    }

    void LogicalAnd()
    {
        Equality();
        while (Take("&&"))
        {
            Equality();
        }
    }

    void Equality()
    {
        Compare();
        while (Take("=="))
        {
            Compare();
        }
    }

    void Compare()
    {
        Add();
        while (Take("<=") || Take(">=") || Take("<") || Take(">"))
        {
            Add();
        }
    }

    void Add()
    {
        Negate();
        while (Take("+") || Take("-"))
        {
            Negate();
        }
    }

    void Negate()
    {
        Take("-");
        Multiply();
    }

    void Multiply()
    {
        Power();
        while (Take("*") || Take("/") || Take("%"))
        {
            Power();
        }
    }

    void Power()
    {
        Base();
        if (Take("^"))
        {
            Power();
        }
    }

    static bool IdentifierCharacter(char value)
    {
        return std::isalnum(static_cast<unsigned char>(value)) || value == '_' ||
               value == '$';
    }

    bool Number()
    {
        std::size_t start = m_position;
        if (m_position < m_text.size() &&
            (m_text[m_position] == '+' || m_text[m_position] == '-'))
        {
            ++m_position;
        }
        bool digits = false;
        while (m_position < m_text.size() &&
               std::isdigit(static_cast<unsigned char>(m_text[m_position])))
        {
            digits = true;
            ++m_position;
        }
        if (m_position < m_text.size() && m_text[m_position] == '.')
        {
            ++m_position;
            while (m_position < m_text.size() &&
                   std::isdigit(static_cast<unsigned char>(m_text[m_position])))
            {
                digits = true;
                ++m_position;
            }
        }
        if (!digits)
        {
            m_position = start;
            return false;
        }
        if (m_position < m_text.size() &&
            (m_text[m_position] == 'e' || m_text[m_position] == 'E'))
        {
            ++m_position;
            if (m_position < m_text.size() &&
                (m_text[m_position] == '+' || m_text[m_position] == '-'))
            {
                ++m_position;
            }
            std::size_t exponent = m_position;
            while (m_position < m_text.size() &&
                   std::isdigit(static_cast<unsigned char>(m_text[m_position])))
            {
                ++m_position;
            }
            if (exponent == m_position)
            {
                Fail("invalid numeric exponent");
            }
        }
        if (m_position < m_text.size() &&
            IdentifierCharacter(m_text[m_position]))
        {
            Fail("invalid numeric literal");
        }
        return true;
    }

    void Base()
    {
        Space();
        if (Take("("))
        {
            LogicalOr();
            if (!Take(")"))
            {
                Fail("unclosed parenthesis");
            }
            return;
        }
        static const std::set<std::string> constants{
            "MEANINGLESS", "E", "LOG2E", "LOG10E", "LN2", "LN10", "PI",
            "PI_2", "PI_4", "1_PI", "2_PI", "2_SQRTPI", "SQRT2",
            "SQRT1_2", "GAMMA", "DEG", "PHI"};
        // Three built-in constant names start with a digit.
        for (const auto &name : {"1_PI", "2_PI", "2_SQRTPI"})
        {
            std::size_t size = std::char_traits<char>::length(name);
            if (m_text.compare(m_position, size, name) == 0 &&
                (m_position + size == m_text.size() ||
                 !IdentifierCharacter(m_text[m_position + size])))
            {
                m_position += size;
                return;
            }
        }
        if (Number())
        {
            return;
        }
        std::size_t start = m_position;
        if (m_position == m_text.size() ||
            (!std::isalpha(static_cast<unsigned char>(m_text[m_position])) &&
             m_text[m_position] != '_' && m_text[m_position] != '$'))
        {
            Fail("expected a number, parameter, or pure function");
        }
        while (m_position < m_text.size() &&
               IdentifierCharacter(m_text[m_position]))
        {
            ++m_position;
        }
        std::string name = m_text.substr(start, m_position - start);
        if (Upper(name) == "X" || Upper(name) == "Y" || Upper(name) == "Z" ||
            Upper(name) == "T")
        {
            Fail("coordinate/time dependence is unsupported");
        }
        Space();
        if (Take("("))
        {
            static const std::map<std::string, int> functions{
                {"abs", 1},   {"asin", 1}, {"acos", 1}, {"atan", 1},
                {"atan2", 2}, {"ang", 2},  {"ceil", 1}, {"cos", 1},
                {"cosh", 1},  {"exp", 1},  {"fabs", 1}, {"floor", 1},
                {"fmax", 2},  {"fmin", 2}, {"fmod", 2}, {"log", 1},
                {"log10", 1}, {"max", 2},  {"min", 2},  {"rad", 2},
                {"sin", 1},   {"sinh", 1}, {"sqrt", 1}, {"tan", 1},
                {"tanh", 1},  {"sign", 1}};
            auto function = functions.find(name);
            if (function == functions.end())
            {
                Fail("unknown or non-pure function " + name +
                     "; random/time-dependent functions are forbidden");
            }
            int arguments = 1;
            LogicalOr();
            while (Take(","))
            {
                ++arguments;
                LogicalOr();
            }
            if (!Take(")") || arguments != function->second)
            {
                Fail("invalid argument count/parentheses for " + name);
            }
        }
        else if (!constants.count(name) && !m_parameters.count(name))
        {
            Fail("unknown parameter " + name +
                 "; parameter spelling is case-sensitive");
        }
    }

    LibUtilities::InterpreterSharedPtr m_interpreter;
    std::set<std::string> m_parameters;
    std::string m_text, m_setting;
    std::size_t m_position = 0;
};

struct LocalXML
{
    bool present = false;
    std::string custom;
    std::string boundaries;
};

LocalXML CollectXML(const SessionReaderSharedPtr &session)
{
    LocalXML result;
    auto root = session->GetDocument().FirstChildElement("NEKTAR");
    if (!root)
    {
        Invalid("missing NEKTAR root");
    }
    for (auto element = root->FirstChildElement(); element;
         element = element->NextSiblingElement())
    {
        if (Upper(element->ValueStr()) == "HITINFLOW")
        {
            result.present = true;
            if (element->ValueStr() != "HITINFLOW" || !result.custom.empty())
            {
                Invalid("require exactly one uppercase NEKTAR/HITINFLOW block");
            }
            result.custom = XML(element);
        }
    }
    TiXmlElement gathered("BOUNDARIES");
    auto conditions = root->FirstChildElement("CONDITIONS");
    auto boundaries =
        conditions ? conditions->FirstChildElement("BOUNDARYCONDITIONS") : nullptr;
    if (boundaries)
    {
        std::set<int> localHitRegions;
        for (auto region = boundaries->FirstChildElement("REGION"); region;
             region = region->NextSiblingElement("REGION"))
        {
            bool hit = false;
            for (auto bc = region->FirstChildElement(); bc;
                 bc = bc->NextSiblingElement())
            {
                const char *type = bc->Attribute("USERDEFINEDTYPE");
                hit = hit || (type && Upper(type) == "HITINFLOW");
            }
            if (hit)
            {
                result.present = true;
                const char *ref = region->Attribute("REF");
                if (!ref)
                {
                    Invalid("HIT-tagged boundary REGION requires REF");
                }
                if (!localHitRegions.insert(RegionId(ref)).second)
                {
                    Invalid("duplicate HIT-tagged boundary region REF");
                }
                gathered.InsertEndChild(*region);
            }
        }
    }
    // Leave ordinary sessions entirely to their existing boundary parser.
    // Custom blocks survive mesh partitioning, even on ranks without an inlet.
    auto regions = conditions ? conditions->FirstChildElement("BOUNDARYREGIONS")
                              : nullptr;
    if (result.present && regions)
    {
        std::set<int> localBoundaryRegions;
        for (auto region = regions->FirstChildElement("B"); region;
             region = region->NextSiblingElement("B"))
        {
            const char *id = region->Attribute("ID");
            if (!id)
            {
                Invalid("boundary region requires ID");
            }
            int regionId = RegionId(id);
            if (!localBoundaryRegions.insert(regionId).second)
            {
                Invalid("duplicate boundary region ID");
            }
            TiXmlElement entry("B");
            entry.SetAttribute("ID", regionId);
            gathered.InsertEndChild(entry);
        }
    }
    result.boundaries = XML(&gathered);
    return result;
}

NekDouble Parameter(const SessionReaderSharedPtr &session,
                    const std::string &name, NekDouble fallback)
{
    NekDouble result = session->DefinesParameter(name)
                           ? session->GetParameter(name)
                           : fallback;
    if (!std::isfinite(result))
    {
        Invalid("parameter " + name + " must be finite");
    }
    return result;
}

void SolverInfo(const SessionReaderSharedPtr &session, const std::string &name,
                const std::string &supported, const std::string &fallback = "")
{
    std::string value = session->DefinesSolverInfo(name)
                            ? session->GetSolverInfo(name)
                            : fallback;
    if (Upper(value) != Upper(supported))
    {
        Invalid("require " + name + "=" + supported + " (received " + value +
                ")");
    }
}

int GridCount(ScalarExpression &scalar, const TiXmlElement *element,
              const std::string &name)
{
    NekDouble value = scalar.Attribute(element, name);
    if (value < 4.0 || value > static_cast<NekDouble>(INT_MAX) ||
        std::floor(value) != value || static_cast<int>(value) % 2)
    {
        Invalid(name + " must be an even integer in 4..INT_MAX");
    }
    return static_cast<int>(value);
}

struct Resolved
{
    std::shared_ptr<HITInflowConfiguration> configuration =
        std::make_shared<HITInflowConfiguration>();
    bool tabulated = false;
    HITInterpolation interpolation = HITInterpolation::LogLinear;
    HITLowWavenumber low = HITLowWavenumber::K4;
    HITHighWavenumber high = HITHighWavenumber::Gaussian;
    NekDouble length = 1.0, amplitude = 1.0;
    std::string table;
};

std::string ReadBytes(const std::string &filename)
{
    std::ifstream input(filename, std::ios::binary | std::ios::ate);
    if (!input)
    {
        Invalid("cannot read spectrum file " + filename +
                " (FILE is relative to the process working directory)");
    }
    auto count = input.tellg();
    if (count < 0 || count > static_cast<std::streamoff>(INT_MAX))
    {
        Invalid("spectrum file size exceeds the communicator's int count limit");
    }
    std::string bytes(static_cast<std::size_t>(count), '\0');
    input.seekg(0);
    if (!bytes.empty())
    {
        input.read(&bytes[0], static_cast<std::streamsize>(bytes.size()));
    }
    if (!input)
    {
        Invalid("failed while reading spectrum file " + filename);
    }
    return bytes;
}

void SetSpectrum(Resolved &resolved)
{
    if (resolved.tabulated)
    {
        std::istringstream input(resolved.table);
        input.imbue(std::locale::classic());
        resolved.configuration->spectrum = std::make_shared<HITSpectrum>(
            HITSpectrum::ReadTable(input, resolved.interpolation, resolved.low,
                                   resolved.high));
    }
    else
    {
        resolved.configuration->spectrum = std::make_shared<HITSpectrum>(
            HITSpectrum::VonKarman(resolved.length, resolved.amplitude));
    }
}

Resolved Resolve(const SessionReaderSharedPtr &session,
                 const std::string &custom,
                 const std::map<int, std::string> &boundaryConditions,
                 const std::set<int> &boundaryRegions)
{
    if (custom.empty())
    {
        Invalid("HITInflow boundary tags require NEKTAR/HITINFLOW");
    }
    TiXmlDocument document;
    document.Parse(custom.c_str());
    if (document.Error())
    {
        Invalid("cannot parse HITINFLOW XML");
    }
    auto block = document.FirstChildElement("HITINFLOW");
    Attributes(block, {"VERSION", "UNITS"});
    if (std::string(block->Attribute("VERSION")) != "1" ||
        std::string(block->Attribute("UNITS")) != "SI")
    {
        Invalid("require HITINFLOW VERSION=\"1\" UNITS=\"SI\"");
    }
    auto region = Children(block, {"REGION"}).at("REGION");
    Attributes(region, {"REF"});
    Resolved result;
    auto &configuration = *result.configuration;
    auto &generator = configuration.generator;
    configuration.regionId = RegionId(region->Attribute("REF"));
    if (!boundaryRegions.count(configuration.regionId))
    {
        Invalid("configured REF is absent from global BOUNDARYREGIONS");
    }
    if (boundaryConditions.size() != 1 ||
        !boundaryConditions.count(configuration.regionId))
    {
        Invalid("require exactly one global HITInflow boundary region matching "
                "HITINFLOW/REGION/@REF");
    }

    SolverInfo(session, "SolverType", "VelocityCorrectionScheme");
    SolverInfo(session, "EQTYPE", "UnsteadyNavierStokes");
    SolverInfo(session, "EvolutionOperator", "Nonlinear", "Nonlinear");
    SolverInfo(session, "Projection", "Continuous", "Continuous");
    SolverInfo(session, "Extrapolation", "Standard", "Standard");
    SolverInfo(session, "Driver", "Standard", "Standard");
    if (session->DefinesSolverInfo("HOMOGENEOUS") ||
        session->DefinesSolverInfo("ModeType"))
    {
        Invalid("Milestone 3 supports full 3D only; Homogeneous1D and Fourier "
                "plane support are scheduled for Milestone 4");
    }
    if (session->GetComm()->IsParallelInTime())
    {
        Invalid("parallel-in-time integration is unsupported");
    }
    if (!session->DefinesTimeIntScheme())
    {
        Invalid("require TIMEINTEGRATIONSCHEME with METHOD=IMEX and ORDER=2");
    }
    const auto &scheme = session->GetTimeIntScheme();
    if (scheme.method != "IMEX" || !scheme.variant.empty() || scheme.order != 2 ||
        !scheme.freeParams.empty())
    {
        Invalid("require IMEX order 2 with empty variant and no free parameters");
    }
    if (session->GetVariables() != std::vector<std::string>{"u", "v", "w", "p"})
    {
        Invalid("require exactly the ordered variables u,v,w,p");
    }
    auto root = session->GetDocument().FirstChildElement("NEKTAR");
    auto conditions = root->FirstChildElement("CONDITIONS");
    if (!conditions)
    {
        Invalid("missing CONDITIONS");
    }
    if (root->FirstChildElement("MOVEMENT") ||
        root->FirstChildElement("MAPPING"))
    {
        Invalid("moving geometry, ALE and coordinate mappings are unsupported");
    }
    auto geometry = root->FirstChildElement("GEOMETRY");
    if (geometry &&
        ((geometry->Attribute("DIM") &&
          std::string(geometry->Attribute("DIM")) != "3") ||
         (geometry->Attribute("SPACE") &&
          std::string(geometry->Attribute("SPACE")) != "3")))
    {
        Invalid("require full three-dimensional geometry");
    }
    bool body = false;
    for (auto forcing = root->FirstChildElement("FORCING"); forcing;
         forcing = forcing->NextSiblingElement("FORCING"))
    {
        for (auto force = forcing->FirstChildElement(); force;
             force = force->NextSiblingElement())
        {
            const char *type = force->Attribute("TYPE");
            if (force->ValueStr() != "FORCE" || !type || Upper(type) != "BODY")
            {
                Invalid("only ordinary Body momentum forcing is supported; "
                        "moving-frame, stochastic, and controller forcing "
                        "cannot define or change the HIT convection speed");
            }
            body = true;
        }
    }
    for (const auto &parameter :
         std::map<std::string, NekDouble>{{"FinTime", 0.0}, {"CFL", 0.0},
                                          {"CFLEnd", 0.0}, {"IO_CheckTime", 0.0},
                                          {"CFLGrowth", 1.0},
                                          {"TimeIncrementFactor", 1.0},
                                          {"Flowrate", 0.0}})
    {
        if (Parameter(session, parameter.first, parameter.second) !=
            parameter.second)
        {
            Invalid("require fixed-step/controller setting " + parameter.first +
                    "=" + std::to_string(parameter.second));
        }
    }
    generator.timeStep = Parameter(session, "TimeStep", 0.0);
    NekDouble numSteps = Parameter(session, "NumSteps", 0.0);
    if (!(generator.timeStep > 0.0) || numSteps < 1.0 || numSteps > INT_MAX ||
        std::floor(numSteps) != numSteps)
    {
        Invalid("require positive TimeStep and positive integral NumSteps <= "
                "INT_MAX");
    }

    ScalarExpression scalar(session, conditions);
    TiXmlDocument boundaryDocument;
    boundaryDocument.Parse(boundaryConditions.at(configuration.regionId).c_str());
    auto boundary = boundaryDocument.FirstChildElement("REGION");
    Attributes(boundary, {"REF"});
    std::map<std::string, const TiXmlElement *> bcs;
    for (auto bc = boundary->FirstChildElement(); bc;
         bc = bc->NextSiblingElement())
    {
        Attributes(bc, {"VAR", "USERDEFINEDTYPE", "VALUE"});
        Children(bc, {});
        std::string variable = bc->Attribute("VAR");
        if ((variable != "u" && variable != "v" && variable != "w" &&
             variable != "p") ||
            !bcs.emplace(variable, bc).second)
        {
            Invalid("HIT inlet must prescribe u,v,w,p exactly once each");
        }
        bool pressure = variable == "p";
        if (bc->ValueStr() != (pressure ? "N" : "D") ||
            std::string(bc->Attribute("USERDEFINEDTYPE")) !=
                (pressure ? "H" : "HITInflow"))
        {
            Invalid("inlet requires three D/HITInflow velocity conditions and "
                    "one N/H pressure condition");
        }
    }
    if (bcs.size() != 4)
    {
        Invalid("inlet requires all three HITInflow velocities and H pressure");
    }
    for (const auto &variable : {"u", "v", "w", "p"})
    {
        NekDouble value = scalar.Attribute(bcs.at(variable), "VALUE");
        if (std::string(variable) == "u")
        {
            if (!(value > 0.0))
            {
                Invalid("mean inlet u must be positive");
            }
            generator.meanSpeed = value;
        }
        else if (value != 0.0)
        {
            Invalid("inlet v,w and H pressure VALUE must be zero");
        }
    }

    auto children = Children(region, {"ORIGIN", "SYNTHESIS", "SPECTRUM", "CUTOFF"});
    auto origin = children.at("ORIGIN");
    auto synthesis = children.at("SYNTHESIS");
    auto spectrum = children.at("SPECTRUM");
    auto cutoff = children.at("CUTOFF");
    Attributes(origin, {"X", "Y", "Z"});
    Attributes(synthesis,
               {"LX", "LY", "LZ", "NX", "NY", "NZ", "SEED", "TIMEORIGIN",
                "RECORDDURATION", "REPEAT"});
    Attributes(cutoff, {"KMAX", "TAPERSTART"});
    for (auto element : {origin, synthesis, spectrum, cutoff})
    {
        Children(element, {});
    }
    const std::array<std::string, 3> coordinates{"X", "Y", "Z"};
    for (int direction = 0; direction < 3; ++direction)
    {
        generator.origin[direction] = scalar.Attribute(origin, coordinates[direction]);
        generator.lengths[direction] =
            scalar.Attribute(synthesis, "L" + coordinates[direction]);
        generator.counts[direction] =
            GridCount(scalar, synthesis, "N" + coordinates[direction]);
    }
    generator.seed = HITGenerator::ParseSeed(synthesis->Attribute("SEED"));
    generator.timeOrigin = scalar.Attribute(synthesis, "TIMEORIGIN");
    generator.recordDuration = scalar.Attribute(synthesis, "RECORDDURATION");
    std::string repeat = synthesis->Attribute("REPEAT");
    if (repeat != "Error" && repeat != "Allow")
    {
        Invalid("SYNTHESIS/@REPEAT must be Error or Allow");
    }
    generator.allowRepeat = repeat == "Allow";
    generator.kMax = scalar.Attribute(cutoff, "KMAX");
    generator.taperStart = scalar.Attribute(cutoff, "TAPERSTART");

    // TYPE and NORMALIZATION are checked before the type-specific attribute set.
    const char *typeAttribute = spectrum->Attribute("TYPE");
    const char *normalizationAttribute = spectrum->Attribute("NORMALIZATION");
    if (!typeAttribute || !normalizationAttribute)
    {
        Invalid("SPECTRUM requires TYPE and NORMALIZATION");
    }
    std::string type = typeAttribute;
    std::string normalization = normalizationAttribute;
    if ((type != "Tabulated" && type != "VonKarman") ||
        (normalization != "Preserve" && normalization != "TargetRMS"))
    {
        Invalid("unsupported SPECTRUM TYPE or NORMALIZATION");
    }
    generator.normalization = normalization == "TargetRMS"
                                  ? HITNormalization::TargetRMS
                                  : HITNormalization::Preserve;
    std::set<std::string> required{"TYPE", "NORMALIZATION"};
    if (generator.normalization == HITNormalization::TargetRMS)
    {
        required.insert("RMS");
    }
    std::string filename, interpolation, low, high;
    if (type == "Tabulated")
    {
        required.insert({"FILE", "INTERPOLATION", "LOW", "HIGH"});
        Attributes(spectrum, required);
        filename = spectrum->Attribute("FILE");
        interpolation = spectrum->Attribute("INTERPOLATION");
        low = spectrum->Attribute("LOW");
        high = spectrum->Attribute("HIGH");
        if ((interpolation != "Linear" && interpolation != "LogLinear") ||
            (low != "Zero" && low != "K4") ||
            (high != "Zero" && high != "Gaussian"))
        {
            Invalid("unsupported tabulated interpolation or low/high tail policy");
        }
        result.tabulated = true;
        result.interpolation = interpolation == "Linear"
                                   ? HITInterpolation::Linear
                                   : HITInterpolation::LogLinear;
        result.low = low == "Zero" ? HITLowWavenumber::Zero : HITLowWavenumber::K4;
        result.high = high == "Zero" ? HITHighWavenumber::Zero
                                     : HITHighWavenumber::Gaussian;
    }
    else
    {
        required.insert("LENGTH");
        if (generator.normalization == HITNormalization::Preserve)
        {
            required.insert("AMPLITUDE");
        }
        Attributes(spectrum, required);
        result.length = scalar.Attribute(spectrum, "LENGTH");
        result.amplitude = generator.normalization == HITNormalization::Preserve
                               ? scalar.Attribute(spectrum, "AMPLITUDE")
                               : 1.0;
    }
    if (generator.normalization == HITNormalization::TargetRMS)
    {
        generator.targetRMS = scalar.Attribute(spectrum, "RMS");
        if (generator.targetRMS < 0.0)
        {
            Invalid("SPECTRUM/@RMS must be nonnegative");
        }
    }
    if (result.tabulated)
    {
        // Resolve runs on communicator rank zero only, before any dependent
        // broadcast. The digest includes headers, comments and line endings.
        result.table = ReadBytes(filename);
        configuration.spectrumDigest = HITInflowConfig::SHA256(result.table);
    }
    SetSpectrum(result);
    std::size_t modeBytes =
        HITGenerator::EstimateModeBytes(generator, *configuration.spectrum);
    long double duration = static_cast<long double>(numSteps) * generator.timeStep;
    long double tolerance = 64.0L * std::numeric_limits<NekDouble>::epsilon() *
                            std::max(1.0L, std::abs(static_cast<long double>(
                                               generator.recordDuration)));
    if (duration > static_cast<long double>(generator.recordDuration) + tolerance)
    {
        Invalid("NumSteps*TimeStep exceeds SYNTHESIS/@RECORDDURATION; job start "
                "and end times must also lie in the declared record interval");
    }

    std::ostringstream description;
    description.imbue(std::locale::classic());
    description << std::setprecision(std::numeric_limits<NekDouble>::max_digits10)
                << "HITInflow VERSION=1 UNITS=SI algorithm=" << HITGenerator::Algorithm()
                << " expansion=Full3D regionREF=" << configuration.regionId
                << " seed=" << generator.seed;
    for (int direction = 0; direction < 3; ++direction)
    {
        description << " L" << coordinates[direction] << '=' << generator.lengths[direction]
                    << " N" << coordinates[direction] << '=' << generator.counts[direction]
                    << " ORIGIN." << coordinates[direction] << '=' << generator.origin[direction];
    }
    description << " meanU=" << generator.meanSpeed << " meanV=0 meanW=0"
                << " timeOrigin=" << generator.timeOrigin
                << " recordDuration=" << generator.recordDuration
                << " TimeStep=" << generator.timeStep
                << " NumSteps=" << static_cast<int>(numSteps)
                << " repeat=" << repeat << " KMAX=" << generator.kMax
                << " TAPERSTART=" << generator.taperStart
                << " spectrum=" << type << " normalization=" << normalization;
    if (generator.normalization == HITNormalization::TargetRMS)
    {
        description << " targetRMS=" << generator.targetRMS;
    }
    if (result.tabulated)
    {
        description << " FILE=" << std::quoted(filename)
                    << " interpolation=" << interpolation << " low=" << low
                    << " high=" << high << " SHA256=" << configuration.spectrumDigest;
    }
    else
    {
        description << " length=" << result.length << " amplitude=" << result.amplitude;
    }
    description << " estimatedModeBytes=" << modeBytes
                << " externalMomentumSource=" << (body ? "Body" : "None");
    configuration.description = description.str();
    return result;
}

void BroadcastResolved(const CommSharedPtr &comm, Resolved &resolved)
{
    auto &configuration = *resolved.configuration;
    auto &generator = configuration.generator;
    std::vector<int> integers(10);
    std::vector<NekDouble> reals(15);
    std::string seed;
    if (comm->GetRank() == 0)
    {
        integers = {configuration.regionId, generator.counts[0], generator.counts[1],
                    generator.counts[2], generator.allowRepeat ? 1 : 0,
                    static_cast<int>(generator.normalization), resolved.tabulated ? 1 : 0,
                    static_cast<int>(resolved.interpolation), static_cast<int>(resolved.low),
                    static_cast<int>(resolved.high)};
        reals = {generator.lengths[0], generator.lengths[1], generator.lengths[2],
                 generator.origin[0], generator.origin[1], generator.origin[2],
                 generator.meanSpeed, generator.timeOrigin, generator.recordDuration,
                 generator.timeStep, generator.kMax, generator.taperStart,
                 generator.targetRMS, resolved.length, resolved.amplitude};
        seed = std::to_string(generator.seed);
    }
    comm->Bcast(integers, 0);
    comm->Bcast(reals, 0);
    BroadcastString(comm, seed, 0);
    BroadcastString(comm, resolved.table, 0);
    BroadcastString(comm, configuration.spectrumDigest, 0);
    BroadcastString(comm, configuration.description, 0);
    configuration.regionId = integers[0];
    generator.counts = {integers[1], integers[2], integers[3]};
    generator.allowRepeat = integers[4] != 0;
    generator.normalization = static_cast<HITNormalization>(integers[5]);
    resolved.tabulated = integers[6] != 0;
    resolved.interpolation = static_cast<HITInterpolation>(integers[7]);
    resolved.low = static_cast<HITLowWavenumber>(integers[8]);
    resolved.high = static_cast<HITHighWavenumber>(integers[9]);
    generator.lengths = {reals[0], reals[1], reals[2]};
    generator.origin = {reals[3], reals[4], reals[5]};
    generator.meanSpeed = reals[6];
    generator.timeOrigin = reals[7];
    generator.recordDuration = reals[8];
    generator.timeStep = reals[9];
    generator.kMax = reals[10];
    generator.taperStart = reals[11];
    generator.targetRMS = reals[12];
    resolved.length = reals[13];
    resolved.amplitude = reals[14];
    generator.seed = HITGenerator::ParseSeed(seed);
    if (comm->GetRank() != 0)
    {
        SetSpectrum(resolved);
    }
}

std::uint32_t RotateRight(std::uint32_t value, int count)
{
    return (value >> count) | (value << (32 - count));
}
} // namespace

std::shared_ptr<const HITInflowConfiguration> HITInflowConfig::Load(
    const LibUtilities::SessionReaderSharedPtr &session)
{
    if (!session)
    {
        Invalid("null session");
    }
    auto comm = session->GetComm();
    LocalXML local;
    std::string error;
    try
    {
        local = CollectXML(session);
    }
    catch (...)
    {
        error = CurrentException();
    }
    // Even an invalid HIT block must take the collective error path. Avoid
    // imposing new checks on sessions which do not use this boundary feature.
    int present = local.present || !error.empty() ? 1 : 0;
    comm->AllReduce(present, LibUtilities::ReduceMax);
    if (!present)
    {
        return {};
    }
    AgreeError(comm, error);

    std::string custom = local.custom;
    BroadcastString(comm, custom, 0);
    error.clear();
    if (custom != local.custom)
    {
        error = "HITInflow: HITINFLOW block differs between session ranks";
    }
    AgreeError(comm, error);

    // XML partition files omit boundary regions not owned by that rank. Rank
    // zero need not own the inlet: merge only the small tagged BC descriptions,
    // retaining global REF rather than a rank-local boundary-array index.
    std::map<int, std::string> boundaryConditions;
    std::set<int> boundaryRegions;
    for (int rank = 0; rank < comm->GetSize(); ++rank)
    {
        std::string boundaries = comm->GetRank() == rank ? local.boundaries : "";
        BroadcastString(comm, boundaries, rank);
        if (comm->GetRank() == 0 && error.empty())
        {
            try
            {
                TiXmlDocument document;
                document.Parse(boundaries.c_str());
                auto collection = document.FirstChildElement("BOUNDARIES");
                if (document.Error() || !collection)
                {
                    Invalid("cannot collect partitioned boundary metadata");
                }
                for (auto region = collection->FirstChildElement(); region;
                     region = region->NextSiblingElement())
                {
                    if (region->ValueStr() == "B")
                    {
                        boundaryRegions.insert(RegionId(region->Attribute("ID")));
                    }
                    else
                    {
                        int ref = RegionId(region->Attribute("REF"));
                        std::string xml = XML(region);
                        auto inserted = boundaryConditions.emplace(ref, xml);
                        if (!inserted.second && inserted.first->second != xml)
                        {
                            Invalid("HIT boundary conditions differ between ranks "
                                    "for region REF=" + std::to_string(ref));
                        }
                    }
                }
            }
            catch (...)
            {
                error = CurrentException();
            }
        }
    }
    AgreeError(comm, error);

    Resolved resolved;
    error.clear();
    if (comm->GetRank() == 0)
    {
        try
        {
            resolved = Resolve(session, custom, boundaryConditions, boundaryRegions);
        }
        catch (...)
        {
            error = CurrentException();
        }
    }
    AgreeError(comm, error);
    try
    {
        BroadcastResolved(comm, resolved);
    }
    catch (...)
    {
        error = CurrentException();
    }
    AgreeError(comm, error);
    return resolved.configuration;
}

std::string HITInflowConfig::SHA256(const std::string &bytes)
{
    // FIPS 180-4 SHA-256. Work on octets and explicit big-endian words so the
    // digest covers the exact file bytes independently of host byte order.
    static constexpr std::array<std::uint32_t, 64> constants{
        0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u, 0x3956c25bu,
        0x59f111f1u, 0x923f82a4u, 0xab1c5ed5u, 0xd807aa98u, 0x12835b01u,
        0x243185beu, 0x550c7dc3u, 0x72be5d74u, 0x80deb1feu, 0x9bdc06a7u,
        0xc19bf174u, 0xe49b69c1u, 0xefbe4786u, 0x0fc19dc6u, 0x240ca1ccu,
        0x2de92c6fu, 0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau, 0x983e5152u,
        0xa831c66du, 0xb00327c8u, 0xbf597fc7u, 0xc6e00bf3u, 0xd5a79147u,
        0x06ca6351u, 0x14292967u, 0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu,
        0x53380d13u, 0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u,
        0xa2bfe8a1u, 0xa81a664bu, 0xc24b8b70u, 0xc76c51a3u, 0xd192e819u,
        0xd6990624u, 0xf40e3585u, 0x106aa070u, 0x19a4c116u, 0x1e376c08u,
        0x2748774cu, 0x34b0bcb5u, 0x391c0cb3u, 0x4ed8aa4au, 0x5b9cca4fu,
        0x682e6ff3u, 0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u,
        0x90befffau, 0xa4506cebu, 0xbef9a3f7u, 0xc67178f2u};
    std::array<std::uint32_t, 8> state{
        0x6a09e667u, 0xbb67ae85u, 0x3c6ef372u, 0xa54ff53au,
        0x510e527fu, 0x9b05688cu, 0x1f83d9abu, 0x5be0cd19u};
    if (bytes.size() > std::numeric_limits<std::uint64_t>::max() / 8)
    {
        throw std::overflow_error("HITInflow: SHA-256 input bit length overflows");
    }
    const auto bitLength = static_cast<std::uint64_t>(bytes.size()) * 8;
    auto compress = [&](const std::array<unsigned char, 64> &block)
    {
        std::array<std::uint32_t, 64> words{};
        for (int word = 0; word < 16; ++word)
        {
            for (int octet = 0; octet < 4; ++octet)
            {
                words[word] = (words[word] << 8) | block[4 * word + octet];
            }
        }
        for (int word = 16; word < 64; ++word)
        {
            auto x = words[word - 15], y = words[word - 2];
            auto sigma0 = RotateRight(x, 7) ^ RotateRight(x, 18) ^ (x >> 3);
            auto sigma1 = RotateRight(y, 17) ^ RotateRight(y, 19) ^ (y >> 10);
            words[word] = words[word - 16] + sigma0 + words[word - 7] + sigma1;
        }
        auto a = state[0], b = state[1], c = state[2], d = state[3];
        auto e = state[4], f = state[5], g = state[6], h = state[7];
        for (int word = 0; word < 64; ++word)
        {
            auto sum1 = RotateRight(e, 6) ^ RotateRight(e, 11) ^ RotateRight(e, 25);
            auto choice = (e & f) ^ (~e & g);
            auto temporary1 = h + sum1 + choice + constants[word] + words[word];
            auto sum0 = RotateRight(a, 2) ^ RotateRight(a, 13) ^ RotateRight(a, 22);
            auto majority = (a & b) ^ (a & c) ^ (b & c);
            auto temporary2 = sum0 + majority;
            h = g;
            g = f;
            f = e;
            e = d + temporary1;
            d = c;
            c = b;
            b = a;
            a = temporary1 + temporary2;
        }
        state[0] += a;
        state[1] += b;
        state[2] += c;
        state[3] += d;
        state[4] += e;
        state[5] += f;
        state[6] += g;
        state[7] += h;
    };
    std::array<unsigned char, 64> block{};
    std::size_t offset = 0;
    while (bytes.size() - offset >= block.size())
    {
        for (std::size_t octet = 0; octet < block.size(); ++octet)
        {
            block[octet] = static_cast<unsigned char>(bytes[offset + octet]);
        }
        compress(block);
        offset += block.size();
    }
    std::size_t remaining = bytes.size() - offset;
    block.fill(0);
    for (std::size_t octet = 0; octet < remaining; ++octet)
    {
        block[octet] = static_cast<unsigned char>(bytes[offset + octet]);
    }
    block[remaining] = 0x80;
    if (remaining >= 56)
    {
        compress(block);
        block.fill(0);
    }
    for (int octet = 0; octet < 8; ++octet)
    {
        block[63 - octet] = static_cast<unsigned char>(bitLength >> (8 * octet));
    }
    compress(block);
    std::ostringstream digest;
    digest.imbue(std::locale::classic());
    digest << std::hex << std::setfill('0');
    for (auto word : state)
    {
        digest << std::setw(8) << word;
    }
    return digest.str();
}
} // namespace Nektar
