// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/WJNativeScalar.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeFailures.h>
#include <SpatialDomains/PointGeom.h>
#include <Eigen/Cholesky>
#include <Eigen/QR>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>

namespace Nektar::WJNative
{
namespace
{
using Triplet = Eigen::Triplet<double>;
constexpr int EdgeModes[4][4] = {
    {0, 4, 8, 12}, {1, 5, 9, 13}, {0, 1, 2, 3}, {4, 5, 6, 7}};

void Check(bool condition, const char *message)
{
    if (!condition)
    {
        throw std::runtime_error(message);
    }
}

double Clock()
{
    return std::chrono::duration<double>(
               std::chrono::steady_clock::now().time_since_epoch()).count();
}

void CheckSparse(const Sparse &matrix, int size)
{
    Check(matrix.rows() == size && matrix.cols() == size,
          "native scalar matrix dimension differs");
    for (int column = 0; column < matrix.outerSize(); ++column)
    {
        for (Sparse::InnerIterator it(matrix, column); it; ++it)
        {
            Check(std::isfinite(it.value()), "nonfinite native scalar matrix");
        }
    }
}

double BernsteinLower(const Eigen::Matrix4d &coefficients,int remaining,double reserve)
{
    const double lower=coefficients.minCoeff()-reserve;
    if(lower>0.||remaining==0)return lower;
    const Eigen::Matrix4d left=(Eigen::Matrix4d()<<
        1.,0.,0.,0.,
        .5,.5,0.,0.,
        .25,.5,.25,0.,
        .125,.375,.375,.125).finished();
    const Eigen::Matrix4d right=left.reverse();
    const Eigen::Matrix4d first=left*coefficients,second=right*coefficients;
    return std::min({BernsteinLower(first*left.transpose(),remaining-1,reserve),
        BernsteinLower(first*right.transpose(),remaining-1,reserve),
        BernsteinLower(second*left.transpose(),remaining-1,reserve),
        BernsteinLower(second*right.transpose(),remaining-1,reserve)});
}
}

Vec Scalar::Evaluate(const NativeRule &rule,const Vec &state,int derivative) const
{
    Check(state.size()==m_size,"scalar evaluation state size differs");
    const Cell *cell=Find(rule);
    Check(cell!=nullptr,"scalar evaluation has no specified native owner");
    return Basis(rule,*cell,derivative,false)*Local(state,cell->ids);
}

std::array<double,2> Scalar::PolynomialBounds(const Vec&state,int subdivisionDepth)const
{
    Check(state.size()==m_size&&state.allFinite(),"scalar polynomial-bound state differs");
    Check(subdivisionDepth>=0&&subdivisionDepth<=6,"scalar Bernstein depth out of range");
    Eigen::Matrix4d transform;
    transform<<1,0,0,0, 2./3,1./3,1./3,-2./3,
               1./3,2./3,1./3,2./3, 0,1,0,0;
    std::array<double,2>bounds{std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity()};
    for(const auto&cell:m_cells){Eigen::Matrix4d local;
        for(int j=0;j<16;++j)local(j%4,j/4)=state[cell.ids[j]]*cell.signs[j];
        Eigen::Matrix4d bernstein=transform*local*transform.transpose();
        double lower=bernstein.minCoeff(),upper=bernstein.maxCoeff();
        if(subdivisionDepth>0)
        {
            // Keep a reserve for the basis transformation and the convex
            // combinations in the dyadic de Casteljau subdivisions.
            const double reserve=128.*(1+subdivisionDepth)*std::numeric_limits<double>::epsilon()*
                local.cwiseAbs().maxCoeff();
            lower=BernsteinLower(bernstein,subdivisionDepth,reserve);
            upper+=reserve;
        }
        bounds[0]=std::min(bounds[0],lower);bounds[1]=std::max(bounds[1],upper);
    }
    return bounds;
}

Scalar::Scalar(const MultiRegions::ContFieldSharedPtr &field, const Flow &flow,
               Function residualScale, std::vector<int> essentialSides)
    : m_field(field), m_flow(flow), m_referenceScale(residualScale),
      m_explicitEssential(!essentialSides.empty())
{
    Check(bool(field) && field->GetComm()->GetSize() == 1,
          "native scalar requires a serial ContField");
    Check(bool(residualScale), "native scalar needs a positive reference scale");
    Check(flow.Mapping().curve == 0,
          "native scalar requires affine chart rectangles");
    Check(field->GetExpSize() > 0, "empty native scalar domain");
    for (int side : essentialSides)
    {
        Check(side >= 0 && side < 4, "invalid scalar essential side");
        m_essentialSides[side] = true;
    }
    auto map = field->GetLocalToGlobalMap();
    m_size = map->GetNumGlobalCoeffs();
    Check(m_size > 0, "empty native scalar coefficient space");
    std::vector<bool> used(m_size, false), essential(m_size, false);
    std::unordered_map<int, const Element *> geometry;
    for (const auto &element : flow.Elements())
    {
        Check(geometry.emplace(element.geometry, &element).second,
              "duplicate native flow geometry ID");
    }
    double xl = std::numeric_limits<double>::infinity(), rl = xl;
    double xr = -xl, rr = -xl, area = 0;
    for (int e = 0; e < field->GetExpSize(); ++e)
    {
        auto expansion = field->GetExp(e);
        const int id = expansion->GetGeom()->GetGlobalID();
        auto found = geometry.find(id);
        Check(found != geometry.end(), "scalar geometry is absent from flow");
        Check(m_geometryToCell.emplace(id, int(m_cells.size())).second,
              "duplicate native scalar geometry ID");
        const auto &element = *found->second;
        Check(expansion->GetNumBases() == 2 && expansion->GetNcoeffs() == 16,
              "native scalar requires tensor Q3");
        Check(expansion->GetGeom()->GetNumVerts() == 4,
              "native scalar requires quadrilateral elements");
        for (int d = 0; d < 2; ++d)
        {
            Check(expansion->GetBasisNumModes(d) == 4 &&
                      expansion->GetBasis(d)->GetBasisType() ==
                          LibUtilities::eModified_A &&
                      element.velocity->GetBasis(d)->GetBasisType() ==
                          LibUtilities::eModified_A,
                  "scalar/rule basis must be native Q3 Modified_A");
        }
        for (int vertex = 0; vertex < 4; ++vertex)
        {
            Check(expansion->GetGeom()->GetVid(vertex) ==
                      element.velocity->GetGeom()->GetVid(vertex),
                  "scalar/flow element orientation differs");
        }
        // Bounds come from the scalar field's own native geometry. The flow
        // element is used only after geometry and local orientation agree.
        double exl = std::numeric_limits<double>::infinity(), erl = exl;
        double exr = -exl, err = -exl;
        for (int vertex = 0; vertex < 4; ++vertex)
        {
            double x = 0, z = 0, dummy = 0;
            expansion->GetGeom()->GetVertex(vertex)->GetCoords(x, z, dummy);
            Check(std::isfinite(x) && std::isfinite(z),
                  "nonfinite native scalar geometry coordinate");
            exl = std::min(exl, x);
            exr = std::max(exr, x);
            const double r = z - flow.Mapping().Ground(x);
            erl = std::min(erl, r);
            err = std::max(err, r);
        }
        const double geometryTolerance = 1e-11 *
            std::max({1., std::abs(exl), std::abs(exr), std::abs(erl),
                      std::abs(err)});
        Check(std::abs(exl - element.xl) <= geometryTolerance &&
                  std::abs(exr - element.xr) <= geometryTolerance &&
                  std::abs(erl - element.rl) <= geometryTolerance &&
                  std::abs(err - element.rr) <= geometryTolerance,
              "scalar/flow native geometry bounds differ");
        Check(exr > exl && err > erl, "nonpositive scalar element area");
        xl = std::min(xl, exl); xr = std::max(xr, exr);
        rl = std::min(rl, erl); rr = std::max(rr, err);
        area += (exr - exl) * (err - erl);
        Cell cell;
        cell.flow = &element;
        for (int j = 0; j < 16; ++j)
        {
            const int local = field->GetCoeff_Offset(e) + j;
            const int global = map->GetLocalToGlobalMap(local);
            const double sign = map->GetLocalToGlobalSign(local);
            Check(global >= 0 && global < m_size &&
                      std::abs(std::abs(sign) - 1.) < 1e-14,
                  "invalid signed native scalar map");
            cell.ids.push_back(global);
            cell.signs.push_back(sign);
            used[global] = true;
        }
        m_cells.push_back(std::move(cell));
    }
    Check(std::all_of(used.begin(), used.end(), [](bool value) { return value; }),
          "incomplete native scalar map");
    Check(std::abs(area - (xr - xl) * (rr - rl)) <= 1e-10 * area,
          "scalar geometry subset is not a complete rectangle");
    RuleFactory factory(flow);
    const double boundaryTolerance = 1e-11 *
        std::max({1., std::abs(xl), std::abs(xr), std::abs(rl), std::abs(rr)});
    for (int e = 0; e < int(m_cells.size()); ++e)
    {
        const auto &cell = m_cells[e];
        const auto &element = *cell.flow;
        const std::array<bool, 4> boundary = {
            std::abs(element.xl - xl) <= boundaryTolerance,
            std::abs(element.xr - xr) <= boundaryTolerance,
            std::abs(element.rl - rl) <= boundaryTolerance,
            std::abs(element.rr - rr) <= boundaryTolerance};
        for (int side = 0; side < 4; ++side)
        {
            if (!boundary[side])
            {
                continue;
            }
            m_faces.push_back({e, side, factory.Face(element, side, 16)});
            if (m_explicitEssential && m_essentialSides[side])
            {
                for (int mode : EdgeModes[side])
                {
                    essential[cell.ids[mode]] = true;
                }
            }
        }
    }
    if (!m_explicitEssential)
    {
        const int nativeDirichlet = map->GetNumGlobalDirBndCoeffs();
        Check(nativeDirichlet >= 0 && nativeDirichlet <= m_size,
              "invalid native scalar Dirichlet count");
        for (int i = 0; i < nativeDirichlet; ++i)
        {
            essential[i] = true;
        }
    }
    m_fullToFree.assign(m_size, -1);
    for (int i = 0; i < m_size; ++i)
    {
        if (essential[i])
        {
            m_essentialIds.push_back(i);
        }
        else
        {
            m_fullToFree[i] = int(m_freeIds.size());
            m_freeIds.push_back(i);
        }
    }
    Check(FreeSize() > 0, "empty homogeneous native scalar space");

    std::vector<Triplet> mass, weighted;
    for (const auto &pointer : factory.Static())
    {
        Check(bool(pointer), "null scalar mass rule");
        const auto &rule = *pointer;
        const Cell *cell = Find(rule);
        if (!cell)
        {
            continue;
        }
        Mat b = Basis(rule, *cell, -1);
        Vec weights(rule.w.size());
        for (int q = 0; q < weights.size(); ++q)
        {
            const double scale = residualScale(rule.x[q], rule.z[q]);
            Check(std::isfinite(scale) && scale > 0,
                  "nonpositive/nonfinite scalar residual reference scale");
            weights[q] = rule.w[q] * scale * scale;
            Check(std::isfinite(weights[q]) && weights[q] > 0,
                  "invalid weighted scalar Riesz quadrature");
        }
        Mat localMass = b.transpose() * rule.w.asDiagonal() * b;
        Mat localWeighted = b.transpose() * weights.asDiagonal() * b;
        for (int i = 0; i < 16; ++i)
        {
            for (int j = 0; j < 16; ++j)
            {
                mass.emplace_back(cell->ids[i], cell->ids[j], localMass(i, j));
                weighted.emplace_back(cell->ids[i], cell->ids[j],
                                      localWeighted(i, j));
            }
        }
        m_volume += rule.w.sum();
    }
    Check(m_volume > 0 && (flow.Mapping().normalMap ||
                          std::abs(m_volume - area) <= 1e-10 * area),
          "native scalar mass rules do not cover its domain");
    m_mass.resize(m_size, m_size);
    m_mass.setFromTriplets(mass.begin(), mass.end());
    Sparse weightedMass(m_size, m_size);
    weightedMass.setFromTriplets(weighted.begin(), weighted.end());
    const Sparse freeWeighted = Gather(weightedMass);
    m_riesz.compute(freeWeighted);
    Check(m_riesz.info() == Eigen::Success &&
              (m_riesz.vectorD().array() > 0).all(),
          "native scalar reference Riesz factor failed");
    m_scale.resize(FreeSize());
    for (int i = 0; i < FreeSize(); ++i)
    {
        Check(freeWeighted.coeff(i, i) > 0,
              "nonpositive scalar reference Riesz diagonal");
        m_scale[i] = 1. / std::sqrt(freeWeighted.coeff(i, i));
    }
}

const Scalar::Cell *Scalar::Find(const NativeRule &rule) const
{
    auto found = m_geometryToCell.find(rule.geometry);
    return found == m_geometryToCell.end() ? nullptr : &m_cells[found->second];
}

Mat Scalar::Basis(const NativeRule &rule, const Cell &cell, int derivative,
                  bool integration) const
{
    Check(rule.geometry == cell.flow->geometry && rule.u == cell.flow->u &&
              rule.us.size() == 16 && rule.x.size() == rule.w.size() &&
              rule.z.size() == rule.w.size() && rule.w.size() > 0 &&
              rule.x.allFinite() && rule.z.allFinite() && rule.w.allFinite() &&
              (integration?(rule.w.array()>0).all():(rule.w.array()>=0).all()),
          "native scalar rule geometry/map/weights differ");
    Mat basis = derivative < 0 ? rule.B : derivative == 0 ? rule.Dx : rule.Dz;
    Check(basis.rows() == rule.w.size() && basis.cols() == 16 && basis.allFinite(),
          "native scalar rule basis dimension/data differ");
    for (int j = 0; j < 16; ++j)
    {
        Check(std::abs(std::abs(rule.us[j]) - 1.) < 1e-14 &&
                  rule.us[j] == cell.flow->us[j],
              "native scalar rule sign differs from canonical velocity");
        basis.col(j) *= cell.signs[j] / rule.us[j];
    }
    return basis;
}

Vec Scalar::Local(const Vec &full, const std::vector<int> &ids) const
{
    Vec local(ids.size());
    for (int j = 0; j < int(ids.size()); ++j)
    {
        local[j] = full[ids[j]];
    }
    return local;
}

void Scalar::Add(Vec &full, const Vec &local, const Cell &cell) const
{
    Check(local.size() == int(cell.ids.size()), "native scalar local load size");
    for (int i = 0; i < int(cell.ids.size()); ++i)
    {
        full[cell.ids[i]] += local[i];
    }
}

Sparse Scalar::Gather(const Sparse &full) const
{
    CheckSparse(full, Size());
    std::vector<Triplet> triplets;
    for (int column = 0; column < full.outerSize(); ++column)
    {
        for (Sparse::InnerIterator it(full, column); it; ++it)
        {
            const int row = m_fullToFree[it.row()], col = m_fullToFree[it.col()];
            if (row >= 0 && col >= 0)
            {
                triplets.emplace_back(row, col, it.value());
            }
        }
    }
    Sparse result(FreeSize(), FreeSize());
    result.setFromTriplets(triplets.begin(), triplets.end());
    return result;
}

Vec Scalar::Project(const Function &function) const
{
    Check(bool(function), "missing scalar projection function");
    Array<OneD, NekDouble> x(m_field->GetTotPoints()), z(x.size()), dummy(x.size()),
        values(x.size()), coefficients(m_field->GetNcoeffs());
    m_field->GetCoords(x, z, dummy);
    for (int q = 0; q < int(x.size()); ++q)
    {
        values[q] = function(x[q], z[q]);
        Check(std::isfinite(values[q]), "nonfinite scalar projection value");
    }
    m_field->FwdTransBndConstrained(values, coefficients);
    if (!m_explicitEssential)
    {
        for (int b = 0; b < int(m_field->GetBndConditions().size()); ++b)
        {
            if (m_field->GetBndConditions()[b]->GetBoundaryConditionType() !=
                SpatialDomains::eDirichlet)
            {
                continue;
            }
            auto edge = m_field->GetBndCondExpansions()[b];
            Array<OneD, NekDouble> ex(edge->GetTotPoints()), ez(ex.size()),
                ed(ex.size()), ev(ex.size());
            edge->GetCoords(ex, ez, ed);
            for (int q = 0; q < int(ex.size()); ++q)
            {
                ev[q] = function(ex[q], ez[q]);
                Check(std::isfinite(ev[q]), "nonfinite scalar boundary projection");
            }
            edge->FwdTransBndConstrained(ev, edge->UpdateCoeffs());
        }
        m_field->ImposeDirichletConditions(coefficients);
    }
    m_field->UpdateCoeffs() = coefficients;
    m_field->AvgAssemble();
    m_field->GlobalToLocal();
    auto map = m_field->GetLocalToGlobalMap();
    Array<OneD, NekDouble> global(Size());
    map->LocalToGlobal(m_field->GetCoeffs(), global, true);
    Vec result(Size());
    for (int i = 0; i < Size(); ++i)
    {
        result[i] = global[i];
    }
    if (m_explicitEssential)
    {
        RuleFactory factory(m_flow);
        const Vec lift = EssentialLift(
            [&function](double px, double pz, const Vector &) {
                return function(px, pz);
            }, factory);
        for (int id : m_essentialIds)
        {
            result[id] = lift[id];
            global[id] = lift[id];
        }
        map->GlobalToLocal(global, m_field->UpdateCoeffs());
    }
    Check(result.allFinite(), "nonfinite projected scalar coefficients");
    return result;
}

Vec Scalar::EssentialLift(const ScalarBoundaryFunction &function,
                          const RuleFactory &factory) const
{
    Check(m_explicitEssential && bool(function),
          "EssentialLift needs explicit scalar boundary sides and a function");
    Vec lift = Vec::Zero(Size());
    std::vector<bool> assigned(Size(), false);
    auto assign = [&](int id, double value, double scale) {
        Check(std::isfinite(value), "nonfinite essential scalar coefficient");
        Check(std::isfinite(scale) && scale > 0, "invalid scalar corner reference scale");
        if (assigned[id])
        {
            const double tolerance = 512 * std::numeric_limits<double>::epsilon() *
                std::max({scale, std::abs(value), std::abs(lift[id])});
            if (std::abs(lift[id] - value) > tolerance)
            {
                throw std::runtime_error(
                    "incompatible scalar essential corner/edge coefficient " +
                    std::to_string(id));
            }
        }
        else
        {
            assigned[id] = true;
            lift[id] = value;
        }
    };
    for (const auto &face : m_faces)
    {
        if (!m_essentialSides[face.side])
        {
            continue;
        }
        const Cell &cell = m_cells[face.cell];
        const Element &element = *cell.flow;
        const NativeRule rule = factory.Face(element, face.side, 16);
        const Mat basis = Basis(rule, cell, -1);
        const int *modes = EdgeModes[face.side];
        std::array<double, 2> endpoints{};
        for (int endpoint = 0; endpoint < 2; ++endpoint)
        {
            const double x = face.side < 2 ?
                (face.side == 0 ? element.xl : element.xr) :
                (endpoint == 0 ? element.xl : element.xr);
            const double r = face.side < 2 ?
                (endpoint == 0 ? element.rl : element.rr) :
                (face.side == 2 ? element.rl : element.rr);
            endpoints[endpoint] = function(x, r + m_flow.Mapping().Ground(x),
                                           factory.FaceNormal(element,face.side,x,r)) /
                                  cell.signs[modes[endpoint]];
            assign(cell.ids[modes[endpoint]], endpoints[endpoint],
                   m_referenceScale(x, r + m_flow.Mapping().Ground(x)));
        }
        Vec values(rule.w.size());
        for (int q = 0; q < values.size(); ++q)
        {
            values[q] = function(rule.x[q], rule.z[q], rule.Outward(q));
            Check(std::isfinite(values[q]), "nonfinite essential scalar edge value");
        }
        values -= basis.col(modes[0]) * endpoints[0] +
                  basis.col(modes[1]) * endpoints[1];
        Mat interior(rule.w.size(), 2);
        interior.col(0) = basis.col(modes[2]);
        interior.col(1) = basis.col(modes[3]);
        const Mat mass = interior.transpose() * rule.w.asDiagonal() * interior;
        const Vec rhs = interior.transpose() * rule.w.asDiagonal() * values;
        Eigen::LDLT<Mat> factor(mass);
        Check(factor.info() == Eigen::Success && factor.isPositive(),
              "essential scalar edge interior mass failed");
        const Vec coefficients = factor.solve(rhs);
        Check(factor.info() == Eigen::Success && coefficients.allFinite(),
              "essential scalar edge projection failed");
        const double scale = m_referenceScale(rule.x.mean(), rule.z.mean());
        assign(cell.ids[modes[2]], coefficients[0], scale);
        assign(cell.ids[modes[3]], coefficients[1], scale);
    }
    for (int id : m_essentialIds)
    {
        Check(assigned[id], "incomplete scalar essential lift");
    }
    return lift;
}

Vec Scalar::Lift(const Vec &fullDirValues) const
{
    Check(fullDirValues.size() == Size() && fullDirValues.allFinite(),
          "native scalar lift dimension/data differ");
    Vec lift = Vec::Zero(Size());
    for (int id : m_essentialIds)
    {
        lift[id] = fullDirValues[id];
    }
    return lift;
}

Vec Scalar::Restrict(const Vec &full) const
{
    Check(full.size() == Size() && full.allFinite(),
          "native scalar full coefficient/load dimension/data differ");
    return Local(full, m_freeIds);
}

Vec Scalar::Expand(const Vec &free, const Vec &lift) const
{
    Check(free.size() == FreeSize() && free.allFinite(),
          "native scalar free coefficient dimension/data differ");
    Vec result = Lift(lift);
    for (int i = 0; i < FreeSize(); ++i)
    {
        result[m_freeIds[i]] = free[i];
    }
    return result;
}

Sparse Scalar::System(double alpha, const Sparse &K) const
{
    Check(std::isfinite(alpha) && alpha >= 0, "invalid scalar mass coefficient");
    CheckSparse(K, Size());
    Sparse full = alpha * m_mass + K;
    return Gather(full);
}

Vec Scalar::LiftedRhs(double alpha, const Sparse &K, const Vec &fullLoad,
                      const Vec &lift) const
{
    Check(std::isfinite(alpha) && alpha >= 0 && fullLoad.size() == Size() &&
              fullLoad.allFinite(), "invalid scalar lifted RHS data");
    CheckSparse(K, Size());
    const Vec essential = Lift(lift);
    return Restrict(fullLoad - alpha * (m_mass * essential) - K * essential);
}

Sparse Scalar::Diffusion(const Rules &rules, const std::vector<Vec> &diffusion,
                         const std::vector<Vec> &reaction) const
{
    Check(diffusion.size() == rules.size() &&
              (reaction.empty() || reaction.size() == rules.size()),
          "scalar diffusion/reaction rule count differs");
    std::vector<Triplet> triplets;
    for (int e = 0; e < int(rules.size()); ++e)
    {
        Check(bool(rules[e]), "null scalar diffusion rule");
        const auto &rule = *rules[e];
        const Cell *cell = Find(rule);
        if (!cell)
        {
            continue;
        }
        Check(diffusion[e].size() == rule.w.size() && diffusion[e].allFinite() &&
                  (diffusion[e].array() > 0).all(),
              "scalar diffusion must be finite and positive at every rule point");
        const Mat dx = Basis(rule, *cell, 0), dz = Basis(rule, *cell, 2);
        const Vec weights = rule.w.cwiseProduct(diffusion[e]);
        Mat local = dx.transpose() * weights.asDiagonal() * dx +
                    dz.transpose() * weights.asDiagonal() * dz;
        if (!reaction.empty())
        {
            // The net epsilon reaction may be negative; no positivity floor.
            Check(reaction[e].size() == rule.w.size() && reaction[e].allFinite(),
                  "invalid scalar reaction samples");
            const Mat b = Basis(rule, *cell, -1);
            const Vec rw = rule.w.cwiseProduct(reaction[e]);
            local.noalias() += b.transpose() * rw.asDiagonal() * b;
        }
        for (int i = 0; i < 16; ++i)
        {
            for (int j = 0; j < 16; ++j)
            {
                triplets.emplace_back(cell->ids[i], cell->ids[j], local(i, j));
            }
        }
    }
    Sparse result(Size(), Size());
    result.setFromTriplets(triplets.begin(), triplets.end());
    return result;
}

Vec Scalar::DiffusionAction(const Rules &rules,
                            const std::vector<Vec> &diffusion,
                            const Vec &state,
                            const std::vector<Vec> &reaction) const
{
    Check(state.size() == Size() && state.allFinite() &&
              diffusion.size() == rules.size() &&
              (reaction.empty() || reaction.size() == rules.size()),
          "invalid scalar diffusion action inventory");
    Vec result = Vec::Zero(Size());
    for (int e = 0; e < int(rules.size()); ++e)
    {
        Check(bool(rules[e]), "null scalar diffusion action rule");
        const auto &rule = *rules[e];
        const Cell *cell = Find(rule);
        if (!cell) { continue; }
        Check(diffusion[e].size() == rule.w.size() && diffusion[e].allFinite() &&
                  (diffusion[e].array() > 0).all(),
              "scalar diffusion action requires positive finite samples");
        const Vec local = Local(state, cell->ids);
        const Mat dx = Basis(rule, *cell, 0), dz = Basis(rule, *cell, 2);
        const Vec weights = rule.w.cwiseProduct(diffusion[e]);
        Vec action = dx.transpose() * weights.cwiseProduct(dx * local) +
                     dz.transpose() * weights.cwiseProduct(dz * local);
        if (!reaction.empty())
        {
            Check(reaction[e].size() == rule.w.size() && reaction[e].allFinite(),
                  "invalid scalar diffusion action reaction samples");
            const Mat b = Basis(rule, *cell, -1);
            action.noalias() += b.transpose() *
                rule.w.cwiseProduct(reaction[e]).cwiseProduct(b * local);
        }
        Add(result, action, *cell);
    }
    return result;
}

Vec Scalar::Source(const Rules &rules, const std::vector<Vec> &source) const
{
    Check(source.size() == rules.size(), "scalar source rule count differs");
    Vec result = Vec::Zero(Size());
    for (int e = 0; e < int(rules.size()); ++e)
    {
        Check(bool(rules[e]), "null scalar source rule");
        const auto &rule = *rules[e];
        const Cell *cell = Find(rule);
        if (!cell)
        {
            continue;
        }
        Check(source[e].size() == rule.w.size() && source[e].allFinite(),
              "invalid scalar source samples");
        const Mat b = Basis(rule, *cell, -1);
        Add(result, b.transpose() * rule.w.cwiseProduct(source[e]), *cell);
    }
    return result;
}

Vec Scalar::Advection(const Rules &rules, const Vec &velocity,
                      const Vec &scalar) const
{
    const int nu = m_flow.VelocitySize();
    Check((velocity.size() == 2 * nu ||
           velocity.size() == 2 * nu + m_flow.PressureSize()) &&
              velocity.allFinite() && scalar.size() == Size() && scalar.allFinite(),
          "native scalar advection state dimension/data differ");
    auto sampleVelocity = [&](const NativeRule &rule, int component) {
        Vec local(16);
        for (int j = 0; j < 16; ++j)
        {
            local[j] = velocity[component * nu + rule.u[j]];
        }
        return Vec(rule.B * local);
    };
    Vec result = Vec::Zero(Size());
    for (const auto &pointer : rules)
    {
        Check(bool(pointer), "null scalar advection rule");
        const auto &rule = *pointer;
        const Cell *cell = Find(rule);
        if (!cell)
        {
            continue;
        }
        const Mat b = Basis(rule, *cell, -1);
        const Mat dx = Basis(rule, *cell, 0), dz = Basis(rule, *cell, 2);
        const Vec weightedScalar = rule.w.cwiseProduct(b * Local(scalar, cell->ids));
        const Vec u = sampleVelocity(rule, 0), v = sampleVelocity(rule, 1);
        Add(result, -dx.transpose() * weightedScalar.cwiseProduct(u) -
                    dz.transpose() * weightedScalar.cwiseProduct(v), *cell);
    }
    // Interior CG face fluxes cancel in the assembled weak form. Only the
    // scalar domain's actual exterior remains, including its essential faces.
    for (const auto &face : m_faces)
    {
        const Cell &cell = m_cells[face.cell];
        const NativeRule &rule = face.rule;
        const Mat b = Basis(rule, cell, -1);
        const Vec normalVelocity = rule.NormalVelocity(sampleVelocity(rule,0),
                                                       sampleVelocity(rule,1));
        const Vec q = b * Local(scalar, cell.ids);
        Add(result, b.transpose() *
            rule.w.cwiseProduct(normalVelocity).cwiseProduct(q), cell);
    }
    return result;
}

Sparse Scalar::TransportMatrix(const Rules&rules,const Vec&velocity)const
{
    const int nu=m_flow.VelocitySize();
    Check((velocity.size()==2*nu||velocity.size()==2*nu+m_flow.PressureSize())&&velocity.allFinite(),
          "invalid frozen scalar transport velocity");
    auto sample=[&](const NativeRule&r,int component){Vec local(16);for(int j=0;j<16;++j)local[j]=velocity[component*nu+r.u[j]];return Vec(r.B*local);};
    std::vector<Eigen::Triplet<double>>triplets;
    auto append=[&](const Cell&cell,const Mat&matrix){for(int i=0;i<16;++i)for(int j=0;j<16;++j)triplets.emplace_back(cell.ids[i],cell.ids[j],matrix(i,j));};
    for(const auto&pointer:rules){Check(bool(pointer),"null scalar transport rule");const auto&r=*pointer;const Cell*cell=Find(r);if(!cell)continue;
        const Mat b=Basis(r,*cell,-1),dx=Basis(r,*cell,0),dz=Basis(r,*cell,2);
        append(*cell,-dx.transpose()*r.w.cwiseProduct(sample(r,0)).asDiagonal()*b-
                     dz.transpose()*r.w.cwiseProduct(sample(r,1)).asDiagonal()*b);
    }
    for(const auto&face:m_faces){const Cell&cell=m_cells[face.cell];const auto&r=face.rule;const Mat b=Basis(r,cell,-1);
        const Vec normal=r.NormalVelocity(sample(r,0),sample(r,1));
        append(cell,b.transpose()*r.w.cwiseProduct(normal).asDiagonal()*b);
    }
    Sparse result(Size(),Size());result.setFromTriplets(triplets.begin(),triplets.end());return result;
}

Vec Scalar::NaturalFlux(const ScalarBoundaryFunction &outwardFlux,
                        const RuleFactory &factory,
                        const std::function<double(const Element&,double)>&cutoff) const
{
    Check(bool(outwardFlux), "missing scalar natural flux callback");
    Vec result = Vec::Zero(Size());
    for (const auto &face : m_faces)
    {
        const Cell &cell = m_cells[face.cell];
        const NativeRule rule = factory.Face(*cell.flow, face.side, 16,cutoff);
        const Mat b = Basis(rule, cell, -1);
        Vec flux(rule.w.size());
        for (int q = 0; q < flux.size(); ++q)
        {
            flux[q] = outwardFlux(rule.x[q], rule.z[q], rule.Outward(q));
            Check(std::isfinite(flux[q]), "nonfinite scalar natural flux");
        }
        Add(result, -b.transpose() * rule.w.cwiseProduct(flux), cell);
    }
    return result;
}

double Scalar::FreeDual(const Vec &residual) const
{
    Check(residual.size() == FreeSize() && residual.allFinite(),
          "invalid scalar free residual");
    Vec represented = m_riesz.solve(residual);
    Check(m_riesz.info() == Eigen::Success && represented.allFinite(),
          "native scalar residual Riesz solve failed");
    const double square = residual.dot(represented);
    Check(std::isfinite(square) && square >= 0,
          "invalid native scalar residual dual norm");
    return std::sqrt(square / m_volume);
}

double Scalar::Dual(const Vec &fullResidual) const
{
    return FreeDual(Restrict(fullResidual));
}

ScalarSolveInfo Scalar::Residual(const Sparse &current, const Vec &rhs,
                                 const Vec &free) const
{
    CheckSparse(current, FreeSize());
    Check(rhs.size() == FreeSize() && free.size() == FreeSize() &&
              rhs.allFinite() && free.allFinite(),
          "native scalar residual state dimension/data differ");
    ScalarSolveInfo info;
    const double began = Clock();
    Vec raw = current * free;
    raw -= rhs;
    info.residual = FreeDual(raw);
    info.residualSeconds = Clock() - began;
    return info;
}

double Scalar::Prepare(const Sparse &current)
{
    const double began = Clock();
    m_factor = std::make_unique<Eigen::SparseLU<Sparse>>();
    m_factor->compute(current);
    RequireLinear(m_factor->info() == Eigen::Success,
          "native scalar current-operator preconditioner factor failed");
    m_initialIterations = m_lastIterations = 0;
    return Clock() - began;
}

Vec Scalar::Precondition(const Vec &rhs) const
{
    Check(bool(m_factor), "missing native scalar preconditioner");
    Vec result = m_factor->solve(rhs);
    RequireLinear(m_factor->info() == Eigen::Success && result.allFinite(),
          "native scalar preconditioner application failed");
    return result;
}

ScalarSolveInfo Scalar::Solve(const Sparse &current, const Vec &rhs, Vec &free,
                             double coefficientChange, double tolerance)
{
    CheckSparse(current, FreeSize());
    Check(rhs.size() == FreeSize() && free.size() == FreeSize() &&
              rhs.allFinite() && free.allFinite() &&
              std::isfinite(coefficientChange) && coefficientChange >= 0 &&
              std::isfinite(tolerance) && tolerance > 0,
          "invalid native scalar solve arguments");
    ScalarSolveInfo info;
    const double began = Clock();
    // Only the preconditioner is retained. Every Krylov product and final
    // residual below uses the supplied CURRENT alpha*M + K operator.
    if (!m_factor || coefficientChange > .25 || m_lastIterations > 30 ||
        (m_initialIterations > 0 && m_lastIterations > 2 * m_initialIterations))
    {
        info.factorSeconds += Prepare(current);
        ++info.rebuilds;
    }
    constexpr int maximumIterations = 400, restart = 30;
    int used = 0;
    double trueResidual = Residual(current, rhs, free).residual;
    while (trueResidual > tolerance && used < maximumIterations)
    {
        // Form the SAME explicit current-matrix residual used by the final
        // acceptance. Fused rhs-A*x and A*x-rhs expressions can round
        // differently at a tight physical tolerance.
        Vec raw = current * free;
        raw -= rhs;
        const Vec residual = -raw.cwiseProduct(m_scale);
        const double beta = residual.norm();
        RequireLinear(std::isfinite(beta) && beta > 0, "invalid scalar Krylov residual");
        const int dimension = std::min(restart, maximumIterations - used);
        Mat basis = Mat::Zero(FreeSize(), dimension + 1);
        Mat preconditioned(FreeSize(), dimension);
        Mat hessenberg = Mat::Zero(dimension + 1, dimension);
        basis.col(0) = residual / beta;
        Vec leastSquaresRhs = Vec::Zero(dimension + 1);
        leastSquaresRhs[0] = beta;
        Vec best;
        int done = 0;
        for (int j = 0; j < dimension; ++j)
        {
            preconditioned.col(j) =
                Precondition(basis.col(j).cwiseQuotient(m_scale));
            Vec action = (current * preconditioned.col(j)).cwiseProduct(m_scale);
            const double actionNorm = action.norm();
            // Reorthogonalised Arnoldi avoids silently losing independence
            // for nearly exact retained LU preconditioners.
            for (int pass = 0; pass < 2; ++pass)
            {
                for (int i = 0; i <= j; ++i)
                {
                    const double value = basis.col(i).dot(action);
                    hessenberg(i, j) += value;
                    action.noalias() -= value * basis.col(i);
                }
            }
            hessenberg(j + 1, j) = action.norm();
            const bool breakdown = hessenberg(j + 1, j) <=
                32 * std::numeric_limits<double>::epsilon() *
                    std::max(actionNorm, std::numeric_limits<double>::min());
            if (!breakdown)
            {
                basis.col(j + 1) = action / hessenberg(j + 1, j);
            }
            best = hessenberg.block(0, 0, j + 2, j + 1)
                       .colPivHouseholderQr().solve(leastSquaresRhs.head(j + 2));
            RequireLinear(best.allFinite(), "nonfinite scalar GMRES least-squares solve");
            ++used;
            done = j + 1;
            const double estimate =
                (leastSquaresRhs.head(j + 2) -
                 hessenberg.block(0, 0, j + 2, j + 1) * best).norm();
            if (estimate <= tolerance * std::sqrt(m_volume) || breakdown)
            {
                break;
            }
        }
        free.noalias() += preconditioned.leftCols(done) * best;
        trueResidual = Residual(current, rhs, free).residual;
    }
    const auto final = Residual(current, rhs, free);
    info.residual = final.residual;
    info.residualSeconds = final.residualSeconds;
    info.iterations = used;
    info.solveSeconds = Clock() - began - info.factorSeconds - info.residualSeconds;
    if(info.residual>tolerance){std::ostringstream error;error<<std::setprecision(17)
        <<"native scalar CURRENT-operator reference-scaled residual failed: residual="<<info.residual
        <<", tolerance="<<tolerance<<", iterations="<<used<<", rebuilds="<<info.rebuilds;
        throw NativeLinearFailure(error.str());}
    if (m_initialIterations == 0 && used > 0)
    {
        m_initialIterations = used;
    }
    m_lastIterations = used;
    return info;
}
}
