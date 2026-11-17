// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/WJNativeFlowTangent.h>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>

namespace Nektar::WJNative
{
namespace
{
void Require(bool condition,const char*message)
{if(!condition)throw std::invalid_argument(message);}
constexpr std::array<int,2> Axis{{0,2}};
void Finite(const FrozenPlanarVelocityTangent&t)
{for(const auto&column:t.dVdG)for(const auto&row:column)for(double value:row)
    Require(std::isfinite(value),"nonfinite planar velocity tangent");}
}

FrozenPlanarVelocityTangent EvaluatePlanarVelocityTangent(double k,double epsilon,
    const Tensor&gradient,double nu,double auxiliary,double r,double H,
    ConstitutiveBackend backend,WallEvaluationStats*statistics)
{
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)
    {
        Require(std::isfinite(gradient[i][j]),"nonfinite planar tangent gradient");
        if(i==1||j==1)Require(gradient[i][j]==0.,"planar velocity tangent requires zero spanwise gradient entries");
    }
    FrozenPlanarVelocityTangent out;
    for(int component=0;component<2;++component)for(int derivative=0;derivative<2;++derivative)
    {
        Tensor direction{};direction[Axis[component]][Axis[derivative]]=1.;
        out.dVdG[2*component+derivative]=ApplyAmendedDirection(k,epsilon,gradient,{0,0,1},
            nu,auxiliary,r,H,0.,0.,direction,0.,0.,backend,statistics).V;
    }
    Finite(out);return out;
}

Sparse AssemblePlanarVelocityTangent(const Rules&rules,int velocitySize,
    const std::vector<OwnedPlanarVelocityTangents>&samples)
{
    Require(velocitySize>0&&velocitySize<=std::numeric_limits<int>::max()/2,
            "invalid planar tangent velocity inventory");
    Require(!rules.empty()&&samples.size()==rules.size(),"planar tangent rule/sample count differs");
    std::vector<Eigen::Triplet<double>>triplets;triplets.reserve(rules.size()*1024);
    for(std::size_t element=0;element<rules.size();++element)
    {
        Require(bool(rules[element]),"null planar tangent rule");
        const auto&r=*rules[element];const auto&t=samples[element];const auto nq=r.w.size();
        Require(r.geometry>=0&&t.geometry==r.geometry&&t.points.size()==std::size_t(nq),
                "planar tangent owner/point inventory differs");
        Require(nq>0&&r.w.allFinite()&&(r.w.array()>0).all()&&r.x.size()==nq&&r.z.size()==nq&&
                r.x.allFinite()&&r.z.allFinite(),"invalid planar tangent integration rule");
        Require(r.u.size()==16&&r.us.size()==16&&std::set<int>(r.u.begin(),r.u.end()).size()==16,
                "planar tangent requires a native Q3 signed velocity map");
        for(int j=0;j<16;++j)Require(r.u[j]>=0&&r.u[j]<velocitySize&&(r.us[j]==1.||r.us[j]==-1.),
                "invalid planar tangent signed global index");
        for(const Mat*matrix:{&r.B,&r.Dx,&r.Dz})Require(matrix->rows()==nq&&matrix->cols()==16&&matrix->allFinite(),
                "planar tangent native basis dimensions/data differ");
        for(const auto&point:t.points)Finite(point);
        const std::array<const Mat*,2>gradient{{&r.Dx,&r.Dz}};
        for(int row=0;row<2;++row)for(int column=0;column<2;++column)
        {
            Mat block=Mat::Zero(16,16);
            for(int testDerivative=0;testDerivative<2;++testDerivative)
                for(int trialDerivative=0;trialDerivative<2;++trialDerivative)
                {
                    Vec weight(nq);
                    for(Eigen::Index q=0;q<nq;++q)weight[q]=r.w[q]*
                        t.points[q].dVdG[2*column+trialDerivative][Axis[row]][Axis[testDerivative]];
                    block.noalias()+=gradient[testDerivative]->transpose()*weight.asDiagonal()*(*gradient[trialDerivative]);
                }
            Require(block.allFinite(),"nonfinite assembled planar velocity tangent");
            for(int i=0;i<16;++i)for(int j=0;j<16;++j)
                triplets.emplace_back(row*velocitySize+r.u[i],column*velocitySize+r.u[j],block(i,j));
        }
    }
    Sparse matrix(2*velocitySize,2*velocitySize);matrix.setFromTriplets(triplets.begin(),triplets.end());
    return matrix;
}
}
