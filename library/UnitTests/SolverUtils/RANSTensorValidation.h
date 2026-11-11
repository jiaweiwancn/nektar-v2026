// SPDX-License-Identifier: MIT
#ifndef NEKTAR_RANS_TENSOR_VALIDATION_H
#define NEKTAR_RANS_TENSOR_VALIDATION_H
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <numeric>
#include <set>

namespace Nektar
{
// Test-only independent native transforms. This never contributes a governing
// load, changes a state/history, or substitutes a test basis after assembly.
inline std::string RANSTensorValidation(RANSWallEvolutionPatch &patch)
{
    namespace J=RANSWallSolveJson;
    using Vector=RANSWallAlgebra::Vector;
    using Matrix=RANSWallAlgebra::DenseMatrix;
    RANSWallAlgebra::Require(patch.Dimension()==2 && !patch.Mapped() && patch.Size()<500,
        "tensor validation requires a small serial affine 2D fixture");
    std::vector<std::string> fields;
    for(int b=0;b<patch.Gauge();++b)
    {
        const auto &field=*patch.fields[b];const auto native=field.Field();
        const auto cg=native->GetLocalToGlobalMap();
        Vector defects(8);int modeChecks=0;Matrix independentMass(field.Size(),Vector(field.Size()));
        for(int e=0;e<native->GetExpSize();++e)
        {
            const auto exp=native->GetExp(e);const auto &el=field.Elements()[e];
            const int nc=exp->GetNcoeffs(),nq=exp->GetTotPoints(),start=native->GetCoeff_Offset(e);
            const double lx=el.upper[0]-el.lower[0],lz=el.upper[2]-el.lower[2];
            std::vector<Array<OneD,NekDouble>> modes;
            for(int j=0;j<nc;++j)
            {
                Array<OneD,NekDouble> c(nc,0.),p(nq),dx(nq),dz(nq),dxz(nq),dzz(nq);c[j]=1;
                exp->BwdTrans(c,p);modes.push_back(p);exp->PhysDeriv(0,p,dx);exp->PhysDeriv(1,p,dz);
                exp->PhysDeriv(1,dx,dxz);exp->PhysDeriv(1,dz,dzz);
                Vector state(field.Size());state[cg->GetLocalToGlobalMap(start+j)]=cg->GetLocalToGlobalSign(start+j);
                for(double x:{-1.,-.73,.11,.66,1.})for(double z:{-1.,-.82,.27,.74,1.})
                {
                    const RANSWallCoordinate physical{el.lower[0]+(x+1)*lx/2,0,el.lower[2]+(z+1)*lz/2};
                    const auto row=field.ElementPoint(e,physical);Array<OneD,NekDouble> local(2);
                    local[0]=row.rawLocalCoordinate[0];local[1]=row.rawLocalCoordinate[2];
                    for(int d=0;d<2;++d)if(std::abs(local[d])>1)local[d]=local[d]<0?-1:1;
                    const auto g=row.Gradient(state);const auto h=field.RadialGradientDerivative(row);
                    auto contract=[&](const Vector &r){double y=0;for(int i=0;i<row.ids.size();++i)y+=r[i]*state[row.ids[i]];return y;};
                    defects[0]=std::max(defects[0],std::abs(row.Value(state)-exp->StdPhysEvaluate(local,p)));
                    defects[1]=std::max(defects[1],lx*std::abs(g[0]-exp->StdPhysEvaluate(local,dx)));
                    defects[2]=std::max(defects[2],lz*std::abs(g[2]-exp->StdPhysEvaluate(local,dz)));
                    defects[3]=std::max(defects[3],lx*lz*std::abs(contract(h[0])-exp->StdPhysEvaluate(local,dxz)));
                    defects[4]=std::max(defects[4],lz*lz*std::abs(contract(h[2])-exp->StdPhysEvaluate(local,dzz)));
                    ++modeChecks;
                }
            }
            Array<OneD,NekDouble> product(nq);
            for(int i=0;i<nc;++i)for(int j=0;j<nc;++j)
            {
                for(int q=0;q<nq;++q)product[q]=modes[i][q]*modes[j][q];
                independentMass[cg->GetLocalToGlobalMap(start+i)][cg->GetLocalToGlobalMap(start+j)]+=
                    cg->GetLocalToGlobalSign(start+i)*cg->GetLocalToGlobalSign(start+j)*exp->Integral(product);
            }
        }
        const auto mass=field.Mass();
        for(int i=0;i<field.Size();++i)for(int j=0;j<field.Size();++j)
            defects[5]=std::max(defects[5],std::abs(mass[i][j]-independentMass[i][j]));
        for(int kind=0;kind<3;++kind)
        {
            const RANSWallFunction f=[kind](const RANSWallCoordinate &p)
            {return kind==0?1.:kind==1?1.+.25*p[0]-.3*p[2]:1.+p[0]*p[2]+.125*p[0]*p[0]*p[2]*p[2];};
            const auto c=field.Project(f);const auto copy=field.RestrictCoefficients(field,c);
            RANSWallAlgebra::Require(c==copy,"signed native self restriction differs");
            for(const auto &q:field.Quadrature(7))
            {
                const auto n=field.NativeEvaluate(c,q.coordinate);const auto g=q.Gradient(c);
                defects[6]=std::max({defects[6],std::abs(q.Value(c)-f(q.coordinate)),
                    std::abs(n[0]-q.Value(c)),std::abs(n[1]-g[0]),std::abs(n[3]-g[2])});
            }
            for(int axis:{0,2})for(int side:{0,1})
            {
                const auto boundary=field.Boundary(f,[&](int a,int s,const RANSWallGlobalField::Element&){return a==axis && s==side;});
                Vector trace(field.Size());for(auto [id,value]:boundary)trace[id]=value;
                for(int e=0;e<field.Elements().size();++e)
                {
                    const auto &el=field.Elements()[e];const int free=axis==0?2:0;
                    for(int j=0;j<=20;++j)
                    {
                        RANSWallCoordinate p{};p[axis]=side?el.upper[axis]:el.lower[axis];
                        p[free]=el.lower[free]+(el.upper[free]-el.lower[free])*j/20.;
                        defects[7]=std::max(defects[7],std::abs(field.ElementPoint(e,p).Value(trace)-f(p)));
                    }
                }
            }
        }
        const auto radial=field.RadialProject([](double z){return std::exp(-.2*z);},{});
        std::vector<std::string> samples;
        for(const auto &el:field.Elements())for(double t:{.13,.41,.83})
        {
            const RANSWallCoordinate p{.23*el.lower[0]+.77*el.upper[0],0,(1-t)*el.lower[2]+t*el.upper[2]};
            samples.push_back(J::Numbers(Vector{p[2],field.Point(p).Value(radial)}));
        }
        const auto constants=field.Project([](const RANSWallCoordinate&){return 1.;});
        const auto mc=RANSWallAlgebra::Multiply(mass,constants);
        const double volume=std::inner_product(constants.begin(),constants.end(),mc.begin(),0.);
        fields.push_back(J::Object({{"block",J::Number(b)},{"size",J::Number(field.Size())},
            {"axis_degrees",J::Numbers(field.Elements().front().degrees)},
            {"mode_checks",J::Number(modeChecks)},{"defects_value_dx_dz_dxz_dzz_mass_polynomial_faces",J::Numbers(defects)},
            {"mass",J::Rows(mass)},{"radial_projection_samples",J::List(samples)},{"integrated_volume",J::Number(volume)}}));
    }
    return J::Object({{"schema",J::Quote("NativeDirectionalTensorChecksV1")},{"fields",J::List(fields)},
        {"scope",J::Quote("independent Nektar BwdTrans/PhysDeriv/Integral, signed CG packing, native faces and restriction; all individual modes; constants, affine and mixed polynomials")}});
}
}
#endif
