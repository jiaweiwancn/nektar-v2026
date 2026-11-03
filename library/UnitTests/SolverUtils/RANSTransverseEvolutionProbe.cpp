///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallTimeIntegration.h>
#include <SolverUtils/Turbulence/EquilibriumAWFProfiles.h>
#include <LibUtilities/LinearAlgebra/Lapack.hpp>
#include <StdRegions/StdSegExp.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Nektar;
using namespace Nektar::SolverUtils;
namespace
{
using Vector=RANSWallAlgebra::Vector;
using Matrix=RANSWallAlgebra::DenseMatrix;
using Engine=RANSWallTimeIntegration;
constexpr double height=2,lower=4,upper=12,ustar=.4,finalTime=.05,spatialDt=.00015625;
int checks=0;
void Check(bool good,const std::string &message)
{++checks;if(!good)throw std::runtime_error("Transverse probe: "+message);}
void Write(std::ostream &out,const Vector &v)
{
    out<<'[';
    for(std::size_t j=0;j<v.size();++j){if(j)out<<',';Check(std::isfinite(v[j]),"nonfinite report value");out<<v[j];}
    out<<']';
}
void Write(std::ostream &out,const Matrix &v)
{
    out<<'[';for(std::size_t j=0;j<v.size();++j){if(j)out<<',';Write(out,v[j]);}out<<']';
}
double Dot(const Vector &a,const Vector &b)
{long double sum=0;for(std::size_t j=0;j<a.size();++j)sum+=static_cast<long double>(a[j])*b[j];return sum;}
Vector Multiply(const Matrix &a,const Vector &b){return RANSWallAlgebra::Multiply(a,b);}
struct Factor
{
    int n;Vector entries;std::vector<int> pivots;
    explicit Factor(const Matrix &matrix):n(matrix.size()),entries(n*n),pivots(n)
    {
        for(int i=0;i<n;++i)for(int j=0;j<n;++j)entries[i+n*j]=matrix[i][j];
        int info=0;Lapack::Dgetrf(n,n,entries.data(),n,pivots.data(),info);Check(info==0,"singular Galerkin system");
    }
    Vector Apply(Vector rhs)const
    {
        int info=0;Lapack::Dgetrs('N',n,1,entries.data(),n,pivots.data(),rhs.data(),n,info);
        Check(info==0,"failed Galerkin solve");return rhs;
    }
};
double MinimumEigenvalue(RANSTensor a)
{
    for(int iteration=0;iteration<40;++iteration)
    {
        int p=0,q=1;
        for(int i=0;i<3;++i)for(int j=i+1;j<3;++j)if(std::abs(a[i][j])>std::abs(a[p][q])){p=i;q=j;}
        if(std::abs(a[p][q])<1e-16)break;
        const double angle=.5*std::atan2(2*a[p][q],a[q][q]-a[p][p]),c=std::cos(angle),s=std::sin(angle);
        const double app=c*c*a[p][p]-2*s*c*a[p][q]+s*s*a[q][q],
            aqq=s*s*a[p][p]+2*s*c*a[p][q]+c*c*a[q][q];
        for(int k=0;k<3;++k)if(k!=p && k!=q)
        {
            const double ap=c*a[k][p]-s*a[k][q],aq=s*a[k][p]+c*a[k][q];
            a[k][p]=a[p][k]=ap;a[k][q]=a[q][k]=aq;
        }
        a[p][p]=app;a[q][q]=aqq;a[p][q]=a[q][p]=0;
    }
    return std::min({a[0][0],a[1][1],a[2][2]});
}
struct Tangent
{
    std::array<std::array<double,2>,2> block{};
    double nut=0,production=0,nuDerivative=0,productionDerivative=0,kDerivative=0,epsilonDerivative=0;
    double covarianceMinimum=0;
};
Tangent At(const RANSNeutralReference &model,double z,double matching,double relativeStep=1e-6,double eta=0)
{
    const auto r=model.AtHeight(z,ustar,height,matching);RANSTensor G{};G[0][2]=r.velocityGradient;
    const auto base=EvaluateRANSClosure(RANSModel::CLSHighRe,r.k,r.epsilon,G,z);
    Tangent result;result.nut=base.eddyViscosity;result.production=base.production;
    result.covarianceMinimum=MinimumEigenvalue(base.reynoldsStress);
    const double h=relativeStep*r.velocityGradient;
    for(int column=0;column<2;++column)
    {
        const int axis=column?2:0;auto plus=G,minus=G;plus[1][axis]+=h;minus[1][axis]-=h;
        const auto a=EvaluateRANSClosure(RANSModel::CLSHighRe,r.k,r.epsilon,plus,z),
            b=EvaluateRANSClosure(RANSModel::CLSHighRe,r.k,r.epsilon,minus,z);
        for(int row=0;row<2;++row)
        {
            const int direction=row?2:0;
            const double anchor=-eta*base.reynoldsStress[1][direction];
            result.block[row][column]=((-a.reynoldsStress[1][direction]-anchor)-
                (-b.reynoldsStress[1][direction]-anchor))/(2*h);
        }
        result.nuDerivative=std::max(result.nuDerivative,std::abs((a.eddyViscosity-b.eddyViscosity)/(2*h)));
        result.productionDerivative=std::max(result.productionDerivative,std::abs((a.production-b.production)/(2*h)));
        result.kDerivative=std::max(result.kDerivative,std::abs((a.kSource-b.kSource)/(2*h)));
        result.epsilonDerivative=std::max(result.epsilonDerivative,std::abs((a.epsilonSource-b.epsilonSource)/(2*h)));
    }
    return result;
}
void Write(std::ostream &out,const Tangent &a)
{
    out<<"{\"block\":";Write(out,Matrix{{a.block[0][0],a.block[0][1]},{a.block[1][0],a.block[1][1]}});
    out<<",\"nut\":"<<a.nut<<",\"production\":"<<a.production<<",\"covariance_minimum\":"<<a.covarianceMinimum
        <<",\"linear_turbulence_forcing\":";Write(out,Vector{a.nuDerivative,a.productionDerivative,a.kDerivative,a.epsilonDerivative});out<<'}';
}
double Envelope(double z)
{
    const double s=(z-lower)/(upper-lower);
    return s<=0 || s>=1?0:std::exp(4-1/(s*(1-s)));
}
struct Point
{
    double z,weight;std::vector<int> ids;Vector phi,derivative;
};
struct Space
{
    int degree,cells,size;std::vector<int> free,index;std::vector<Point> points;
    std::vector<std::vector<int>> maps;
    Space(int p,int elements):degree(p),cells(elements),size(elements*p+1),index(size,-1)
    {
        for(int i=0;i<size;++i)if(i!=0 && i!=cells){index[i]=free.size();free.push_back(i);}
        const int order=std::max(24,2*p+6);
        StdRegions::StdSegExp segment({LibUtilities::eModified_A,static_cast<std::size_t>(p+1),
            {static_cast<std::size_t>(order),LibUtilities::eGaussGaussLegendre}});
        const auto nodes=segment.GetBasis(0)->GetZ(),weights=segment.GetBasis(0)->GetW();
        std::vector<Vector> phi(order,Vector(p+1)),derivative(order,Vector(p+1));
        for(int mode=0;mode<=p;++mode)
        {
            Array<OneD,NekDouble> coefficients(p+1,0.),value(order),gradient(order);coefficients[mode]=1;
            segment.BwdTrans(coefficients,value);segment.PhysDeriv(0,value,gradient);
            for(int q=0;q<order;++q){phi[q][mode]=value[q];derivative[q][mode]=gradient[q];}
        }
        const double width=(upper-lower)/cells;
        for(int e=0;e<cells;++e)
        {
            std::vector<int> ids{e,e+1};for(int j=2;j<=p;++j)ids.push_back(cells+1+e*(p-1)+j-2);maps.push_back(ids);
            for(int q=0;q<order;++q)
            {
                auto dz=derivative[q];for(double &v:dz)v*=2/width;
                points.push_back({lower+width*(e+(nodes[q]+1)/2),weights[q]*width/2,ids,phi[q],std::move(dz)});
            }
        }
    }
    double Value(const Vector &coefficients,const Point &q,bool gradient=false,int part=0)const
    {
        long double value=0;const auto &basis=gradient?q.derivative:q.phi;
        for(std::size_t i=0;i<q.ids.size();++i)if(index[q.ids[i]]>=0)
            value+=static_cast<long double>(basis[i])*coefficients[part*free.size()+index[q.ids[i]]];
        return value;
    }
};
struct Operator
{
    std::string name;Space space;double wave,matching;
    Matrix mass,real,imaginary,pairedMass,load;Vector initial;
    double productionDerivative=0,nuDerivative=0,kDerivative=0,epsilonDerivative=0,covarianceMinimum=1e100;
    Operator(const RANSNeutralReference &model,int p,int cells,double xi,double match,double eta=0)
        :name("p"+std::to_string(p)+"-n"+std::to_string(cells)+"-xi"+std::to_string(static_cast<int>(xi))),
         space(p,cells),wave(xi),matching(match)
    {
        const int n=space.free.size();mass=real=imaginary=Matrix(n,Vector(n));Vector rhs(n);
        for(const auto &q:space.points)
        {
            const auto tangent=At(model,q.z,matching,1e-6,eta),base=model.AtHeight(q.z,ustar,height,matching);
            covarianceMinimum=std::min(covarianceMinimum,tangent.covarianceMinimum);
            productionDerivative=std::max(productionDerivative,tangent.productionDerivative);
            nuDerivative=std::max(nuDerivative,tangent.nuDerivative);kDerivative=std::max(kDerivative,tangent.kDerivative);
            epsilonDerivative=std::max(epsilonDerivative,tangent.epsilonDerivative);
            for(std::size_t i=0;i<q.ids.size();++i)if(space.index[q.ids[i]]>=0)
            {
                const int row=space.index[q.ids[i]];rhs[row]+=q.weight*q.phi[i]*Envelope(q.z);
                for(std::size_t j=0;j<q.ids.size();++j)if(space.index[q.ids[j]]>=0)
                {
                    const int column=space.index[q.ids[j]];const double m=q.weight*q.phi[i]*q.phi[j];
                    mass[row][column]+=m;
                    real[row][column]-=wave*wave*tangent.block[0][0]*m+q.weight*q.derivative[i]*tangent.block[1][1]*q.derivative[j];
                    imaginary[row][column]+=wave*q.weight*(q.phi[i]*tangent.block[0][1]*q.derivative[j]-
                        q.derivative[i]*tangent.block[1][0]*q.phi[j]-q.phi[i]*base.velocity*q.phi[j]);
                }
            }
        }
        initial=Factor(mass).Apply(rhs);initial.resize(2*n,0.);
        pairedMass=load=Matrix(2*n,Vector(2*n));
        for(int i=0;i<n;++i)for(int j=0;j<n;++j)
        {
            pairedMass[i][j]=pairedMass[i+n][j+n]=mass[i][j];
            load[i][j]=load[i+n][j+n]=real[i][j];load[i][j+n]=-imaginary[i][j];load[i+n][j]=imaginary[i][j];
        }
    }
    double Energy(const Vector &state)const{return .5*Dot(state,Multiply(pairedMass,state));}
    double EnergyRate(const Vector &state)const{return Dot(state,Multiply(load,state));}
    std::array<double,2> Covariance(const RANSNeutralReference &model,const Vector &state)const
    {
        double minimum=1e100,velocity=0;
        for(const auto &q:space.points)
        {
            const auto r=model.AtHeight(q.z,ustar,height,matching);
            const double c=space.Value(state,q),s=space.Value(state,q,false,1),
                dc=space.Value(state,q,true),ds=space.Value(state,q,true,1);
            for(double phase:{0.,.25*std::acos(-1.),.5*std::acos(-1.),.75*std::acos(-1.)})
            {
                RANSTensor g{};g[0][2]=r.velocityGradient;
                g[1][0]=-wave*(c*std::sin(phase)+s*std::cos(phase));g[1][2]=dc*std::cos(phase)-ds*std::sin(phase);
                const auto closure=EvaluateRANSClosure(RANSModel::CLSHighRe,r.k,r.epsilon,g,q.z);
                minimum=std::min(minimum,MinimumEigenvalue(closure.reynoldsStress));
                velocity=std::max(velocity,std::abs(c*std::cos(phase)-s*std::sin(phase)));
            }
        }
        return {minimum,velocity};
    }
};
void Describe(std::ostream &out,const Operator &op)
{
    out<<"{\"id\":\""<<op.name<<"\",\"degree\":"<<op.space.degree<<",\"cells\":"<<op.space.cells
        <<",\"xi\":"<<op.wave<<",\"size\":"<<op.space.size<<",\"free_ids\":";
    Write(out,Vector(op.space.free.begin(),op.space.free.end()));out<<",\"maps\":[";
    for(std::size_t e=0;e<op.space.maps.size();++e){if(e)out<<',';Write(out,Vector(op.space.maps[e].begin(),op.space.maps[e].end()));}
    out<<"],\"mass\":";Write(out,op.mass);out<<",\"real_load\":";Write(out,op.real);
    out<<",\"imaginary_load\":";Write(out,op.imaginary);out<<",\"initial_unit_coefficients\":";Write(out,op.initial);
    out<<",\"initial_unit_energy\":"<<op.Energy(op.initial)<<",\"initial_unit_energy_rate\":"<<op.EnergyRate(op.initial)
        <<",\"base_covariance_minimum\":"<<op.covarianceMinimum<<",\"linear_turbulence_forcing\":";
    Write(out,Vector{op.nuDerivative,op.productionDerivative,op.kDerivative,op.epsilonDerivative});
    out<<",\"native_rows\":[";
    int count=0;
    for(std::size_t j=0;j<op.space.points.size();++j)if(j%24==0 || j%24==12)
    {
        const auto &q=op.space.points[j];if(count++)out<<',';
        out<<"{\"z\":"<<q.z<<",\"weight\":"<<q.weight<<",\"ids\":";Write(out,Vector(q.ids.begin(),q.ids.end()));
        out<<",\"phi\":";Write(out,q.phi);out<<",\"derivative\":";Write(out,q.derivative);out<<'}';
    }
    out<<"]}";
}
void Trial(std::ostream &out,const RANSNeutralReference &model,const Operator &op,int order,double dt,double amplitude)
{
    Vector initial=op.initial;for(double &v:initial)v*=amplitude;
    std::map<double,std::shared_ptr<Factor>> factors;Vector stageHistory;double stageTime=0,stageLambda=0,stageError=0;
    Engine engine(order,dt,[&](const Vector &history,double time,double lambda)
    {
        if(!factors.count(lambda))
        {
            auto matrix=op.pairedMass;for(std::size_t i=0;i<matrix.size();++i)for(std::size_t j=0;j<matrix.size();++j)matrix[i][j]-=lambda*op.load[i][j];
            factors[lambda]=std::make_shared<Factor>(matrix);
        }
        const auto value=factors.at(lambda)->Apply(Multiply(op.pairedMass,history));
        auto difference=value;for(std::size_t i=0;i<value.size();++i)difference[i]-=history[i];
        auto residual=Multiply(op.pairedMass,difference);const auto load=Multiply(op.load,value);stageError=0;
        for(std::size_t i=0;i<residual.size();++i)stageError=std::max(stageError,std::abs(residual[i]-lambda*load[i]));
        stageHistory=history;stageTime=time;stageLambda=lambda;return value;
    },[](const Vector &state,double){return Vector(state.size(),0.);});
    engine.Initialise(initial,0);const auto before=op.Covariance(model,initial);
    out<<"{\"operator\":\""<<op.name<<"\",\"order\":"<<order<<",\"dt\":"<<dt<<",\"amplitude\":"<<amplitude
        <<",\"initial_state\":";Write(out,initial);out<<",\"initial_covariance\":";Write(out,Vector{before[0],before[1]});out<<",\"trajectory\":[";
    const int steps=std::lround(finalTime/dt);Check(std::abs(steps*dt-finalTime)<1e-13,"nonintegral final time");
    Vector final=initial;
    for(int step=0;step<steps;++step)
    {
        final=engine.Advance();if(step)out<<',';
        out<<"{\"time\":"<<engine.Time()<<",\"energy\":"<<op.Energy(final)<<",\"energy_rate\":"<<op.EnergyRate(final)
            <<",\"stage_time\":"<<stageTime<<",\"lambda\":"<<stageLambda<<",\"stage_residual\":"<<stageError<<",\"stage_history\":";
        Write(out,stageHistory);out<<",\"state\":";Write(out,final);out<<'}';
    }
    const auto after=op.Covariance(model,final),snapshot=engine.Capture();
    out<<"],\"final_covariance\":";Write(out,Vector{after[0],after[1]});
    out<<",\"native_history\":";Write(out,snapshot.history);out<<",\"native_history_times\":";Write(out,snapshot.times);
    out<<",\"gain\":"<<op.Energy(final)/op.Energy(initial)<<'}';
}
} // namespace

int main(int argc,char **argv)
{
    try
    {
        Check(argc==3 && std::string(argv[1])=="--output","use --output REPORT.json");
        const auto model=RANSNeutralReference::CLS();
        AWFProfileOptions options;options.quadratureOrder=32;
        const AWFProfileGeometry geometry{1.5e-5,height,.1,1.,1.};
        const double K=model.Q()*ustar*ustar;
        const EquilibriumAWFProfiles profile(geometry,K,model,{},options);
        const double matching=profile.Reference().matchingVelocity;
        std::ofstream out(argv[2]);Check(static_cast<bool>(out),"cannot open report");out<<std::setprecision(17);
        out<<"{\"schema\":1,\"experiment\":\"reduced linearized Fourier-x/native-CG-z transverse original CLS\",\"parameters\":";
        Write(out,Vector{height,lower,upper,ustar,K,model.Q(),model.Cmu(),model.Kappa(),matching,finalTime,spatialDt});
        out<<",\"support\":\"v=0 at band endpoints and identically zero outside [4,12], strictly above H=2; no wall/interface evolution\""
            <<",\"phase_convention\":\"real physical perturbation is Re((c+i*s)*exp(i*xi*x))*e_y\""
            <<",\"energy_convention\":\"one half integral of squared complex Fourier amplitude over z\""
            <<",\"time_partition\":\"all advection and original tensor diffusion implicit; explicit callback zero\""
            <<",\"constitutive_samples\":[";
        int count=0;
        for(double z:{4.,8.,12.})for(double step:{1e-3,1e-4,1e-5,1e-6})
        {
            if(count++)out<<',';
            out<<"{\"z\":"<<z<<",\"relative_step\":"<<step<<",\"value\":";Write(out,At(model,z,matching,step));out<<'}';
        }
        out<<"],\"operators\":[";std::vector<std::unique_ptr<Operator>> operators;
        for(const auto &mesh:std::vector<std::array<int,2>>{{3,2},{5,2},{7,2},{3,4},{3,8},{7,8}})
            for(double wave:{1.,2.,4.})operators.push_back(std::make_unique<Operator>(model,mesh[0],mesh[1],wave,matching));
        for(std::size_t i=0;i<operators.size();++i){if(i)out<<',';Describe(out,*operators[i]);}
        out<<"],\"trials\":[";count=0;
        for(const auto &op:operators)
        {if(count++)out<<',';Trial(out,model,*op,2,spatialDt,1e-4);}
        const auto &reference=*operators.back();
        for(int order:{1,2})for(double dt:{.0025,.00125,.000625})
        {if(count++)out<<',';Trial(out,model,reference,order,dt,1e-4);}
        for(double amplitude:{1e-5,1e-3})
        {if(count++)out<<',';Trial(out,model,reference,2,spatialDt,amplitude);}
        out<<"],\"restricted_controls\":[";count=0;
        for(double chi:{0.,.5,1.})for(double eta:{0.,1.})
        {
            options.requestedWeight=chi;const EquilibriumAWFProfiles control(geometry,K,model,{},options);
            const Operator candidate(model,7,8,4.,control.Reference().matchingVelocity,eta);
            double error=0;for(std::size_t i=0;i<candidate.load.size();++i)for(std::size_t j=0;j<candidate.load.size();++j)
                error=std::max(error,std::abs(candidate.load[i][j]-reference.load[i][j]));
            if(count++)out<<',';
            out<<"{\"chi\":"<<chi<<",\"eta\":"<<eta<<",\"matching_velocity\":"<<control.Reference().matchingVelocity
                <<",\"maximum_operator_difference\":"<<error<<'}';
        }
        out<<"],\"control_scope\":\"equal-input restricted outer tangent; chi uses original AWF momentum matching, eta subtracts fixed base stress; not a nonlinear whole-patch identity\""
            <<",\"checks\":"<<checks<<",\"success\":true}\n";
        Check(static_cast<bool>(out),"cannot finish report");std::cout<<"Reduced native transverse trials completed\n";
    }
    catch(const std::exception &error){std::cerr<<error.what()<<'\n';return 1;}
    return 0;
}
