// SPDX-License-Identifier: MIT
#ifndef WJ_H26_CHECKS_H
#define WJ_H26_CHECKS_H
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include <LibUtilities/Foundations/ManagerAccess.h>
namespace H26
{
Eigen::Matrix3d Matrix(const Tensor&T){Eigen::Matrix3d A;for(int i=0;i<3;++i)for(int j=0;j<3;++j)A(i,j)=T[i][j];return A;}
Tensor TensorOf(const Eigen::Matrix3d&A){Tensor T{};for(int i=0;i<3;++i)for(int j=0;j<3;++j)T[i][j]=A(i,j);return T;}
void Check(bool condition,const char*message){if(!condition)throw std::runtime_error(message);}
void AdapterChecks(std::ostream&out)
{
    double k=.008809015145515192,e=.002315524054799254,nu=6.423173803526447e-5,nuA=.0007,H=.15;
    double worst=0,minEigen=1e100;int count=0;double start=Clock();
    auto error=[&](double v){worst=std::max(worst,v);Check(v<1e-10,"tensor adapter identity failed");};
    Eigen::Matrix3d rotation=Eigen::AngleAxisd(.63,Eigen::Vector3d(1,2,-1).normalized()).toRotationMatrix();
    for(double shear:{-.05,.02,.1})for(double r:{.01,.08,.12,.15,.16})
    {
        Tensor G{};G[0][2]=shear;auto w=Adapt(k,e,G,{0,0,1},nu,nuA,r,H);double q=k*shear/e;
        // Independent scalar cubic root, from the planar appendix (not WJ Root).
        double lo=1.8,hi=10.;for(int j=0;j<100;++j){double N=(lo+hi)/2,f=N*N*N-1.8*N*N-.35*q*q*N-1.8*q*q;if(f>0)hi=N;else lo=N;}double N=(lo+hi)/2,C=3*N/(5*(N*N+q*q)),d=C*q*q/N;
        error(std::abs(w.R[0][0]/k-(2./3+d)));error(std::abs(w.R[1][1]/k-2./3));error(std::abs(w.R[2][2]/k-(2./3-d)));
        error(std::abs(w.R[0][2]+((1-w.chi)*nuA+w.chi*C*k*k/e)*shear)/k);
        Eigen::Matrix3d R=Matrix(w.R);error(std::abs(R.trace()-2*k)/k);minEigen=std::min(minEigen,Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>(R).eigenvalues().minCoeff()/k);
        double work=-(R.array()*Matrix(G).array()).sum();error(std::abs(work-w.production)/e);
        auto rotatedG=TensorOf(rotation*Matrix(G)*rotation.transpose());Eigen::Vector3d rn=rotation*Eigen::Vector3d(0,0,1);auto rw=Adapt(k,e,rotatedG,{rn[0],rn[1],rn[2]},nu,nuA,r,H);
        error((Matrix(rw.R)-rotation*R*rotation.transpose()).norm()/k);
        error((Matrix(w.V)-w.mu*(Matrix(G)+Matrix(G).transpose())-Matrix(w.remainder)).norm()/k);
        if(r>=.12){auto full=SolverUtils::EvaluateWJClosure(k,e,G,{},SolverUtils::WJGradientPolicy::DevStrainRawWorkV1,0);error((R-Matrix(full.reynoldsStress)).norm()/k);}
        ++count;
    }
    Tensor G{{{.012,.004,.03},{-.002,-.003,.001},{-.01,.002,-.008}}};Vector n{.3,.4,std::sqrt(.75)};auto inside=Adapt(k,e,G,n,nu,nuA,.14,H),outside=Adapt(k,e,G,n,nu,nuA,H+.01,H);
    error((Matrix(inside.R)-Matrix(outside.R)).norm()/k);double p=.02,Pi=p+2*k/3;
    Eigen::Vector3d rotatedN=rotation*Eigen::Vector3d(n[0],n[1],n[2]);auto generalRotated=Adapt(k,e,TensorOf(rotation*Matrix(G)*rotation.transpose()),{rotatedN[0],rotatedN[1],rotatedN[2]},nu,nuA,.14,H);error((Matrix(generalRotated.R)-rotation*Matrix(inside.R)*rotation.transpose()).norm()/k);error(std::abs(Matrix(inside.R).trace()-2*k)/k);error(std::abs(inside.production+(Matrix(inside.R).array()*Matrix(G).array()).sum())/e);minEigen=std::min(minEigen,Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>(Matrix(inside.R)).eigenvalues().minCoeff()/k);
    Eigen::Vector3d normal(n[0],n[1],n[2]);error(((Matrix(inside.V)-Pi*Eigen::Matrix3d::Identity())*normal-(nu*(Matrix(G)+Matrix(G).transpose())-Matrix(inside.R)-p*Eigen::Matrix3d::Identity())*normal).norm()/k);
    Check(minEigen>=0,"selected adapter fixture is unrealizable; no clipping allowed");
    out<<"{\"status\":\"PASS\",\"fixtures\":"<<count<<",\"maximum_dimensionless_identity_error\":"<<worst<<",\"minimum_R_eigenvalue_over_k\":"<<minEigen<<",\"seconds\":"<<Clock()-start<<",\"checks\":\"independent planar cubic/C/d, trace, work, covariance, full-WJ collar, positive realizability, finite-nu total traction\"}";
}
double TraceValue(const Element&e,const Vec&physical,double xi,double eta)
{
    Array<OneD,NekDouble>point(2),values(physical.size());point[0]=xi;point[1]=eta;for(int q=0;q<physical.size();++q)values[q]=physical[q];return e.velocity->StdPhysEvaluate(point,values);
}
void WallTrace(Flow&flow,std::ostream&out)
{
    // Smooth prescribed tensor fixture: project 0.01 times a trigonometric
    // velocity to native Q3. Gradients retain distinct element limits at H.
    MMS m;m.solution="trigonometric";m.geo=flow.Mapping();m.xl=flow.Elements().front().xl;m.Lx=flow.Elements().back().xr-m.xl;m.zl=m.geo.d0;m.Lz=flow.Elements().back().rr;Vec state=Exact(flow,m.Functions(0));state.head(2*flow.VelocitySize())*=.01;
    double k=.008809015145515192,e=.002315524054799254,nu=6.423173803526447e-5,jump=0,area=0,pressureMismatch=0,minEigen=1e100;Eigen::Matrix<double,2,3>moments=Eigen::Matrix<double,2,3>::Zero();int faces=0;
    // Native trace interpolation, with an independent Gauss face rule supplied
    // by StdRegions' public Points manager (not volume endpoint weights).
    LibUtilities::PointsKey key(20,LibUtilities::eGaussGaussLegendre);auto points=LibUtilities::PointsManager()[key];auto nodes=points->GetZ(),weights=points->GetW();
    auto trace=[&](const Element&el,double xi,double eta){Vec u(16),v(16),p(9);for(int j=0;j<16;++j){u[j]=state[el.u[j]];v[j]=state[flow.VelocitySize()+el.u[j]];}for(int j=0;j<9;++j)p[j]=state[2*flow.VelocitySize()+el.p[j]];Tensor G{};G[0][0]=TraceValue(el,el.Dx*u,xi,eta);G[0][2]=TraceValue(el,el.Dz*u,xi,eta);G[2][0]=TraceValue(el,el.Dx*v,xi,eta);G[2][2]=TraceValue(el,el.Dz*v,xi,eta);return std::make_pair(G,TraceValue(el,el.P*p,xi,eta));};
    for(const auto&lower:flow.Elements())if(std::abs(lower.rr-m.geo.H)<1e-10)
    {
        const Element*upper=nullptr;for(const auto&el:flow.Elements())if(std::abs(el.rl-m.geo.H)<1e-10&&std::abs(el.xl-lower.xl)<1e-10&&std::abs(el.xr-lower.xr)<1e-10){upper=&el;break;}Check(upper,"missing two-sided native matching-height trace");++faces;
        for(int q=0;q<20;++q){double xi=nodes[q],x=(lower.xl+lower.xr)/2+(lower.xr-lower.xl)*xi/2,s=m.geo.Slope(x),norm=std::sqrt(1+s*s),w=weights[q]*(lower.xr-lower.xl)*norm/2;Vector normal{-s/norm,0,1/norm};auto a=trace(lower,xi,1),b=trace(*upper,xi,-1);auto Ra=Adapt(k,e,a.first,normal,nu,.0007,m.geo.H*(1-1e-12),m.geo.H),Rb=Adapt(k,e,b.first,normal,nu,.0007,m.geo.H*(1+1e-12),m.geo.H);pressureMismatch=std::max(pressureMismatch,std::abs(a.second-b.second));minEigen=std::min({minEigen,Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>(Matrix(Ra.R)).eigenvalues().minCoeff()/k,Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>(Matrix(Rb.R)).eigenvalues().minCoeff()/k});
            for(int i=0;i<2;++i){int c=i?2:0;double delta=0;for(int j=0;j<3;++j)delta+=(nu*(a.first[c][j]+a.first[j][c]-b.first[c][j]-b.first[j][c])-Ra.R[c][j]+Rb.R[c][j]-(c==j?a.second-b.second:0))*normal[j];jump+=w*delta*delta;for(int j=0;j<3;++j)moments(i,j)+=w*delta*std::pow((x-m.xl)/m.Lx,j);}area+=w;
        }
    }Check(faces>0&&pressureMismatch<1e-11&&minEigen>=0,"C0 WJ trace fixture failed; no clipping allowed");out<<"{\"status\":\"PASS_FIXTURE_REQUIRES_REFINEMENT_COMPARISON\",\"scope\":\"native C0 projection of smooth prescribed velocity, common positive k/E and full-WJ collar, distinct element gradients; no C1 constraint\",\"faces\":"<<faces<<",\"L2_traction_jump_over_k\":"<<std::sqrt(jump/area)/k<<",\"weak_traction_moments_over_k\":"<<moments.norm()/area/k<<",\"pressure_trace_mismatch\":"<<pressureMismatch<<",\"minimum_R_eigenvalue_over_k\":"<<minEigen<<'}';
}
struct Constitutive
{
    std::vector<std::vector<WallStress>>values;double low=1e100,high=0,seconds=0;
    Constitutive(Flow&flow,const Vec&u)
    {
        double began=Clock();auto geo=flow.Mapping();int n=flow.VelocitySize();
        for(const auto&el:flow.Elements())
        {
            values.emplace_back();Vec a(16),b(16);for(int i=0;i<16;++i){a[i]=u[el.u[i]];b[i]=u[n+el.u[i]];}Vec ux=el.Dx*a,uz=el.Dz*a,vx=el.Dx*b,vz=el.Dz*b;
            for(int q=0;q<el.x.size();++q){double x=el.x[q],z=el.z[q],s=geo.Slope(x),r=z-geo.Ground(x);Tensor G{};G[0][0]=ux[q];G[0][2]=uz[q];G[2][0]=vx[q];G[2][2]=vz[q];Vector normal{-s/std::sqrt(1+s*s),0,1/std::sqrt(1+s*s)};auto w=Adapt(.008809015145515192,.002315524054799254,G,normal,6.423173803526447e-5,.0007*(1+r),r,geo.H);low=std::min(low,w.mu);high=std::max(high,w.mu);values.back().push_back(w);}
        }seconds=Clock()-began;
    }
    std::vector<Vec> Mu()const{std::vector<Vec>result;for(const auto&el:values){Vec v(el.size());for(int q=0;q<v.size();++q)v[q]=el[q].mu;result.push_back(v);}return result;}
    std::vector<std::vector<Tensor>> Stress(bool complete=false)const{std::vector<std::vector<Tensor>>result;for(const auto&el:values){result.emplace_back();for(const auto&w:el)result.back().push_back(complete?w.V:w.remainder);}return result;}
};
void Nonlinear(Flow&flow,std::ostream&out)
{
    double began=Clock(),alpha=300,awf=1,constitutive=0;Vec target=Exact(flow,Polynomial(0)),lift=flow.Lift(target);Constitutive reference(flow,target);constitutive+=reference.seconds;
    Sparse Kr=flow.StrainSamples(reference.Mu());Vec sr=flow.StressSamples(reference.Stress());Vec body=flow.ApplyFull(alpha,Kr,awf,target)-sr;
    Vec current=.8*target;current.tail(flow.PressureSize()).setZero();for(int c=0;c<2;++c)current.segment(c*flow.VelocitySize(),flow.DirichletSize())=lift.segment(c*flow.VelocitySize(),flow.DirichletSize());
    std::vector<SolveInfo>linear;std::vector<double>residuals;double moment=1,continuity=1;int outer=0;double weakBlend=0,low=0,high=0;
    for(;outer<60;++outer)
    {
        Constitutive data(flow,current);constitutive+=data.seconds;Sparse K=flow.StrainSamples(data.Mu());Vec load=body+flow.StressSamples(data.Stress());Sparse A=flow.System(alpha,K,awf);Vec rhs=flow.LiftedRhs(alpha,K,awf,load,lift),x=flow.Restrict(current);linear.push_back(flow.Solve(A,rhs,x,0,1e-12));current=flow.Expand(x,lift);
        Constitutive fresh(flow,current);constitutive+=fresh.seconds;Sparse latestK=flow.StrainSamples(fresh.Mu());Vec latestStress=flow.StressSamples(fresh.Stress());Sparse latestA=flow.System(alpha,latestK,awf);Vec latestRhs=flow.LiftedRhs(alpha,latestK,awf,body+latestStress,lift);auto residual=flow.Residual(latestA,latestRhs,x);moment=residual.momentum;continuity=residual.continuity;residuals.push_back(moment);low=fresh.low;high=fresh.high;
        Vec all=Vec::Zero(current.size());all.head(2*flow.VelocitySize())=latestK*current.head(2*flow.VelocitySize());all-=latestStress;all+=flow.StressSamples(fresh.Stress(true));weakBlend=flow.VelocityDual(all);
        if(std::max(moment,continuity)<1e-10)break;
    }
    Check(outer<60&&moment<=1e-9&&continuity<=1e-9,"nonlinear frozen-scalar flow failed");Check(weakBlend<1e-10,"complete blended stress weak identity failed");double error=(current-target).norm()/target.norm();Check(error<1e-8,"nonlinear flow direct target agreement failed");
    out<<"{\"status\":\"PASS\",\"scope\":\"flow-only algebraic nonlinear refresh with prescribed positive k and epsilon, full adapter and current AWF pressure; not scalar transport\",\"outer_iterations\":"<<outer+1<<",\"momentum\":"<<moment<<",\"continuity\":"<<continuity<<",\"solution_error\":"<<error<<",\"complete_blend_weak_identity\":"<<weakBlend<<",\"mu_range\":["<<low<<','<<high<<"],\"constitutive_seconds\":"<<constitutive<<",\"total_seconds\":"<<Clock()-began<<",\"residuals\":[";for(int i=0;i<residuals.size();++i){if(i)out<<',';out<<residuals[i];}out<<"],\"linear\":[";for(int i=0;i<linear.size();++i){if(i)out<<',';InfoJSON(out,linear[i]);}out<<"]}";
}
}
#endif
