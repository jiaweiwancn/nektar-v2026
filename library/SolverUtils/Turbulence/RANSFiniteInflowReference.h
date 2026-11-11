// SPDX-License-Identifier: MIT
#ifndef NEKTAR_SOLVERUTILS_RANS_FINITE_INFLOW_REFERENCE_H
#define NEKTAR_SOLVERUTILS_RANS_FINITE_INFLOW_REFERENCE_H
#include <SolverUtils/Turbulence/RANSRobustVariant.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace Nektar::SolverUtils
{
/** Immutable, prescribed incoming-frame reference, independent of solution,
 * wall distance and chart. SI table columns: zeta,U,dU/dzeta,k,epsilon.
 * U uses C1 cubic Hermite interpolation and its exact derivative. k/epsilon
 * use positive log-linear interpolation. Beyond the final zero-slope knot all
 * fields are constant. No moving-flow fitting, clipping or isotropy is used.
 * Content bytes and interpretation belong in the caller's physics signature.
 */
class RANSFiniteInflowReference
{
public:
    struct Knot { double z,U,dU,k,epsilon; };
    struct Value
    {
        double zeta=0,U=0,dU=0,k=0,epsilon=0,nuT=0,gamma=0;
        RANSTensor strain{};
    };
    RANSFiniteInflowReference(std::vector<Knot> knots,RANSVector origin,
        RANSVector flow,RANSVector normal,double scale,
        const RANSParameters &parameters={})
        : m_knots(std::move(knots)),m_origin(origin),m_flow(flow),
          m_normal(normal),m_scale(scale),m_parameters(parameters)
    {
        Check(m_knots.size()>=2,"at least two finite-reference knots required");
        double ff=0,nn=0,fn=0;
        for(int a=0;a<3;++a)
        {
            Check(std::isfinite(origin[a]) && std::isfinite(flow[a]) && std::isfinite(normal[a]),"nonfinite incoming frame");
            ff+=flow[a]*flow[a];nn+=normal[a]*normal[a];fn+=flow[a]*normal[a];
        }
        Check(std::abs(ff-1)<1e-13 && std::abs(nn-1)<1e-13 && std::abs(fn)<1e-13,"incoming frame must be orthonormal");
        Check(std::isfinite(scale) && scale>=0,"invalid prescribed gamma scale");
        Check(m_knots.front().z==0 && m_knots.back().dU==0,"finite reference requires zeta=0 origin and a zero-shear final knot");
        double last=-1;
        for(const auto &q:m_knots)
        {
            Check(std::isfinite(q.z) && std::isfinite(q.U) && std::isfinite(q.dU) && std::isfinite(q.k) && std::isfinite(q.epsilon),"nonfinite reference knot");
            Check(q.z>last && q.k>0 && q.epsilon>0,"ordered knots and positive prescribed k/epsilon required");last=q.z;
        }
        m_parameters.Validate();
    }
    static std::vector<Knot> Parse(const std::string &text)
    {
        std::istringstream in(text);in.imbue(std::locale::classic());
        std::string line;std::getline(in,line);
        Check(line=="FiniteInflowV1 SI:z,U,dU_dz,k,epsilon","unknown finite-reference format/units");
        std::vector<Knot> knots;
        while(std::getline(in,line))
        {
            if(line.empty())continue;
            std::istringstream row(line);row.imbue(std::locale::classic());Knot q{};std::string extra;
            Check(bool(row>>q.z>>q.U>>q.dU>>q.k>>q.epsilon) && !(row>>extra),"invalid finite-reference row");knots.push_back(q);
        }
        return knots;
    }
    Value At(const RANSVector &physical) const
    {
        return AtHeight(Height(physical));
    }
    double Height(const RANSVector &physical) const
    {double z=0;for(int a=0;a<3;++a)z+=(physical[a]-m_origin[a])*m_normal[a];return z;}
    Value AtHeight(double z) const
    {
        Check(std::isfinite(z) && z>=0,"incoming-frame height outside finite reference");
        Value v;v.zeta=z;
        if(z>=m_knots.back().z)
        {const auto &q=m_knots.back();v.U=q.U;v.dU=0;v.k=q.k;v.epsilon=q.epsilon;}
        else
        {
            auto b=std::upper_bound(m_knots.begin(),m_knots.end(),z,[](double x,const Knot &q){return x<q.z;});
            const auto &a=*(b-1);const double h=b->z-a.z,t=(z-a.z)/h,t2=t*t,t3=t2*t;
            v.U=(2*t3-3*t2+1)*a.U+(t3-2*t2+t)*h*a.dU+(-2*t3+3*t2)*b->U+(t3-t2)*h*b->dU;
            v.dU=((6*t2-6*t)*a.U+(-6*t2+6*t)*b->U)/h+(3*t2-4*t+1)*a.dU+(3*t2-2*t)*b->dU;
            v.k=std::exp((1-t)*std::log(a.k)+t*std::log(b->k));
            v.epsilon=std::exp((1-t)*std::log(a.epsilon)+t*std::log(b->epsilon));
        }
        RANSTensor g{};
        for(int a=0;a<3;++a)for(int b=0;b<3;++b)
        {g[a][b]=v.dU*m_flow[a]*m_normal[b];v.strain[a][b]=v.dU*(m_flow[a]*m_normal[b]+m_normal[a]*m_flow[b]);}
        // Eddy viscosity has no Yap-distance dependence. A positive dummy
        // distance permits zeta=0 without interpreting the incoming height as
        // physical wall distance. No other closure outputs are consumed here.
        v.nuT=EvaluateRANSClosure(RANSModel::CLSHighRe,v.k,v.epsilon,g,1.,m_parameters,0).eddyViscosity;
        v.gamma=m_scale*v.nuT;
        Check(std::isfinite(v.gamma) && v.gamma>=0,"invalid immutable reference viscosity");return v;
    }
    std::string Interpretation() const
    {
        std::ostringstream s;s.imbue(std::locale::classic());s<<std::setprecision(17)
          <<"FiniteInflowV1:SI:HermiteU-exactDerivative:loglinear-k-epsilon:constant-final-zero-shear:gamma=Gamma*OriginalCLSnuT:gamma0=0:Gamma="<<m_scale;
        for(const auto &v:{m_origin,m_flow,m_normal})for(double x:v)s<<':'<<x;
        return s.str();
    }
    const std::vector<Knot> &Knots() const { return m_knots; }
private:
    static void Check(bool yes,const char *message)
    {if(!yes)throw std::invalid_argument(std::string("Finite inflow reference: ")+message);}
    const std::vector<Knot> m_knots;
    const RANSVector m_origin,m_flow,m_normal;
    const double m_scale;
    const RANSParameters m_parameters;
};
}
#endif
