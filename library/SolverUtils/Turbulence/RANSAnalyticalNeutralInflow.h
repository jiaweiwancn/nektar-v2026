// SPDX-License-Identifier: MIT
#ifndef NEKTAR_RANS_ANALYTICAL_NEUTRAL_INFLOW_H
#define NEKTAR_RANS_ANALYTICAL_NEUTRAL_INFLOW_H
#include <SolverUtils/Turbulence/RANSFiniteInflowReference.h>
namespace Nektar::SolverUtils
{
/** Exact, immutable incoming-frame law. No displaced height, wall-distance
 * substitution, table interpolation or freestream completion. The caller owns
 * the inner AWF completion below H; this law is defined only at positive zeta. */
class RANSAnalyticalNeutralInflow
{
public:
    using Value=RANSFiniteInflowReference::Value;
    RANSAnalyticalNeutralInflow(const RANSNeutralReference &model,double us,
        double H,double z0,double Gamma,RANSVector origin,RANSVector flow,RANSVector normal)
        : m_model(model),m_us(us),m_H(H),m_z0(z0),m_Gamma(Gamma),
          m_origin(origin),m_flow(flow),m_normal(normal)
    {
        if(!std::isfinite(us+H+z0+Gamma) || model.Model()!=RANSModel::CLSHighRe || !(us>0 && H>z0 && z0>0 && Gamma>=0))
            throw std::invalid_argument("AnalyticalNeutralInflowV1 invalid CLS reference");
        double ff=0,nn=0,fn=0;
        for(int i=0;i<3;++i)
        {if(!std::isfinite(origin[i]+flow[i]+normal[i]))throw std::invalid_argument("nonfinite neutral frame");
         ff+=flow[i]*flow[i];nn+=normal[i]*normal[i];fn+=flow[i]*normal[i];}
        if(std::abs(ff-1)>1e-13 || std::abs(nn-1)>1e-13 || std::abs(fn)>1e-13)
            throw std::invalid_argument("neutral frame is not orthonormal");
        (void)AtHeight(H);
    }
    double Height(const RANSVector &p) const
    {double z=0;for(int i=0;i<3;++i)z+=(p[i]-m_origin[i])*m_normal[i];return z;}
    Value At(const RANSVector &p) const {return AtHeight(Height(p));}
    Value AtHeight(double z) const
    {
        const auto n=m_model.AtHeight(z,m_us,m_H,m_us/m_model.Kappa()*std::log(m_H/m_z0));
        Value v;v.zeta=z;v.U=n.velocity;v.dU=n.velocityGradient;v.k=n.k;
        v.epsilon=n.epsilon;v.nuT=n.eddyViscosity;v.gamma=m_Gamma*v.nuT;
        for(int i=0;i<3;++i)for(int j=0;j<3;++j)
            v.strain[i][j]=v.dU*(m_flow[i]*m_normal[j]+m_normal[i]*m_flow[j]);
        return v;
    }
    std::string Interpretation() const
    {
        std::ostringstream s;s.imbue(std::locale::classic());s<<std::setprecision(17)
          <<"AnalyticalNeutralInflowV1:SI:exact-log-unshifted-global-height:exact-1/z-epsilon:gamma=Gamma*kappa_CLS*u_tau*z:no-tail:no-displacement:";
        for(double x:{m_us,m_H,m_z0,m_Gamma,m_model.Q(),m_model.Kappa()})s<<':'<<x;
        for(const auto &v:{m_origin,m_flow,m_normal})for(double x:v)s<<':'<<x;
        return s.str();
    }
private:
    const RANSNeutralReference m_model;
    const double m_us,m_H,m_z0,m_Gamma;
    const RANSVector m_origin,m_flow,m_normal;
};
}
#endif
