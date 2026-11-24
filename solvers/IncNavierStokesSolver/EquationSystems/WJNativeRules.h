// SPDX-License-Identifier: MIT
#ifndef NEKTAR_WJ_NATIVE_RULES_H
#define NEKTAR_WJ_NATIVE_RULES_H
#include <IncNavierStokesSolver/EquationSystems/WJNativeFlow.h>
namespace Nektar::WJNative
{
struct NativeRule
{
    int geometry=0;
    std::vector<int>u,p;
    std::vector<double>us,ps;
    Vec x,z,w; // Chart labels (a, d0+rho); w uses the physical measure.
    Mat B,Dx,Dz,P,Px,Pz;
    Vector normal{};
    std::vector<Vector> outward, wallNormals;
    Vec matchingCurvature;
    Vector Outward(int q) const { return outward.empty() ? normal : outward[q]; }
    Vec NormalVelocity(const Vec &u, const Vec &v) const
    {
        if (outward.empty())
            return normal[0]*u + normal[2]*v;
        Vec result(u.size());
        for (int q=0; q<u.size(); ++q)
            result[q] = outward[q][0]*u[q] + outward[q][2]*v[q];
        return result;
    }
    Vector WallNormal(int q) const
    {
        return wallNormals.empty() ? Vector{0,0,1} : wallNormals[q];
    }
    double MatchingCurvature(int q) const
    {
        return matchingCurvature.size() ? matchingCurvature[q] : 0.;
    }
};
using Rules=std::vector<std::shared_ptr<const NativeRule>>;
struct CutoffPartition
{
    std::array<double,2>radiusRange;
    std::function<std::vector<double>(const Element&)>xBreaks;
};
/** Native tensor basis factors interpolated from Nektar B/Db data to split
 * Gauss rules in affine chart cells, with physical weights and gradients.
 */
class RuleFactory
{
public:
    RuleFactory(const Flow&flow,int nx=8,int nz=12);
    const Rules&Static()const{return m_static;}
    Rules Cutoff(const std::function<double(const Element&,double)>&cutoff,
        const CutoffPartition*partition=nullptr)const;
    NativeRule Face(const Element&element,int side,int order=16,
        const std::function<double(const Element&,double)>&cutoff={})const;
    NativeRule Points(const Element&element,const Vec&x,const Vec&z)const;
    Vector FaceNormal(const Element &element, int side, double a, double rho) const;
    /** One-sided horizontal native edge. Its reference height is exactly
     * -1/+1; physical subtraction must not turn a material face into an
     * interior point on a thin cell. side is 2 (lower) or 3 (upper).
     */
    NativeRule HorizontalTrace(const Element&element,const Vec&x,int side)const;
private:
    NativeRule Volume(const Element&element,const std::function<double(const Element&,double)>&cutoff,
        const std::vector<double>&xBreaks={})const;
    NativeRule Evaluate(const Element&element,const Vec&xi,const Vec&eta,const Vec&weights)const;
    const Flow&m_flow;
    int m_nx,m_nz;
    Rules m_static;
};
}
#endif
