// SPDX-License-Identifier: MIT
#ifndef NEKTAR_WJ_NATIVE_INITIAL_FIELDS_H
#define NEKTAR_WJ_NATIVE_INITIAL_FIELDS_H

#include <Eigen/Core>
#include <array>
#include <map>
#include <string>
#include <vector>

namespace Nektar::WJNative
{
enum class InitialFieldRole { QualifiedReference, Candidate };

struct InitialFieldGeometry
{
    double d0=0,H=0,h=0,left=0,right=0,top=0,curve=0;
    std::string normalMap;
};

struct InitialFieldMetadata
{
    std::string schema,stateSchema,model,sourceIdentity,referenceIdentity;
    std::string basis,mmsVariant,packetSha256,stateSha256,outletPolicy;
    bool manufactured=false,gauge=false,qualifiedReference=false;
    InitialFieldRole role=InitialFieldRole::Candidate;
};

struct InitialFieldElement
{
    int geometry=-1;
    double xl=0,xr=0,rl=0,rr=0;
    std::vector<int>u,p,k,epsilon;
    std::vector<double>us,ps,ks,epsilonSigns;
};

/** Immutable step-zero coefficient packet, never a cold-history restart.
 * Values use the original signed Modified_A Q3/Q2 maps. An explicit owner is
 * preserved. Without one, the lowest geometry ID containing the point and
 * owning the requested field is chosen; epsilon at H uses an outer owner.
 * Only geometric roundoff is tolerated. No scalar, coefficient or state is
 * clipped. Projection, current Gamma reconstruction, fresh advection/history,
 * admissibility, and compatibility policy remain responsibilities of callers.
 */
class InitialFields
{
public:
    explicit InitialFields(const std::string&jsonPath);
    const InitialFieldMetadata&Metadata()const{return m_metadata;}
    const InitialFieldGeometry&GeometryData()const{return m_geometry;}
    // Native packet order: Nu, Np, Nk, Ne (epsilon has only outer coefficients).
    const std::array<int,4>&Sizes()const{return m_sizes;}
    const Eigen::VectorXd&Flow()const{return m_flow;}
    const Eigen::VectorXd&K()const{return m_k;}
    const Eigen::VectorXd&Epsilon()const{return m_epsilon;}
    const std::vector<InitialFieldElement>&Elements()const{return m_elements;}
    // 0=u, 1=v (canonical signed u map), 2=k, 3=outer epsilon, 4=Pi.
    double Evaluate(unsigned fieldIndex,double x,double physicalZ,int ownerGeometry=-1)const;
    void ValidateBindings(const std::string&expectedSource,const std::string&expectedReference,
        const std::string&expectedModel,bool expectedManufactured,const std::string&expectedVariant,
        const InitialFieldGeometry&expectedGeometry)const;
    // Geometry IDs, signed maps and vector dimensions must match. Bounds may
    // differ only by the declared 64-epsilon geometric coincidence tolerance.
    bool SameInventory(const std::vector<InitialFieldElement>&elements,
                       const std::array<int,4>&sizes)const;

private:
    double Value(const InitialFieldElement&element,unsigned field,double x,double r)const;
    bool Contains(const InitialFieldElement&element,double x,double r)const;
    void ValidateInventory();
    void ValidateContinuity()const;
    InitialFieldMetadata m_metadata;
    InitialFieldGeometry m_geometry;
    std::array<int,4>m_sizes{};
    Eigen::VectorXd m_flow,m_k,m_epsilon;
    std::vector<InitialFieldElement>m_elements;
    std::map<int,std::size_t>m_index;
    std::vector<double>m_xBreaks;
    std::vector<std::vector<std::size_t>>m_slabs;
};
}
#endif
