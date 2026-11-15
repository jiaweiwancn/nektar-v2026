// SPDX-License-Identifier: MIT
#ifndef NEKTAR_RANS_WALL_FLOW_VECTOR_CAPTURE_H
#define NEKTAR_RANS_WALL_FLOW_VECTOR_CAPTURE_H
#include <IncNavierStokesSolver/EquationSystems/RANSWallMonitor.h>
#include <cmath>
#include <set>

namespace Nektar::RANSWallFlowCapture
{
// Passive serializer only: all vectors are increments already computed by the
// governing sweep. In particular bGz is NOT inferred from predictor-du.
inline std::string Write(const std::string &stem,double time,double alpha,
    const std::vector<double> &base,const std::vector<double> &flow,
    const std::vector<double> &history,const std::vector<int> &rows,std::size_t nv,
    const std::vector<double> &predictor,const std::vector<double> &bGz,
    const std::vector<double> &du,const std::vector<double> &dz,
    const std::string &executionIdentity,const std::string &summaries)
{
    namespace D=RANSWallR3;namespace J=RANSWallSolveJson;
    D::Check(sizeof(double)==8 && std::numeric_limits<double>::is_iec559,
        "flow vector capture requires IEEE binary64");
    const std::uint16_t endian=1;
    D::Check(*reinterpret_cast<const unsigned char*>(&endian)==1,"flow vector capture requires little endian");
    D::Check(nv>0 && rows.size()>nv && predictor.size()==nv && bGz.size()==nv &&
        du.size()==nv && dz.size()==rows.size()-nv && base.size()==flow.size() &&
        history.size()==base.size(),"invalid passive flow-vector shapes");
    std::set<int> seen;
    for(int id:rows)D::Check(id>=0 && std::size_t(id)<base.size() && seen.insert(id).second,
        "invalid passive native flow row map");
    D::Check(std::isfinite(time) && std::isfinite(alpha) && time>0 && alpha>0,
        "invalid passive flow time/alpha");
    const auto step=std::llround(time*alpha);
    D::Check(step>0 && std::abs(time-step/alpha)<=128*std::numeric_limits<double>::epsilon()*std::max(1.,std::abs(time)),
        "passive flow capture requires an aligned IMEX1 step");
    std::string bytes;bytes.reserve((3*nv+dz.size())*sizeof(double));
    std::vector<std::string> arrays;
    for(const auto &item:std::vector<std::pair<std::string,const std::vector<double>*>>{
        {"predictor",&predictor},{"bGz",&bGz},{"du",&du},{"dz",&dz}})
    {
        const auto offset=bytes.size();const auto &v=*item.second;
        bytes.append(reinterpret_cast<const char*>(v.data()),v.size()*sizeof(double));
        arrays.push_back(J::Object({{"name",J::Quote(item.first)},{"byte_offset",J::Number(offset)},
            {"count",J::Number(v.size())},{"bytes",J::Number(v.size()*sizeof(double))},
            {"sha256",J::Quote(D::Hash(v))}}));
    }
    const auto payload=stem+".vectors.bin",packet=stem+".json";
    RANSWallMonitor::Atomic(payload,bytes);
    RANSWallMonitor::Atomic(packet,J::Object({{"schema",J::Quote("WJMappedFlowVectorsV1")},
        {"status",J::Quote("FLOW_PROVISIONAL_NOT_TIME_COMMIT")},{"absolute_step",J::Number(step)},
        {"time",J::Number(time)},{"time_bits",J::Quote(RANSWallMonitor::Bits(time))},
        {"alpha_mass",J::Number(alpha)},{"alpha_bits",J::Quote(RANSWallMonitor::Bits(alpha))},
        {"base_sha256",J::Quote(D::Hash(base))},{"flow_state_sha256",J::Quote(D::Hash(flow))},
        {"callback_history_sha256",J::Quote(D::Hash(history))},
        {"attempt_id",J::Quote(D::Hash(base)+":"+D::Hash(history)+":"+RANSWallMonitor::Bits(time))},
        {"execution_identity",executionIdentity},{"payload",D::Pin(payload)},
        {"encoding",J::Quote("little-endian IEEE754 binary64; exact native increment vectors")},
        {"arrays",J::List(arrays)},{"velocity_rows",J::Numbers(std::vector<int>(rows.begin(),rows.begin()+nv))},
        {"pressure_rows",J::Numbers(std::vector<int>(rows.begin()+nv,rows.end()))},
        {"pressure_velocity_sign",J::Quote("actual bGz; signed correction is -bGz")},
        {"existing_residual_summaries",summaries},{"additional_governing_operations","0"}}));
    return D::Pin(packet);
}
}
#endif
