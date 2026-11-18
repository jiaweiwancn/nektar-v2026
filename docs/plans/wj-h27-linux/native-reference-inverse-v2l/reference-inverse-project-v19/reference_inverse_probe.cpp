// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/WJNativeReference.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <boost/property_tree/json_parser.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>

using Nektar::WJNative::Reference;
namespace J=Nektar::RANSWallSolveJson;
namespace pt=boost::property_tree;

int main(int argc,char **argv)
{
    try
    {
        if(argc!=4)throw std::invalid_argument("usage: reference_inverse_probe table.json mutations.json new-output.json");
        if(std::filesystem::exists(argv[3]))throw std::invalid_argument("immutable probe output exists");
        Reference reference(argv[1]);const auto &p=reference.Parameters();
        const auto range=reference.CutoffRadiusRange(),bounds=reference.LogCutoffDerivativeBounds();
        std::vector<std::string> roundtrips,inverses,rejections,mutations;
        for(int i=0;i<=128;++i)
        {
            const double K=i==0?p.KMinimum:i==128?p.KMaximum:p.KMinimum+(p.KMaximum-p.KMinimum)*i/128;
            const double radius=i==0?range[0]:i==128?range[1]:range[0]+(range[1]-range[0])*i/128;
            const double c=reference.Cutoff(K);
            roundtrips.push_back(J::Object({{"K",J::Number(K)},{"c",J::Number(c)},
                {"inverse",J::Number(reference.KAtCutoff(c))},{"dc",J::Number(reference.CutoffDerivative(K))}}));
            inverses.push_back(J::Object({{"r",J::Number(radius)},{"K",J::Number(reference.KAtCutoff(radius))}}));
        }
        auto reject=[&](const std::string &name,const auto &operation,const std::string &expected)
        {
            std::string reason;bool rejected=false;
            try{operation();}catch(const std::exception &error){rejected=true;reason=error.what();}
            if(!rejected||reason.find(expected)==std::string::npos)
                throw std::runtime_error("missing/wrong inverse rejection: "+name+"; "+reason);
            return J::Object({{"name",J::Quote(name)},{"reason",J::Quote(reason)}});
        };
        for(const auto &item:std::vector<std::pair<std::string,double>>{
            {"below inverse interval",std::nextafter(range[0],0.)},
            {"above inverse interval",std::nextafter(range[1],std::numeric_limits<double>::infinity())},
            {"negative radius",-1e-20},{"zero radius",0.}})
            rejections.push_back(reject(item.first,[&]{reference.KAtCutoff(item.second);},"outside qualified inverse interval"));
        for(const auto &item:std::vector<std::pair<std::string,double>>{
            {"NaN radius",std::numeric_limits<double>::quiet_NaN()},
            {"positive infinity radius",std::numeric_limits<double>::infinity()},
            {"negative infinity radius",-std::numeric_limits<double>::infinity()}})
            rejections.push_back(reject(item.first,[&]{reference.KAtCutoff(item.second);},"nonfinite"));
        pt::ptree manifest;pt::read_json(argv[2],manifest);
        if(manifest.get<std::string>("schema")!="WJH27ReferenceInverseMutationsV1")throw std::invalid_argument("wrong mutation schema");
        for(const auto &entry:manifest.get_child("mutations"))
        {
            const auto &m=entry.second;const auto path=m.get<std::string>("file.path"),sha=m.get<std::string>("file.sha256");
            const auto name=m.get<std::string>("name");
            // Only the diagnostic's explicit hash admits synthetic bytes far
            // enough to exercise the native constructor guard. Applications
            // retain their default immutable qualified-table identity.
            mutations.push_back(J::Object({{"name",J::Quote(name)},
                {"default_identity",reject(name,[&]{Reference invalid(path);},"immutable table SHA256 differs")},
                {"constructor_guard",reject(name,[&]{Reference invalid(path,sha);},m.get<std::string>("expected_reason"))}}));
        }
        const double a=reference.KAtCutoff(range[0]),b=reference.KAtCutoff(range[1]);
        const double again=reference.KAtCutoff(range[0]);
        if(a!=p.KMaximum||b!=p.KMinimum||a!=again)throw std::runtime_error("endpoint/query-order mismatch");
        const auto result=J::Object({{"schema",J::Quote("WJH27NativeReferenceInverseProbeV1")},
            {"status",J::Quote("RECORDED")},{"source_identity",J::Quote(WJ_H27_SOURCE_ID)},
            {"table_sha256",J::Quote(reference.Identity())},{"radius_range",J::Numbers(range)},
            {"log_cutoff_slope_bounds",J::Numbers(bounds)},{"roundtrips",J::List(roundtrips)},
            {"inverses",J::List(inverses)},{"rejections",J::List(rejections)},
            {"mutations",J::List(mutations)},{"exact_endpoints_and_return_query",J::Boolean(true)},
            {"full_H27_qualification",J::Boolean(false)},{"TASK3_READY",J::Boolean(false)}});
        std::ofstream output(argv[3]);output<<result<<'\n';if(!output.good())throw std::runtime_error("write failed");
        std::cout<<"SOURCE19_REFERENCE_INVERSE_RECORDED\n";
    }
    catch(const std::exception &error){std::cerr<<error.what()<<'\n';return 2;}
}
