// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/RANSWallFlowVectorCapture.h>
#include <iostream>
int main(int argc,char **argv)
{
    using namespace Nektar;
    try
    {
        RANSWallR3::Check(argc==2,"fixture requires an isolated output directory");
        const std::string stem=std::string(argv[1])+"/flow";
        const std::vector<double> base{0.,-0.,2.,3.,4.,5.},flow{.5,-.5,2.,3.,4.25,5.};
        const std::vector<int> rows{0,1,4};
        const std::vector<double> predictor{.75,-0.},bGz{.25,.5},du{.5,-.5},dz{.25};
        const auto result=RANSWallFlowCapture::Write(stem,.005,200.,base,flow,base,rows,2,
            predictor,bGz,du,dz,"{\"fixture\":true}","{\"unchanged\":true}");
        auto packet=RANSWallMonitor::Read(stem+".json");
        RANSWallR3::Check(packet.get<std::string>("schema")=="WJMappedFlowVectorsV1","wrong schema");
        std::ifstream f(stem+".vectors.bin",std::ios::binary);
        std::vector<double> got(7);f.read(reinterpret_cast<char*>(got.data()),56);
        const std::vector<double> expected{.75,-0.,.25,.5,.5,-.5,.25};
        RANSWallR3::Check(f.gcount()==56 && RANSWallR3::Hash(got)==RANSWallR3::Hash(expected),"increment/zero bytes differ");
        bool refused=false;
        try{RANSWallFlowCapture::Write(stem,.005,200.,base,flow,base,rows,2,predictor,bGz,du,dz,"{}","{}");}
        catch(const std::exception &){refused=true;}
        RANSWallR3::Check(refused,"immutable overwrite accepted");
        refused=false;
        try{RANSWallFlowCapture::Write(stem+"-bad",.005,200.,base,flow,base,{0,0,4},2,predictor,bGz,du,dz,"{}","{}");}
        catch(const std::exception &){refused=true;}
        RANSWallR3::Check(refused,"duplicate row accepted");
        std::cout<<"flow vector serializer/sign/zero/row/immutability checks PASS; zero governing work "<<result<<std::endl;
    }
    catch(const std::exception &e){std::cerr<<e.what()<<std::endl;return 1;}
}
