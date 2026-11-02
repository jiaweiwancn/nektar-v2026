///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
using namespace Nektar;
int main(int argc,char **argv)
{
    try
    {
        if(argc!=3)throw std::invalid_argument("usage: RANSWallExpressionProbe SESSION.xml OUTPUT.json");
        auto session=LibUtilities::SessionReader::CreateInstance(2,argv);session->InitSession();
        int dimension=2;session->LoadParameter("RANSWallExpressionDimension",dimension,2);
        std::ofstream out(argv[2]);out<<std::setprecision(17)<<"{\"schema\":1,\"samples\":[";int count=0;
        for(double x:{-3.5,-2.8,-2.4,-1.5,0.5})for(double r:{0.,1e-5,1e-4,.01,.5,1.,2.,3.})
        {
            if(count++){out<<',';}
            out<<"{\"point\":["<<x<<','<<r<<",0.3],\"values\":[";int entry=0;
            for(const auto &name:{"InitialConditionsInner","InitialConditionsOuter","ScalarSourceInner","ScalarSourceOuter"})
                for(const auto &var:{"k","eps_tilde"})if(session->DefinesFunction(name,var))
                {
                    if(entry++){out<<',';}
                    const auto equation=session->GetFunction(name,var);
                    out<<"{\"name\":\""<<name<<"/"<<var<<"\",\"value\":"
                        <<equation->Evaluate(x,dimension==2?r:.37,dimension==2?0:r,.3)<<'}';
                }
            out<<"]}";
        }
        out<<"]}\n";if(!out)throw std::runtime_error("cannot write expression samples");
        session->Finalise();return 0;
    }
    catch(const std::exception &error){std::cerr<<error.what()<<'\n';return 1;}
}
