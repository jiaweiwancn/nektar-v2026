// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/RANSWallTiming.h>
#include <stdexcept>

int main()
{
    using namespace Nektar::RANSWallTiming;
    Scope process("process");
    {
        Context context(17,.125);Scope step("native.step");
        {Scope child("normal");}
        {Scope child("manual_stop");child.Stop();child.Stop();}
        for(int i=0;i<5;++i){Scope aggregate("buffered",true,true);Scope nested("buffered_child",true,true);}
        try{Scope child("exception");throw std::runtime_error("test");}
        catch(const std::runtime_error&){}
    }
    return 0;
}
