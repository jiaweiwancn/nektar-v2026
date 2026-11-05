// Exercise the actual R12 admission/scan schedule used by CoupledDirectionR4.
#include "RANSWallCutBudget.h"
#include <cassert>
#include <iostream>
using Nektar::RANSWallAlgebra::CutBudgetR12;
struct Result {bool accepted;int scans,admissions;};
Result solve(int limit,int required)
{
    CutBudgetR12 budget{limit};
    int admitted=0,scans=0;
    for(int pass=0;pass<budget.Scans();++pass)
    {
        ++scans;
        // Reconstruct from the current working set before invoking the oracle.
        if(admitted>=required)return {true,scans,admitted};
        if(!budget.CanAdmit(admitted))return {false,scans,admitted};
        ++admitted;
    }
    assert(false);return {};
}
int main()
{
    const auto completed=solve(3,3);
    assert(completed.accepted && completed.scans==4 && completed.admissions==3);
    const auto incomplete=solve(3,4);
    assert(!incomplete.accepted && incomplete.scans==4 && incomplete.admissions==3);
    const auto zero=solve(0,0);
    assert(zero.accepted && zero.scans==1 && zero.admissions==0);
    std::cout<<"R12 final-cut completion/rejection/no-extra-admission passed\n";
}
