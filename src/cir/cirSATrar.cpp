#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <utility>
#include <cassert>
#include "cirMgr.h"
#include "cirGate.h"
#include "cirMA.h" // define all the classes used in function SATRar
#include "Solver.h"

using namespace std;

/*

*/
void CirMgr::SATRar()
{
   // For Minisat based approach
   cout << "Begin SATRar Process\n";
   cout << "1. Find MA for w_t and g_d\n\n";
   CirMA MAw_t(this);
   CirMA MAg_d(this);

   pair<unsigned, unsigned> w_t = make_pair(12, 13);

   cout << "  Decide w_t(" << w_t.first << ", " << w_t.second << ")\n  MA(w_t)\n";
   MAw_t.computeSATMA(w_t.first, w_t.second, true);
   MAw_t.printSATMA(w_t.first, 2);
   MAg_d.setCounterpartSolver(&MAw_t, w_t.first);
   cout << "\n  For each g_d belonging w_t's dominators\n";

   int i = 0;
   unsigned gid = w_t.second;
   do {
   // only backtrack to level 𝑖 to find MAs. If no conflict, generate 𝑀𝐴(𝑔_𝑑𝑖)
      cout << "    MA(g_d = " << gid << ")\n";
      MAg_d.computeSATMA(gid, 0, (i == 0));
      MAg_d.printSATMA(gid, 4);
      cout << endl;

      gid = MAg_d._dominators.front();
      MAg_d._dominators.pop_front();
      ++i;

   // 2-way RAR
   
   } while (!MAg_d._dominators.empty());
}
