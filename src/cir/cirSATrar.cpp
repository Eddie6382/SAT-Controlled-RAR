#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <utility>
#include <cassert>
#include "cirMgr.h"
#include "cirGate.h"
#include "cirMA.h"   // define all the classes used in function SATRar
#include "Solver.h"

using namespace std;

void
CirMgr::SATRar()
{
// For Minisat based approach
   cout << "Begin SATRar Process\n";
   cout << "1. Find MA for w_t and g_d\n\n";
   CirMA MAw_t(this);
   CirMA MAg_d(this);

   pair<unsigned, unsigned> w_t = make_pair(12, 13);

   cout << "  Decide w_t("<< w_t.first << ", " << w_t.second << ")\n  MA(w_t)\n";
   MAw_t.computeSATMA(w_t.first, w_t.second, true);
   MAw_t.printSATMA(2);
   MAg_d.setCounterpartSolver(&MAw_t);
   cout << "\n  For each g_d belonging w_t's dominators\n";

   // (_dominators[0], _dominators[1]) are two end of w_t
   for (size_t i=1; i<MAw_t.nDominators(); ++i) {
      unsigned gid = MAw_t.getDominators()[i];
      cout << "    MA(g_d=" << gid << ")\n";
      MAg_d.computeSATMA(gid, 0, (i==1));
      MAg_d.printDominators();
      MAg_d.printSATMA(4);
      cout << endl;
   }


// For handmade MA
   // pair<unsigned, unsigned> w_t (12, 13);
   // unsigned g_d = 7;

   // int conflict_g_d = MAw_t.computeMA(w_t.first, w_t.second);
   // MAg_d.computeMA(g_d, 0);
}

