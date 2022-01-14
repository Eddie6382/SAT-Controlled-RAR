#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <utility>
#include <cassert>
#include <queue>
#include "cirMgr.h"
#include "cirGate.h"
#include "cirMA.h" // define all the classes used in function SATRar
#include "cirDef.h"
#include "Solver.h"

using namespace std;

/*

*/
void CirMgr::SATRar()
{
   // For Minisat based approach
   cout << "Begin SATRar Process\n";
   cout << "1. Given w_t 12\n";
   CirMA MAw_t(this);
   CirMA MAg_d(this);

   pair<unsigned, unsigned> w_t = make_pair(12, 13);

   cout << "  Decide w_t(" << w_t.first << ", " << w_t.second << ")\n  MA(w_t)\n";
   MAw_t.computeSATMA(w_t.first, w_t.second, 1, 0);
   MAw_t.printSATMA(w_t.first, 2);
   MAg_d.setCounterpartSolver(&MAw_t, w_t.first);
   cout << "\n  2. For each g_d belonging w_t's dominators (i=k down to 0)\n";

   int i = 0;
   unsigned g_d = w_t.second;
   do {
   // only backtrack to level 𝑖 to find MAs. If no conflict, generate 𝑀𝐴(𝑔_𝑑𝑖)
      cout << endl;
      cout << "    MA(g_d = " << g_d << ")\n";
      Var conflict_var = MAg_d.computeSATMA(g_d, 0, (i == 0), (i == 0));
      MAg_d.printSATMA(g_d, 4);
      if (conflict_var == -1) cout << "    No conflict\n";
      else cout << "    conflict_var " << conflict_var << "\n";

   // 2-way RAR, make it a wire between conflict_gate and g_d
      if (conflict_var != -1) {
         _numDecl[VARS] += 1;
         unsigned conflict_gid = MAw_t.var2Gid(conflict_var);
         assert(conflict_gid != 0);
         GateList& v_FOs = getFanouts(conflict_gid);
         CirGate* g_add = new CirAigGate(_numDecl[VARS], 0);
         v_FOs.push_back(g_add);

         static_cast<CirAigGate*>(g_add)->setIn0(conflict_gid);
         static_cast<CirAigGate*>(g_add)->setIn1(g_d);
         GateList& gd_FOs = getFanouts(g_d);
         _totGateList[_numDecl[VARS]] = g_add;
         GateList& g_add_FOs = getFanouts(_numDecl[VARS]);
         g_add_FOs.assign(gd_FOs.begin(), gd_FOs.end());
         gd_FOs.clear();
         gd_FOs.push_back(g_add);
      }

   // Decision making, select g_s belonging Phi_wt - Phi_i
      genDfsList();
      cout << "\n      3. Decide g_s\n";
      for (const auto& gate: _dfsList) {
         unsigned gid = gate->getGid();
         if (MAw_t.isInPhi(gid) && !MAg_d.isInPhi(gid)) {
            cout << "        Select " << MAw_t.phase(gid) << gid << "\n";
         }
      }

   // update g_d for next iteration
      g_d = MAg_d._dominators.front();
      MAg_d._dominators.pop_front();
      ++i;
   } while (!MAg_d._dominators.empty());
}

void 
CirMgr::genFanoutCone(CirGate *g)
{
   CirGate::setGlobalRef();
   queue<unsigned> Q;
   Q.push(g->getGid());
   while (!Q.empty()) {
      unsigned gid = Q.front();
      Q.pop();
      GateList& fanouts = getFanouts(gid);
      size_t nFanouts = fanouts.size();
      for (size_t i=0; i<nFanouts; ++i) {
         CirGate *fanout = fanouts[i];
         if (!(fanout->isGlobalRef())) {
            Q.push(fanout->getGid());
            fanout->setToGlobalRef();
         }
      }
   }
}
