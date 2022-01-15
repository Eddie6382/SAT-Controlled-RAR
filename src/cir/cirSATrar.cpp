#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <utility>
#include <cassert>
#include <queue>
#include <sstream>
#include <fstream>
#include <string>
#include <regex>
#include "cirMgr.h"
#include "cirGate.h"
#include "cirMA.h" // define all the classes used in function SATRar
#include "cirDef.h"
#include "Solver.h"

using namespace std;

/*_________________________________________________________________________________________________
________________________________________________________________________________________________@*/
void CirMgr::SATRar()
{
   // For Minisat based approach
   _twoWayRepairList.clear();
   _rarWireRePairList.clear();
   _rarGateRepairList.clear();
   _w_tList.clear();
   
   
   CirMA MAw_t(this);
   CirMA MAg_d(this);

   int Ntar = 0;

   genDfsList();
   for (auto& in: _dfsList) {
      if (in->getType() != AIG_GATE) continue;
      GateList& fanouts = getFanouts(in->getGid());
      for (auto& out: fanouts) {
         if (out->getType() != AIG_GATE) continue;
         pair<unsigned, unsigned> w_t = make_pair(in->getGid(), out->getGid());
         _w_tList.push_back(w_t);
         if (SATRarOnWt(w_t, _w_tList.size()-1, MAw_t, MAg_d))
            Ntar++;
      }
   }

   cout << "\n----  Statistic  ----\n";
   cout << "2-Way rar alternatives: " << _twoWayRepairList.size() << "\n";
   cout << "RAR wire alternatives: " << _rarWireRePairList.size() << "\n";
   cout << "RAR gate alternatives: " << _rarGateRepairList.size() << "\n";
   cout << "#tar: " << Ntar << "\n";
   cout << "#alt: " << _twoWayRepairList.size() + _rarWireRePairList.size() + _rarGateRepairList.size() << endl;
}

void CirMgr::SATRarWrite(string& dir)
{
   cout << "writing files...";
   size_t n = 0;
   for (auto& repair: _rarGateRepairList) {
      string text = "";
      RARGateRepair(repair, text);
      string file_name = dir + "/repair" + to_string(n) + ".aag";
      ofstream outfile(file_name);
      outfile << text;
      ++n;
   }
   cout << "done\n";
}

/*_________________________________________________________________________________________________
________________________________________________________________________________________________@*/
bool CirMgr::SATRarOnWt(pair<unsigned, unsigned> w_t, int w_tIdx, CirMA& MAw_t, CirMA& MAg_d)
{
   bool return_val = false;

   cout << "\n  w_t(" << w_t.first << ", " << w_t.second << ")\n";
   MAw_t.computeSATMA(w_t.first, w_t.second, 1, 0);
   genFanInCone(getGate(w_t.first));
   // MAw_t.printSATMA(w_t.first, 2);
   MAg_d.setCounterpartSolver(&MAw_t, w_t.first);

   int i = 0;
   unsigned g_d = w_t.second;
   do {
   // only backtrack to level 𝑖 to find MAs. If no conflict, compute 𝑀𝐴(𝑔_𝑑𝑖)
      cout << "    g_d = " << g_d << "\n";
      Var conflict_var = MAg_d.computeSATMA(g_d, 0, (i == 0), (i == 0));
      // MAg_d.printSATMA(g_d, 4);

   // 2-way RAR, make it a wire between conflict_gate and g_d
      if (conflict_var != -1) {
         unsigned conflict_gid = MAw_t.var2Gid(conflict_var);
         bool phase = MAw_t.isNeg(conflict_gid);
         cout << "    conflict gid " << conflict_var << "\n";
         _twoWayRepairList.push_back(make_pair(hashWtIdxGd(w_tIdx, g_d), conflict_gid*2 + phase));
         return_val = true;
         break;
      }

   // Decision making, select g_s belonging Phi_wt - Phi_i
      vector<unsigned> decisions;
      unsigned nDecision = 0;

      for (const auto& gate: _dfsList) {
         unsigned gid = gate->getGid();
         if ( (gate->getType() == AIG_GATE) && \
                !gate->isInFanin() && \
                MAw_t.isInPhi(gid) && \
                !MAg_d.isInPhi(gid)) {

            assert(MAw_t.isDet(gid));
            
            ++nDecision;
            decisions.push_back(gid*2 + MAw_t.isNeg(gid));
            cout << "      decide " << MAw_t.phase(gid) << gid << "\n";
            conflict_var = MAg_d.decisionGs(gid, MAw_t.isPos(gid));

            if (conflict_var != -1) {
               unsigned conflict_gid = MAw_t.var2Gid(conflict_var);
               cout << "        conflict gid " << conflict_gid << "\n";
               if (find(decisions.begin(), decisions.end(), conflict_gid) != decisions.end()) { // RAR wire 
                  _rarWireRePairList.push_back(make_pair(hashWtIdxGd(w_tIdx, g_d), decisions));
               } else {
                  decisions.push_back(conflict_gid*2 + MAw_t.isNeg(conflict_gid));
                  _rarGateRepairList.push_back(make_pair(hashWtIdxGd(w_tIdx, g_d), decisions));
               }
               return_val = true;
               break;
            }
         }
      }
      MAg_d.backtrace(nDecision);

   // update g_d for next iteration
      g_d = MAg_d._dominators.front();
      MAg_d._dominators.pop_front();
      ++i;
   } while (!MAg_d._dominators.empty());

   MAg_d.resetSolver();
   MAw_t.resetSolver();

   return return_val;
}

void 
CirMgr::genFanInCone(CirGate *g)
{
   CirGate::setGlobalRef();
   queue<CirGate* > Q;
   Q.push(g);
   g->setToGlobalRef();
   while (!Q.empty()) {
      CirGate* g = Q.front();
      Q.pop();
      if (g->getType() != AIG_GATE) continue;

      CirGate* in0 = g->getIn0Gate();
      CirGate* in1 = g->getIn1Gate();
      CirGate* fanins[] = {in0, in1};
      for (size_t i=0; i<2; ++i) {
         CirGate *fanin = fanins[i];
         if (!(fanin->isGlobalRef())) {
            Q.push(fanin);
            fanin->setToGlobalRef();
         }
      }
   }
}

void 
CirMgr::RARGateRepair(pair<BigNum, vector<unsigned>>& repair, string& text)
{
   unsigned w_tFirst = _w_tList[hashToWtIdx(repair.first)].first;
   unsigned w_tSecond = _w_tList[hashToWtIdx(repair.first)].second;
   unsigned g_d = hashToGd(repair.first);
   vector<unsigned>& decisions = repair.second;

   size_t nAig = 0;

   stringstream outfile;
   for (size_t i = 0, n = _dfsList.size(); i < n; ++i)
      if (_dfsList[i]->isAig()) ++nAig;
   outfile << "aag " << _numDecl[VARS] << " " << _numDecl[PI] << " "
           << _numDecl[LATCH] << " " << _numDecl[PO] << " "
           << nAig << endl;
   for (size_t i = 0, n = _numDecl[PI]; i < n; ++i)
      outfile << (getPi(i)->getGid()*2) << endl;
   for (size_t i = 0, n = _numDecl[PO]; i < n; ++i)
      outfile << getPo(i)->getIn0().litId() << endl;
   for (size_t i = 0, n = _dfsList.size(); i < n; ++i) {
      CirGate *g = _dfsList[i];
      if (!g->isAig()) continue;
      outfile << g->getGid()*2 << " " << g->getIn0().litId() << " "
           << g->getIn1().litId() << endl;
   }
   for (size_t i = 0, n = _numDecl[PI]; i < n; ++i)
      if (getPi(i)->getName())
         outfile << "i" << i << " " << getPi(i)->getName() << endl;
   for (size_t i = 0, n = _numDecl[PO]; i < n; ++i)
      if (getPo(i)->getName())
         outfile << "o" << i << " " << getPo(i)->getName() << endl;
   outfile << "c" << endl;
   outfile << "aag repair, (w_t, g_d) = (" << w_tFirst*2 << ", " << g_d*2 << ")" << endl;

// replace header
   text = outfile.str();
   string newVARS = to_string(_numDecl[VARS] + decisions.size());
   string newAIG = to_string(nAig + decisions.size());
   regex e0 ("(aag)(\\s\\d+)(\\s\\d+)(\\s\\d+)(\\s\\d+)(\\s\\d+)");
   text = regex_replace(text, e0, "$1 "+newVARS+"$3$4$5 "+newAIG);

// replace w_t.first with const 1
   bool inv;
   if (getGate(w_tSecond)->getIn0().gate() == getGate(w_tFirst))
      inv = getGate(w_tSecond)->getIn0().isInv();
   else {
      assert(getGate(w_tSecond)->getIn1().gate() == getGate(w_tFirst));
      inv = getGate(w_tSecond)->getIn1().isInv();
   }

   regex e1 ("(\n"+to_string(w_tSecond*2)+")(\\s\\d+)\\s("+to_string(w_tFirst*2+(int)inv)+")");
   regex e2 ("(\n"+to_string(w_tSecond*2)+")\\s("+to_string(w_tFirst*2+(int)inv)+")");
   string const1 = to_string(1);
   text = regex_replace(text, e1, "$1$2 "+const1);
   text = regex_replace(text, e2, "$1 "+const1);

// modify g_d structure, and construct g_n
   unsigned gids[decisions.size()];
   for (int i=0; i<(int)decisions.size(); ++i) {
      gids[i] = _numDecl[VARS] + 1 + i;
   }
   regex e3 ("("+to_string(g_d*2)+")(\\s\\d+\\s\\d+\n)");
   string n_lit_gd= to_string(gids[0]*2);
   text = regex_replace(text, e3, n_lit_gd+"$2");

   string addcircuit = "\n";
   string s_in0;
   string s_in1;
   for (int i=1; i<(int)decisions.size(); ++i) {
      if (i==1) {
         s_in0 = to_string(decisions[0]);
         s_in1 = to_string(decisions[1]);
      }
      else {
         s_in0 = to_string(gids[i-1]);
         s_in1 = to_string(decisions[i]);
      }
      addcircuit += to_string(gids[i]*2) + " " + s_in0 + " " + s_in1 + "\n";
   }
   addcircuit += to_string(g_d*2) + " " + n_lit_gd + " " + to_string(gids[decisions.size()-1]*2 + 1);

   regex e4 ("("+n_lit_gd+")(\\s\\d+)(\\s\\d+)");
   text = regex_replace(text, e4, "$1$2$3"+addcircuit);

}

void 
CirMgr::findGlobalDominators()
{
   // here dominators includes gid
   _globalDominators.clear();
   
   for (int i=_dfsList.size()-1; i>=0; --i) {
      const CirGate* g = _dfsList[i]; 
      unsigned gid = g->getGid();

      if (g->getType() == PO_GATE)
         _globalDominators[gid] = vector<unsigned> {gid};
      else if ((g->getType() == PI_GATE) || (g->getType() == AIG_GATE)) {
         vector<unsigned> gidDom;
         GateList& fanouts = getFanouts(gid);
         size_t nFanouts = fanouts.size();
         assert(nFanouts);

         gidDom = _globalDominators[fanouts[0]->getGid()];
            
         for (size_t j=1; j<fanouts.size(); ++i) {
            vector<unsigned> foDom = _globalDominators[fanouts[j]->getGid()];

            vector<unsigned> vec_intersect;
            set_intersection(gidDom.begin(), gidDom.end(), foDom.begin(), foDom.end(), back_inserter(vec_intersect));
            gidDom = vec_intersect;
         }
         gidDom.push_back(gid);
         sort(gidDom.begin(), gidDom.end());  // important!!
         _globalDominators[gid] = gidDom;
      }
   }
}
