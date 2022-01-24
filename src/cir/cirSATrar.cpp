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
#include <sys/resource.h>

using namespace std;

long get_mem_usage() {
    struct rusage mysuage;
    getrusage(RUSAGE_SELF, &mysuage);
    return mysuage.ru_maxrss;
}

/*_________________________________________________________________________________________________
________________________________________________________________________________________________@*/
void CirMgr::SATRar(int verb)
{
   clock_t start = clock();
   _verbose = verb;
   
   // For Minisat based approach
   _RepairList.clear();
   _w_tList.clear();
   _repairCat[0] = 0;
   _repairCat[1] = 0;
   _repairCat[2] = 0;
   
   
   if (!_MAw_t) _MAw_t = new CirMA (this);
   if (!_MAg_d) _MAg_d = new CirMA (this);

   int Ntar = 0;
   int Nwire = 0;
    
   cout << "\n----  Preprocess ----\n";
   genDfsList();
   findGlobalDominators();
   findTransitiveClosure();
   cout << "Done\n----  Run SATrar ----\n";
   
   for (auto& in: _dfsList) {
      if (in->getType() != AIG_GATE) continue;
      GateList& fanouts = getFanouts(in->getGid());
      for (auto& out: fanouts) {
         if (out->getType() != AIG_GATE) continue;
         Nwire++;
         pair<unsigned, unsigned> w_t = make_pair(in->getGid(), out->getGid());
         _w_tList.push_back(w_t);
         if (SATRarOnWt(w_t, _w_tList.size()-1, *_MAw_t, *_MAg_d))
            Ntar++;
      }
   }

   cout << "\n----  Statistic  ----\n";
   cout << "Total time usage: " << double(clock()-start)/CLOCKS_PER_SEC << " s\n";
   cout << "#wire: "  << Nwire << "\n";
   cout << "2-Way rar alternatives: " << _repairCat[0] << "\n";
   cout << "RAR wire alternatives: " << _repairCat[1] << "\n";
   cout << "RAR gate alternatives: " << _repairCat[2] << "\n";
   cout << "#tar: " << Ntar << "\n";
   cout << "#alt: " << _RepairList.size() << endl;
}

void CirMgr::SATRarWrite(string& dir)
{
   cout << "writing files...";
   size_t n = 0;
   for (auto& repair: _RepairList) {
      string text = "";
      SATRARRepair(repair, text);
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

   cout << "  w_t(" << w_t.first << ", " << w_t.second << ")\n";
   MAw_t.computeSATMA(w_t.first, w_t.second, 1, 0, 0);
   genFanInCone(getGate(w_t.first));
   if (_verbose > 1) MAw_t.printSATMA(w_t.first, 2);
   MAg_d.setCounterpartSolver(&MAw_t, w_t.first);

   int i = 0;
   unsigned g_d = w_t.second;

   do {
   // only backtrack to level 𝑖 to find MAs. If no conflict, compute 𝑀𝐴(𝑔_𝑑𝑖)
      if (_verbose > 0) cout << "    g_d = " << g_d << "\n";
      Var conflict_var = MAg_d.computeSATMA(g_d, 0, (i == 0), (i == 0), i);
      if (_verbose > 1) MAg_d.printSATMA(g_d, 4);

   // 2-way RAR, make it a wire between conflict_gate and g_d
      
      vector<unsigned> decisions;

      if (conflict_var != -1) {
         unsigned conflict_gid = MAw_t.var2Gid(conflict_var);
         bool phase = MAw_t.isNeg(conflict_gid);
         decisions.push_back(conflict_gid*2 + phase);
         cout << "    conflict gid " << conflict_gid << "\n";
         _RepairList.push_back(make_pair(hashWtIdxGd(w_tIdx, g_d), decisions));
         return_val = true;
         _repairCat[0]++;
         break;
      }

   // Decision making, select g_s belonging Phi_wt - Phi_i
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
            if (_verbose > 1) cout << "      decide " << MAw_t.phase(gid) << gid << "\n";
            conflict_var = MAg_d.decisionGs(gid, MAw_t.isPos(gid));

            if (conflict_var != -1) {
               unsigned conflict_gid = MAw_t.var2Gid(conflict_var);
               cout << "        conflict gid " << conflict_gid << "\n";
               if (find(decisions.begin(), decisions.end(), conflict_gid) != decisions.end()) { // RAR wire 
                  _RepairList.push_back(make_pair(hashWtIdxGd(w_tIdx, g_d), decisions));
                  _repairCat[1]++;
               } else { // RAR gate
                  decisions.push_back(conflict_gid*2 + MAw_t.isNeg(conflict_gid));
                  _RepairList.push_back(make_pair(hashWtIdxGd(w_tIdx, g_d), decisions));
                  _repairCat[2]++;
               }
               return_val = true;
               break;
            }
         }
      }
      MAg_d.backtrace(nDecision);

   // update g_d for next iteration
      if (MAg_d._dominators.empty()) break;
      g_d = MAg_d._dominators[i];
      ++i;
   } while (MAg_d._assump.size() >= 2);  // No dominators, break;

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
CirMgr::SATRARRepair(pair<BigNum, vector<unsigned>>& repair, string& text)
{
   unsigned w_tFirst = _w_tList[hashToWtIdx(repair.first)].first;
   unsigned w_tSecond = _w_tList[hashToWtIdx(repair.first)].second;
   unsigned g_d = hashToGd(repair.first);
   vector<unsigned>& decisions = repair.second;

   size_t nAig = 0;

// main body
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
   outfile << "c\naag repair, redundant wire = (" << w_tFirst << ", " << w_tSecond << ")" << "\n";
   outfile << "            g_d = " << g_d << "\n";
   outfile << "            conflict decisions = ";
   for (int i=0; i<(int)decisions.size(); ++i)
      outfile << decisions[i]/2 << " ";
   outfile << endl;

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
   string const_in = (inv)? to_string(1): to_string(0);
   text = regex_replace(text, e1, "$1$2 "+const_in);
   text = regex_replace(text, e2, "$1 "+const_in);

// modify g_d structure, and construct g_n
   unsigned gids[decisions.size()];
   for (int i=0; i<(int)decisions.size(); ++i) {
      gids[i] = _numDecl[VARS] + 1 + i;
   }
   regex e3 ("("+to_string(g_d*2)+")(\\s\\d+\\s\\d+\n)");
   string n_lit_gd= to_string(gids[0]*2);
   text = regex_replace(text, e3, n_lit_gd+"$2");

   string addcircuit = "\n";
   if (decisions.size() >= 2) {
      string s_in0;
      string s_in1;
      for (int i=1; i<(int)decisions.size(); ++i) {
         if (i==1) {
            s_in0 = to_string(decisions[0]);
            s_in1 = to_string(decisions[1]);
         }
         else {
            s_in0 = to_string(gids[i-1]*2);
            s_in1 = to_string(decisions[i]);
         }
         addcircuit += to_string(gids[i]*2) + " " + s_in0 + " " + s_in1 + "\n";
      }
      addcircuit += to_string(g_d*2) + " " + n_lit_gd + " " + to_string(gids[decisions.size()-1]*2 + 1);
   }
   else {   // <= 1
      // TODO
      bool isNeg = decisions[0] % 2;
      unsigned g_v = (isNeg) ? decisions[0]-1 : decisions[0]+1;
      addcircuit += to_string(g_d*2) + " " + n_lit_gd + " " + to_string(g_v);
   }

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
            
         for (size_t j=1; j<fanouts.size(); ++j) {
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

void 
CirMgr::findTransitiveClosure() {
   if (_graph != 0) return;
   if (_dfsList.empty()) genDfsList();
      
   _graph = new Graph(getNumTots());
   for (auto u: _dfsList) {
      GateList& fanouts = getFanouts(u->getGid());
      for (auto v: fanouts) {
         _graph->addEdge(u->getGid(), v->getGid());
      }
   }
   _graph->transitiveClosure();
}