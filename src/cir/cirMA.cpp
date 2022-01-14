#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <cassert>
#include <deque>
#include <map>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include "cirMgr.h"
#include "cirGate.h"
#include "cirMA.h"   // define all the classes used in function SATRar
#include "cirDef.h"
#include "Solver.h"

using namespace std;

void
CirMA::constructCNF()
{
// intialize variables
   for (size_t i=0; i<_mgr->getNumTots(); ++i)
      if (_mgr->getGate(i) != 0)
         _gid2Var[i] = newVar();

// Add circuit to sat solver
   for (size_t i=0; i<_mgr->getNumTots(); ++i) {
      CirGate* g = _mgr->getGate(i);
      if (g==0) continue;

      switch (g->getType()) {
         case AIG_GATE: {
            Var in0 = _gid2Var[g->getIn0Gate()->getGid()];
            Var in1 = _gid2Var[g->getIn1Gate()->getGid()];
            Var out = _gid2Var[i];
            bool isInv0 = g->getIn0().isInv();
            bool isInv1 = g->getIn1().isInv();
            addAigCNF(out, in0, isInv0, in1, isInv1);
            break;
         }
         case PO_GATE: {
            Var in0 = _gid2Var[g->getIn0Gate()->getGid()];
            Var out = _gid2Var[i];
            bool isInv0 = g->getIn0().isInv();
            addBufCNF(out, in0, isInv0);
            break;
         }
         default:
            break;
      }
   }
}

/*_________________________________________________________________________________________________
|  Description:
|    compute MA using Minisat solver
|  
|  Input:
|    unsigned g1, g2. If (g1, g2) are both specified, than MA(w_t) is calculated, otherwise if only specified g1, calculated MA(g_d)
|    bool init, whether find new set of dominators
|    bool copy, whether to copy implication from pointer c_solver
|________________________________________________________________________________________________@*/
Var
CirMA::computeSATMA(unsigned g1, unsigned g2=0, bool init=true, bool copy=false)
{
   Var conflict_var = -1;

   if (init && !copy) {
      unordered_map<unsigned, unsigned> _nodeAppear;

   // Compute and assign SA0 gate and its dominators, use partial _MA as assumptions
      _initMA.clear();
      _initMA.push_back(make_pair(g1, true));
      if (g2 != 0)
         findDominators(g1, g2, _nodeAppear);
      else
         findDominators(0, g1, _nodeAppear);
      // printDominators();

      for (auto node: _dominators) {
         CirGate *g = _mgr->getGate(node);
         if (g->getType() == AIG_GATE) {
            unsigned In0Gid = g->getIn0Gate()->getGid();
            unsigned In1Gid = g->getIn1Gate()->getGid();
            bool isInv0 = g->getIn0().isInv();
            bool isInv1 = g->getIn1().isInv();

            if (_nodeAppear.find(In0Gid) == _nodeAppear.end())
               _initMA.push_back(make_pair(In0Gid, (isInv0 != 1)));
            if (_nodeAppear.find(In1Gid) == _nodeAppear.end())
               _initMA.push_back(make_pair(In1Gid, (isInv1 != 1)));
         }
      }
      
      assumeRelease();
      for (int i=_initMA.size()-1; i>=0; --i)
         assumeProperty(_gid2Var[_initMA[i].first], _initMA[i].second);
      conflict_var = _solver->oneStepMA(_assump, init);
   }
   else if (init && copy) {
      assert(c_cirMA != NULL);
      assumeRelease();
      for (int i=c_cirMA->_initMA.size()-1; i>=0; --i) {
         pair<unsigned, bool> t = c_cirMA->_initMA[i];
         assumeProperty(_gid2Var[t.first],t.second);
      }
      _assump.pop(); _assump.pop();
      assumeProperty(_gid2Var[g1], true);
      conflict_var = _solver->oneStepMA(_assump, init);

      _dominators.clear();
      for (auto it: c_cirMA->_dominators)
         _dominators.push_back(it);
      _dominators.pop_front();
   }
   else {
      assert(g2 == 0);
      _assump.pop(); _assump.pop();
      assumeProperty(_gid2Var[g1], true);
      conflict_var = _solver->oneStepMA(_assump, init);
   }

   computePhiSet(g1);

   return conflict_var;
}

/*_________________________________________________________________________________________________
|  Description:
|    compute MA using Minisat solver
|  
|  Input:
|    unsigned g0, src. If both are specified, MA(w_t) is calculated and src is included in dominators
|                     If g0 not specified, find dominators of g_d, src is not included
|    bool init, whether find new set of dominators
|________________________________________________________________________________________________@*/
void
CirMA::findDominators(unsigned g0, unsigned src, unordered_map<unsigned, unsigned>& _nodeAppear)
{
   _allpaths.clear();
   _nodeAppear.clear();
   _dominators.clear();
   if (g0 != 0) { _nodeAppear[g0] = 1; }
   else { _nodeAppear[src] = 1; }

   unsigned dst;
   for (unsigned i = 0, n = _mgr->getNumPOs(); i < n; ++i) {
      dst = _mgr->getPo(i)->getGid();
      IdList path;
      if (g0 != 0) { path.push_back(src); }
      DFS(src, dst, path);
   }

   for (const auto& path : _allpaths) {
      for (const auto& node : path) {
         if (_nodeAppear.find(node) != _nodeAppear.end())
            _nodeAppear[node] += 1;
         else
            _nodeAppear[node] = 1;
      }
   }
   if (!_allpaths.size()) return;
   for (const auto& node: _allpaths[0])
      if (_nodeAppear[node] == _allpaths.size())
         _dominators.push_back(node);

}

void 
CirMA::DFS(unsigned src, unsigned dst, IdList& path)
{
   if (src == dst)
      _allpaths.push_back(path);
   else {
      GateList& fanouts = _mgr->getFanouts(src);
      for (size_t i=0; i<fanouts.size(); ++i) {
         path.push_back(fanouts[i]->getGid());
         DFS(fanouts[i]->getGid(), dst, path);
         path.pop_back();
      }
   }
}

void
CirMA::printPath()
{
   for (const auto& path : _allpaths) {
      cout << "Path : ";
      for (const auto& node : path)
            cout << node << " ";
      cout << endl;
   }
}

/*
print MA in decision order
   (unsigned) gid:     Target gate ID to compute MA
   (unsigned) tabsize: Indent
*/
void
CirMA::printSATMA(unsigned gid, unsigned tabsize)
{  
   string tab = string(tabsize, ' ');
   for (auto it: _decisionMA) {
      cout << tab << "level" << setw(3) << it.first << ": ";
      for (auto i: it.second) {
         if (!isInPhi(i)) continue;
         char neg = (_solver->value(_gid2Var[i])==l_True?' ': \
                (_solver->value(_gid2Var[i])==l_False?'!':'X'));
         cout << " " << neg << i;
      }
      cout << "\n";
   }
}

void 
CirMA::computePhiSet(unsigned gid)
{
   _decisionMA.clear();
   _Phi.clear();
   for (size_t i=0; i<_mgr->getNumTots(); ++i) {  // iterate gate gid
      if (_mgr->getGate(i) == 0) continue;
      if (_solver->value(_gid2Var[i]) == l_Undef) continue;
      int level = _solver->atLevel(_gid2Var[i]);
      if (level != -1) {
         if (_decisionMA.find(level) == _decisionMA.end())
            _decisionMA[level] = {i};
         else
            _decisionMA[level].push_back(i);
      }
   }

   unordered_set<unsigned> t;
   t.insert(gid);
   for (auto it: _dominators)
      t.insert(it);

   for (auto it: _decisionMA) {
      for (auto i: it.second) {
         if (t.find(i) != t.end()) continue;
         _Phi.insert(i);
      }
   }
}

void
CirMA::printDominators()
{
   cout << "Dominators:";
   for (size_t i=0; i<_dominators.size(); ++i)
      cout << " " << _dominators[i];
   cout << "\n";
}

void
CirMA::printMA()
{
   vector<bool> temp;
   cout << "MA: ";
   for (auto it: _MA) {
      cout << setw(3) << it.first;
      temp.push_back(it.second);
   }
   cout << "\n    ";
   for (auto it: temp)
      cout << setw(3) << it;
   cout <<"\n";
}

bool
CirMA::redundancyCheck(unsigned gid, bool value)
{
   if (_MA.find(gid) == _MA.end()) {
      _MA[gid] = value;
      return false;
   }
   else {
      if (_MA[gid] != value)
         return true;
      else
         return false;
   }
}