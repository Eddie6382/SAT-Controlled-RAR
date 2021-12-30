#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <cassert>
#include <deque>
#include <map>
#include <string>
#include <unordered_map>
#include "cirMgr.h"
#include "cirGate.h"
#include "cirMA.h"   // define all the classes used in function SATRar
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
|    unsigned g1, g2. If (g1, g2) are both specified, than MA(w_t) is calculated
|                     If only specified g1, calculated MA(g_d)
|    bool init, whether find new set of dominators
|________________________________________________________________________________________________@*/
int
CirMA::computeSATMA(unsigned g1, unsigned g2=0, bool init=true)
{
   if (init) {
      _vecMA.clear();
      _vecMA.push_back(make_pair(g1, true));
      unordered_map<unsigned, unsigned> _nodeAppear;
   // Compute and assign Dominator, use partial _MA as assumptions
      if (g2 != 0)
         findDominators(g1, g2, _nodeAppear);
      else
         findDominators(0, g1, _nodeAppear);

      for (auto node: _dominators) {
         CirGate *g = _mgr->getGate(node);
         if (g->getType() == AIG_GATE) {
            unsigned In0Gid = g->getIn0Gate()->getGid();
            unsigned In1Gid = g->getIn1Gate()->getGid();
            bool isInv0 = g->getIn0().isInv();
            bool isInv1 = g->getIn1().isInv();

            if (_nodeAppear.find(In0Gid) == _nodeAppear.end())
               _vecMA.push_back(make_pair(In0Gid, (isInv0 != 1)));
            if (_nodeAppear.find(In1Gid) == _nodeAppear.end())
               _vecMA.push_back(make_pair(In1Gid, (isInv1 != 1)));
         }
      }
      
      assumeRelease();
      for (int i=_vecMA.size()-1; i>=0; --i)
         assumeProperty(_gid2Var[_vecMA[i].first], _vecMA[i].second);
      _solver->oneStepMA(_assump, init);
   }
   else {
      _solver->oneStepMA(_assump, init);
   }

   return 0;
}

int
CirMA::computeMA(unsigned g1, unsigned g2=0) 
{
   _MA.clear();
   _MA[g1] = true;
   unordered_map<unsigned, unsigned> _nodeAppear;
   
// Compute and assign Dominator
   if (g2 != 0)
      findDominators(g1, g2, _nodeAppear);
   else
      findDominators(0, g1, _nodeAppear);


   for (auto node: _dominators) {
      CirGate *g = _mgr->getGate(node);
      if (g->getType() == AIG_GATE) {
         unsigned In0Gid = g->getIn0Gate()->getGid();
         unsigned In1Gid = g->getIn1Gate()->getGid();
         bool isInv0 = g->getIn0().isInv();
         bool isInv1 = g->getIn1().isInv();

         if (_nodeAppear.find(In0Gid) == _nodeAppear.end())
            if (redundancyCheck(In0Gid, (isInv0 != 1))) return (int)In0Gid ;
         if (_nodeAppear.find(In1Gid) == _nodeAppear.end())
            if (redundancyCheck(In1Gid, (isInv1 != 1))) return (int)In1Gid;
      }
   }
// Iterately calculate sensitivity and propagativity
   size_t preMAsize = _MA.size();
   
   do {
      // sentivity
      deque<unsigned> Q;
      for (auto it: _MA) { 
         // only when an aig having "1's" can its input be forced assigned
         // Q is a subset of MA
         if (it.second)
            Q.push_back(it.first); 
      }
      vector<bool> IsVisit (_mgr->getNumTots(), false);

      while (!Q.empty()) {
         CirGate *g = _mgr->getGate(Q.front()); // g must has value 1
         Q.pop_front();
         if ((g->getType() == PI_GATE) || (g->getType() == CONST_GATE)) continue;
         vector<unsigned> ins = {g->getIn0Gate()->getGid(), g->getIn1Gate()->getGid()};
         vector<bool> isInvs = {g->getIn0().isInv(), g->getIn1().isInv()};
         for (size_t i=0; i<ins.size(); ++i) {
            // update Q
            if (!IsVisit[ins[i]]) {
               IsVisit[ins[i]] = true;
               if (isInvs[i] != 1) Q.push_back(ins[i]);  
            }
            // update M, and check if there is conflict
            if (redundancyCheck(ins[i], (isInvs[i] != 1)))
               return (int)ins[i];
         }
      }

      // propagativity
      for (auto it: _MA) {
         // Every element sholud be added, ga
         // It depends on fanout gate type and phase
         Q.push_back(it.first);
      }
      IsVisit.clear();
      for (size_t i=0; i<_mgr->getNumTots(); ++i) IsVisit.push_back(false);

      while (!Q.empty()) {
         unsigned gid = Q.front();
         Q.pop_front();
         GateList& fanouts = _mgr->getFanouts(gid);
         size_t nFanouts = fanouts.size();
         for (size_t i=0; i<nFanouts; ++i) {
            CirGate *fanout = fanouts[i];
            if (fanout->getType() == PO_GATE) continue;
            CirGateV in0 = fanout->getIn0();
            CirGateV in1 = fanout->getIn1();
            assert (in0.gate() != 0);
            bool fanoutInv = false;
            bool anotherInv;
            CirGate* anotherGate;
            // decide whether v's In0 or in1 conntects to u
            if (in0.gate() == _mgr->getGate(gid)) {
               fanoutInv = in0.isInv();
               anotherInv = in1.isInv();
               anotherGate = in1.gate();
            }
            else {
               assert (in1.gate() == _mgr->getGate(gid));
               fanoutInv = in1.isInv();
               anotherInv = in0.isInv();
               anotherGate = in0.gate();
            }
            
            // two case for a AND gate to progate
            // 1. both input having assigned already;  2. only one input, and is assigned 0
            bool canProp = false;
            if (_MA.find(anotherGate->getGid()) != _MA.end())
               canProp = true;
            else
               canProp = (_MA[gid] == fanoutInv);

            // update Q
            if (!IsVisit[fanout->getGid()]) {
               IsVisit[fanout->getGid()] = true;
               if (canProp) Q.push_back(fanout->getGid());
            }
            // update M
            if (canProp) {
               bool value = (_MA.find(anotherGate->getGid()) != _MA.end()) ? \
                  (bool)((_MA[anotherGate->getGid()]!=anotherInv) & (_MA[gid]!=fanoutInv)) : 0;
               if (redundancyCheck(fanout->getGid(), value)) {
                  printMA();
                  return (int)fanout->getGid();
               }
            }
         }
      }
      preMAsize = _MA.size();
   } while (preMAsize < _MA.size());

   printMA();
   return -1;
}

void
CirMA::findDominators(unsigned g0, unsigned src, unordered_map<unsigned, unsigned>& _nodeAppear)
{
   _allpaths.clear();
   _nodeAppear.clear();
   _dominators.clear();
   if (g0 != 0) {
      _nodeAppear[g0] = 1;
      _dominators.push_back(g0);
   }

   unsigned dst;
   for (unsigned i = 0, n = _mgr->getNumPOs(); i < n; ++i) {
      dst = _mgr->getPo(i)->getGid();
      IdList path;
      path.push_back(src);
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
*/
void
CirMA::printSATMA(int tabsize)
{
   map<int, vector<size_t> > decisionMA;
   for (size_t i=0; i<_mgr->getNumTots(); ++i) {  // iterate gate gid
      if (_mgr->getGate(i) == 0) continue;
      if (_solver->value(_gid2Var[i]) == l_Undef) continue;
      int level = _solver->atLevel(_gid2Var[i]);
      if (level != -1) {
         if (decisionMA.find(level) == decisionMA.end())
            decisionMA[level] = {i};
         else
            decisionMA[level].push_back(i);
      }
   }
   string tab = string(tabsize, ' ');
   
   for (auto it: decisionMA) {
      cout << tab << "level" << setw(3) << it.first << ":";
      for (auto i: it.second) {
         char neg = (_solver->value(_gid2Var[i])==l_True?' ': \
                (_solver->value(_gid2Var[i])==l_False?'!':'X'));
         cout << " " << neg << i;
      }
      cout << "\n";
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