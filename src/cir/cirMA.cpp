#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <cassert>
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
|    k, from g_dk to g_d0
|________________________________________________________________________________________________@*/
Var
CirMA::computeSATMA(unsigned g1, unsigned g2=0, bool init=true, bool copy=false, int k=0)
{
   Var conflict_var = -1;

   if (init && !copy) {
      unordered_map<unsigned, unsigned> _nodeAppear;

   // Compute and assign SA0 gate and its dominators, use partial _MA as assumptions, dominators does not include source
      _initMA.clear();
      _initMA.push_back(make_pair(g1, true));
      // /*
      _dominators.clear();
      unordered_set<unsigned> isDom;
      if (g2 != 0) {
         for (auto it: _mgr->getDominators(g2))
            isDom.insert(it);
      } else {
         for (auto it: _mgr->getDominators(g1))
            if (it != g1) isDom.insert(it);
      }
      unsigned it_g = (g2!=0) ? g2 : g1;
      while(1) {
         GateList& fanouts = _mgr->getFanouts(it_g);
         size_t nFanouts = fanouts.size();
         if (isDom.find(it_g) != isDom.end())
            _dominators.push_back(it_g);
         if (nFanouts == 0) break;
         it_g = fanouts[0]->getGid();
      }
      for (auto node: _dominators) {
         CirGate *g = _mgr->getGate(node);
         if (g->getType() == AIG_GATE) {
            unsigned In0Gid = g->getIn0Gate()->getGid();
            unsigned In1Gid = g->getIn1Gate()->getGid();
            bool isInv0 = g->getIn0().isInv();
            bool isInv1 = g->getIn1().isInv();

            if (!_mgr->isTransitive(g1, In0Gid))
               _initMA.push_back(make_pair(In0Gid, (isInv0 != 1)));
            if (!_mgr->isTransitive(g1, In1Gid))
               _initMA.push_back(make_pair(In1Gid, (isInv1 != 1)));
         }
      }
      // */
      // printDominators();
      
      assumeRelease();
      for (int i=_initMA.size()-1; i>=0; --i)
         assumeProperty(_gid2Var[_initMA[i].first], _initMA[i].second);
      conflict_var = _solver->oneStepMA(_assump, init);

      computeSet(g1, k);
   }
   else if (init && copy) {
      assert(c_cirMA != NULL);
      assumeRelease();
      if (c_cirMA->_initMA.size() > 2) {
         for (int i=c_cirMA->_initMA.size()-1; i>=0; --i) {
            pair<unsigned, bool> t = c_cirMA->_initMA[i];
            assumeProperty(_gid2Var[t.first],t.second);
         }
         _assump.pop(); _assump.pop();
      }
      assumeProperty(_gid2Var[g1], true);
      conflict_var = _solver->oneStepMA(_assump, init);

      _dominators.clear();
      if (c_cirMA->_dominators.size() > 1) {
         for (size_t i=1; i<c_cirMA->_dominators.size(); ++i)
            _dominators.push_back(c_cirMA->_dominators[i]);
      }

      computeSet(g1, k);
   }
   else {
      _assump.pop(); _assump.pop();
      assumeProperty(_gid2Var[g1], true);
      assert(_assump.size());
      conflict_var = _solver->oneStepMA(_assump, init);   // It will cancel decision (backtrack)

      computeSet(g1, k);
   }

   

   return conflict_var;
}

/*_________________________________________________________________________________________________
|  Description:
|    make decision
|    <NOTE> no bracktrace on deciosn level; _assump is self defined var in class MA, no need to be modified
|  
|  Input:
|    unsigned gid
|    bool phase
|________________________________________________________________________________________________@*/
int 
CirMA::decisionGs(unsigned gs, bool phase) {
   Lit lit_gs = (phase) ? Lit(_gid2Var[gs]) : ~Lit(_gid2Var[gs]);
   return _solver->makeDecision(lit_gs);
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
CirMA::computeSet(unsigned gid, int k=0)
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
   for (; k<_dominators.size(); ++k)
      t.insert(_dominators[k]);

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