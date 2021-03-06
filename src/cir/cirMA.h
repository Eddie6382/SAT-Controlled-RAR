#ifndef CIR_MA_H
#define CIR_MA_H

#include <vector>
#include <unordered_map>
#include <map>
#include <unordered_set>
#include <iostream>
#include "cirGate.h"
#include "cirMgr.h"
#include "Solver.h"
#include <bits/stdc++.h>

using namespace std;

/*
   A class responsible for Mandatory Assignment
*/
class CirMA
{
   friend class CirMgr;
public:
   CirMA(CirMgr* _cirmgr): _mgr(_cirmgr), _solver(0), c_cirMA(0) {
      initialize();
      constructCNF();
   };
   ~CirMA() { if (_solver) { delete _solver; _solver=0; } }

   void setCounterpartSolver(CirMA* cir_wt, unsigned w_t) {
      c_cirMA = cir_wt;
      _solver->setCounterpartSolver(cir_wt->_solver, w_t, cir_wt->_dominators, cir_wt->_gid2Var);
   }
   int computeSATMA(unsigned, unsigned, bool, bool, int);
   int decisionGs(unsigned, bool);
   void computeSet(unsigned, int);
   void backtrace(int step) { _solver->cancelUntil(_solver->decisionLevel() - step); }
   void resetSolver() { _solver->cancelUntil(0); }

// print (debug)
   void printSATMA(unsigned, unsigned);
   void printMA();
   void printDominators();
// suport
   bool isInPhi(unsigned gid) const { return (_Phi.find(gid) != _Phi.end()); }
   char phase(unsigned i) { return (_solver->value(_gid2Var[i])==l_True?' ': \
                (_solver->value(_gid2Var[i])==l_False?'!':'X')); }
   bool isPos(unsigned i) { return _solver->value(_gid2Var[i])==l_True; }
   bool isNeg(unsigned i) { return _solver->value(_gid2Var[i])==l_False; }
   int nDecision() { return _solver->decisionLevel(); }
   bool isDet(unsigned i) { return _solver->value(_gid2Var[i])!=l_Undef; }
   unsigned var2Gid(Var var) {
      unsigned gid = 0;
      for (auto it: _gid2Var)
         if (it.second == var)
            gid = it.first;
      return gid;
   }

// Iterate
   inline size_t nDominators() { return _dominators.size(); }

   vector<unsigned> _dominators;
private:
// SAT related
//
   void reset() {
      if (_solver) { delete _solver; };
      _solver = new Solver();
      _gid2Var.clear();
      _assump.clear();
      _curVar = 0;
      _MA.clear();
      _initMA.clear();
      _dominators.clear();
   }
   void initialize() {
      reset();
      if (_curVar == 0) { _solver->newVar(); ++_curVar; }
   }
   // Constructing proof model
   // Return the Var ID of the new Var
   inline Var newVar() { _solver->newVar(); return _curVar++; }
   // fa/fb = true if it is inverted
   void addBufCNF(Var vf, Var va, bool fa) {
      vec<Lit> lits;
      Lit lf = Lit(vf);
      Lit la = fa? ~Lit(va): Lit(va);
      lits.push(la); lits.push(~lf);
      _solver->addClause(lits); lits.clear();
      lits.push(~la); lits.push(lf);
      _solver->addClause(lits); lits.clear();
   }
   // fa/fb = true if it is inverted
   void addAigCNF(Var vf, Var va, bool fa, Var vb, bool fb) {
      vec<Lit> lits;
      Lit lf = Lit(vf);
      Lit la = fa? ~Lit(va): Lit(va);
      Lit lb = fb? ~Lit(vb): Lit(vb);
      lits.push(la); lits.push(~lf);
      _solver->addClause(lits); lits.clear();
      lits.push(lb); lits.push(~lf);
      _solver->addClause(lits); lits.clear();
      lits.push(~la); lits.push(~lb); lits.push(lf);
      _solver->addClause(lits); lits.clear();
   }
   // For incremental proof, use "assumeSolve()"
   void assumeRelease() { _assump.clear(); }
   void assumeProperty(Var prop, bool val) {
      _assump.push(val? Lit(prop): ~Lit(prop));
   }
   bool assumpSolve() { return _solver->solve(_assump); }

   // For one time proof, use "solve"
   void assertProperty(Var prop, bool val) {
      _solver->addUnit(val? Lit(prop): ~Lit(prop));
   }
   bool solve() { _solver->solve(); return _solver->okay(); }

   // Functions about Reporting
   // Return 1/0/-1; -1 means unknown value
   int getValue(Var v) const {
      return (_solver->modelValue(v)==l_True?1:
               (_solver->modelValue(v)==l_False?0:-1)); }
   void printStats() const { const_cast<Solver*>(_solver)->printStats(); }

   void constructCNF();

   CirMgr* _mgr;
   unordered_map<unsigned, bool> _MA;     // These params are initial assignments of MA
   vector<pair<unsigned, bool>> _initMA;   // k (w_t's direct fanout) --> 0 (POs)

   // for sat solver
   Solver            *_solver;   // Pointer to a Minisat solver
   const CirMA* c_cirMA;
   Var               _curVar;    // Variable currently
   vec<Lit>          _assump;    // Assumption List for assumption solve

   // for final MA set
   map<int, vector<size_t> > _decisionMA;
   unordered_set<unsigned> _Phi; 

   unordered_map<unsigned, Var> _gid2Var;

};

#endif // CIR_MA_H