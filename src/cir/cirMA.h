#ifndef CIR_MA_H
#define CIR_MA_H

#include <vector>
#include <unordered_map>
#include <iostream>
#include <deque>
#include "cirGate.h"
#include "cirMgr.h"
#include "Solver.h"

using namespace std;

/*
   A class responsible for Mandatory Assignment
*/
class CirMA
{
public:
   CirMA(CirMgr* _cirmgr): _mgr(_cirmgr), _solver(0), c_cirMA(0) {
      initialize();
      constructCNF();
   };
   ~CirMA() { if (_solver) delete _solver; }

   void setCounterpartSolver(CirMA* cir_wt, unsigned w_t) {
      c_cirMA = cir_wt;
      _solver->setCounterpartSolver(cir_wt->_solver, w_t, cir_wt->_dominators, cir_wt->_gid2Var);
   }
   int computeSATMA(unsigned, unsigned, bool, bool);

// print (debug)
   void printPath();
   void printSATMA(unsigned, unsigned);
   void printMA();
   void printDominators();

// Iterate
   inline size_t nDominators() { return _dominators.size(); }

   deque<unsigned> _dominators;
private:
// MA
   void findDominators(unsigned, unsigned, unordered_map<unsigned, unsigned>&);
   void DFS(unsigned, unsigned, IdList&);
   bool redundancyCheck(unsigned, bool);
// SAT related
//
   void reset() {
      if (_solver) delete _solver;
      _solver = new Solver();
      _gid2Var.clear();
      _assump.clear();
      _curVar = 0;
      _allpaths.clear();
      _MA.clear();
      _vecMA.clear();
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

   const CirMgr* _mgr;
   unordered_map<unsigned, bool> _MA;
   vector<pair<unsigned, bool>> _vecMA;   // k (w_t's direct fanout) --> 0 (POs)

   // for sat solver
   Solver            *_solver;   // Pointer to a Minisat solver
   const CirMA* c_cirMA;
   Var               _curVar;    // Variable currently
   vec<Lit>          _assump;    // Assumption List for assumption solve

   // for dominators
   vector<IdList> _allpaths;

   unordered_map<unsigned, Var> _gid2Var;

};

#endif // CIR_MA_H