/****************************************************************************
  FileName     [ cirMgr.h ]
  PackageName  [ cir ]
  Synopsis     [ Define circuit manager ]
  Author       [ Chung-Yang (Ric) Huang ]
  Copyright    [ Copyleft(c) 2008-present LaDs(III), GIEE, NTU, Taiwan ]
****************************************************************************/

#ifndef CIR_MGR_H
#define CIR_MGR_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>

using namespace std;

// TODO: Feel free to define your own classes, variables, or functions.

#include "cirDef.h"
#include "cirMA.h"

extern CirMgr *cirMgr;

class CirMgr
{
   enum CirMgrFlag { NO_FEC = 0x1 };
   enum ParsePorts { VARS = 0, PI, LATCH, PO, AIG, TOT_PARSE_PORTS };

public:
   CirMgr(): _flag(0), _reserved_storage(10), _piList(0), _poList(0), _totGateList(0),
             _fanoutInfo(0), _simLog(0), _verbose(0) { }
   ~CirMgr() { deleteCircuit(); }

   // Access functions
   CirGate* operator [](unsigned gid) const { return _totGateList[gid]; }
   // return '0' if "gid" corresponds to an undefined gate.
   CirGate* getGate(unsigned gid) const {
      if (gid >= getNumTots()) { return 0; } return _totGateList[gid]; }
   CirGate* litId2Gate(unsigned litId) const { return getGate(litId/2); }
   CirGateV litId2GateV(unsigned litId) const;
   bool isFlag(CirMgrFlag f) const { return _flag & f; }
   void setFlag(CirMgrFlag f) const { _flag |= f; }
   void unsetFlag(CirMgrFlag f) const { _flag &= ~f; }
   void resetFlag() const { _flag = 0; }
   unsigned getNumPIs() const { return _numDecl[PI]; }
   unsigned getNumPOs() const { return _numDecl[PO]; }
   unsigned getNumTots() const { return _numDecl[VARS] + _numDecl[PO] + 1; }
   CirPiGate* getPi(unsigned i) const { return _piList[i]; }
   CirPoGate* getPo(unsigned i) const { return _poList[i]; }
   GateList& getFanouts(unsigned i) const { return _fanoutInfo[i]; }

   // Member function about SAT-Rar
   void SATRar(int);
   bool SATRarOnWt(pair<unsigned, unsigned>, int, CirMA&, CirMA&);
   void SATRarWrite(string& );
   void SATRARRepair(pair<BigNum, vector<unsigned>>& , string& );
   void findGlobalDominators();
   inline BigNum hashWtIdxGd(unsigned w_tIdx, unsigned g_d) { return (BigNum)w_tIdx * (BigNum)getNumTots() + (BigNum)g_d; }
   inline unsigned hashToWtIdx(BigNum hashVal) { return (unsigned)(hashVal / (BigNum)getNumTots()); }
   inline unsigned hashToGd(BigNum hashVal) { return (unsigned)(hashVal % (BigNum)getNumTots()); }
   IdList& getDominators(unsigned i) { return _globalDominators.at(i); }


   // Member functions about circuit construction
   bool readCircuit(const string&);
   void deleteCircuit();
   void genConnections();
   void genDfsList();
   // Member functions about circuit fanout cone
   void genFanInCone(CirGate* g);
   void genFanOutCone(CirGate* g);
//   void updateUndefList();
   void checkFloatList();
   void checkUnusedList();
   size_t checkConnectedGate(size_t);

   // Member functions about circuit optimization
//   void sweep();
//   void optimize();
//   bool checkAigOptimize(CirGate*, const CirGateV&, const CirGateV&,
//                         CirGateV&) const;
//   void deleteAigGate(CirGate *);
//   void deleteUndefGate(CirGate *);


   // Member functions about simulation
//   void randomSim();
//   void fileSim(ifstream&);
   void setSimLog(ofstream *logFile) { _simLog = logFile; }

   // Member functions about fraig
//   void strash();
   IdList* getFECGrps(size_t i) const { return _fecGrps[i]; }
//   void printFEC() const;
//   void fraig();

   // Member functions about circuit reporting
   void printSummary() const;
   void printNetlist() const;
   void printPIs() const;
   void printPOs() const;
   void printFloatGates() const;
//   void printFECPairs() const;
   void writeAag(ostream&) const;
//   void writeGate(ostream&, CirGate*) const;

   // Member functions about flags

   // for hidden command
//   bool createMiter(CirMgr*, CirMgr*);

   static CirGate *_const0;
   GateList            _dfsList;

private:
   unsigned            _numDecl[TOT_PARSE_PORTS];
   mutable unsigned    _flag;
   unsigned            _reserved_storage;
   PiArray             _piList;
   PoArray             _poList;
   // IDs in _undefList are NOT sorted!!
   IdList              _undefList;
   // Make sure the IDs of the following lists are sorted!!
   IdList              _floatList;  // gates with fanin(s) undefined
   IdList              _unusedList; // gates defined but not used
   GateArray           _totGateList;
   GateList           *_fanoutInfo;
   vector<IdList*>     _fecGrps;  // store litId; FECHash<GatePValue, IdList*>
//   SimVector           _fecVector;
   ofstream           *_simLog;
   CirMA              *_MAw_t;
   CirMA              *_MAg_d;
   int                 _verbose;

   // SATRar repair
   vector<pair<BigNum, vector<unsigned>>> _RepairList;
   vector<pair<unsigned, unsigned>> _w_tList;
   int _repairCat[3];

   unordered_map<unsigned, vector<unsigned>>  _globalDominators;

   // private member functions for circuit parsing
   bool parseHeader(ifstream&);
   bool parseInput(ifstream&);
   bool parseLatch(ifstream&);
   bool parseOutput(ifstream&);
   bool parseAig(ifstream&);
   bool parseSymbol(ifstream&);
   bool parseComment(ifstream&);
   bool checkId(unsigned&, const string&);
   CirGate* checkGate(unsigned&, ParsePorts, const string&);

   // private member functions for circuit optimization
   CirGateV constSimplify(CirGate *, const CirGateV&, const CirGateV&) const;

   // private member functions about simulation
//   void setRandPPattern() const;
//   void setPPattern(unsigned *const patterns) const;
//   unsigned gatherPatterns(ifstream&, unsigned *const, size_t);
//   void pSim1Pattern() const;
//   void outputSimLog(size_t nPatterns = 32);

   // private member functions about FRAIG
//   bool simAndCheckFEC();
//   bool initFEC();
//   bool checkFEC();
//   bool checkFECRecur(const IdList&, vector<IdList*>&);
//   void finalizeFEC();
//   void simplifyFECGrps();
//   void clearFECGrps();
//   void initProofModel(SatSolver&, const GateList&);
//   bool satCheckConst(SatSolver&, CirGate *, bool, unsigned *const);
//   bool satCheckFEC(SatSolver&, CirGate *, CirGate *, bool, unsigned *const);
//   void getSatAssignment(SatSolver&, unsigned *const) const;
//   void simplifyByEQ();
//   void updateFECbySatPattern(unsigned *const);
};

#endif // CIR_MGR_H
