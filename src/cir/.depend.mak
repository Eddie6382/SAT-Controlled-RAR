cirGate.o: cirGate.cpp cirGate.h cirDef.h ../../include/Solver.h \
 ../../include/SolverTypes.h ../../include/Global.h \
 ../../include/VarOrder.h ../../include/Heap.h ../../include/Proof.h \
 ../../include/File.h cirMgr.h ../../include/util.h ../../include/rnGen.h \
 ../../include/myUsage.h
cirMA.o: cirMA.cpp cirMgr.h cirDef.h cirGate.h ../../include/Solver.h \
 ../../include/SolverTypes.h ../../include/Global.h \
 ../../include/VarOrder.h ../../include/Heap.h ../../include/Proof.h \
 ../../include/File.h cirMA.h
cirMgr.o: cirMgr.cpp cirMgr.h cirDef.h cirGate.h ../../include/Solver.h \
 ../../include/SolverTypes.h ../../include/Global.h \
 ../../include/VarOrder.h ../../include/Heap.h ../../include/Proof.h \
 ../../include/File.h ../../include/util.h ../../include/rnGen.h \
 ../../include/myUsage.h
cirSATrar.o: cirSATrar.cpp cirMgr.h cirDef.h cirGate.h \
 ../../include/Solver.h ../../include/SolverTypes.h \
 ../../include/Global.h ../../include/VarOrder.h ../../include/Heap.h \
 ../../include/Proof.h ../../include/File.h cirMA.h
cirCmd.o: cirCmd.cpp cirMgr.h cirDef.h cirGate.h ../../include/Solver.h \
 ../../include/SolverTypes.h ../../include/Global.h \
 ../../include/VarOrder.h ../../include/Heap.h ../../include/Proof.h \
 ../../include/File.h cirCmd.h ../../include/cmdParser.h \
 ../../include/cmdCharDef.h ../../include/util.h ../../include/rnGen.h \
 ../../include/myUsage.h
