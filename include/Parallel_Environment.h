#define CPP_FILE

#ifndef __Parallel_Environment_h_wanghan__
#define __Parallel_Environment_h_wanghan__

#include "mpi.h"
#include "MDException.h"


namespace Parallel {
  class MDExcptDimsNotConsistentWithNProc : public MDException {};
  
  class Environment 
  {
    MPI_Comm commCart;
    int dims[3];
    int myRank;
    int numProc;
public:
    Environment (int * argc, char *** argv);
    ~Environment();
public:
    void init (const int division[3]);
  };


}

#endif

