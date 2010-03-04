#define CPP_FILE

#ifndef __Parallel_Environment_h_wanghan__
#define __Parallel_Environment_h_wanghan__

#include "mpi.h"
#include "MDException.h"


namespace Parallel {
  class MDExcptDimsNotConsistentWithNProc : public MDException {};

  enum coordIndex {
    CoordXIndex =	0,
    CoordYIndex	=	1,
    CoordZIndex	=	2
  };
  
  class Environment 
  {
    MPI_Comm commCart;
    int dims[3];
    int myRank_;
    int numProc_;
    bool inited;
public:
    Environment (int * argc, char *** argv);
    ~Environment();
public:
    void init (const int division[3]);
    int myRank  () const {return myRank_;}
    int numProc () const {return numProc_;}
    // const MPI_Comm * communicator () const {return commCart;}
    const MPI_Comm & communicator () const {return commCart;}
    void cartCoordToRank (const int & ix,
			  const int & iy,
			  const int & iz,
			  int & rank) const;
    void randToCartCoord (const int & rank,
			  int & ix,
			  int & iy,
			  int & iz) const;
    void numProcDim (int & nx,
		     int & ny,
		     int & nz) const;
    // void numProcDim (int * dims);
    
  };


}

#endif

