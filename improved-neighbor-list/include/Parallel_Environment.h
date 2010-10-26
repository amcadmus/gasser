#ifndef __Parallel_Environment_h_wanghan__
#define __Parallel_Environment_h_wanghan__

#define MPI_CODE
#include "mpi.h"
#include "common.h"
// #include "MDException.h"

namespace Parallel {
  class MDExcptDimsNotConsistentWithNProc : public MDException {};
  class MDExcptNoActiveProcess : public MDException {};
  
  enum coordIndex {
    CoordXIndex =	0,
    CoordYIndex	=	1,
    CoordZIndex	=	2
  } ;
  
  class Environment 
  {
    static MPI_Comm commActive;
    static MPI_Comm commCart;
    // static int rank_worldRoot;
    static int dims[3];
    static int myRank_;
    static int numProc_;
    static int active;
    static int inited;
    static unsigned cellCapacity;
private:
    static void initCart (const int & nx,
			  const int & ny,
			  const int & nz);
public:
    static void init_mpi (int * argc, char *** argv);
    static void init_env (const unsigned & cellCapacity = 64,
			  const char * deviceName = "Device Emulation (CPU)",
			  const int & nx = 0,
			  const int & ny = 0,
			  const int & nz = 0);
    static void finalize ();
public:
    static int isActive () {return active;}    
    static int myRank   () {return myRank_;}
    static int numProc  () {return numProc_;}
    static const MPI_Comm & communicator () {return commCart;}
    static void cartCoordToRank (const int & ix,
				 const int & iy,
				 const int & iz,
				 int & rank) ;
    static void rankToCartCoord (const int & rank,
				 int & ix,
				 int & iy,
				 int & iz) ;
    static void numProcDim (int & nx,
			    int & ny,
			    int & nz) ;
    static void neighborProcIndex (int direction,
				   int displacement,
				   int & src,
				   int & dest);
    static unsigned getCellCapacity () {return cellCapacity;}
public:    
    static void barrier ();
    static void abort (int errorCode);
    
  };
}

 
#endif

