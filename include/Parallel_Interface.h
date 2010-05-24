#ifndef __Parallel_interface_h__
#define __Parallel_interface_h__

namespace Parallel{
  namespace Interface{
    void initMPI         (int * argc, char *** argv);
    void initEnvironment (const unsigned & cellCapacity = 64,
			  const char * deviceName = "Tesla C1060",
			  const int & nx = 0,
			  const int & ny = 0,
			  const int & nz = 0);
    void finalizeEnvironment ();

    int isActive ();
    int myRank ();
    int numProc ();
    void cartCoordToRank (const int & ix,
			  const int & iy,
			  const int & iz,
			  int & rank );
    void rankToCartCoord (const int & rank,
			  int & ix,
			  int & iy,
			  int & iz) ;
    void numProcDim (int & nx,
		     int & ny,
		     int & nz);
    void barrier ();
    void abort (int errorCode);

    unsigned numThreadsInCell  () ;

    void shiftNeighbor (int direction,
			int displacement,
			int & src,
			int & dest);
    
    // unsigned numThreadsForAtom ();
// namespace Transfer{
    //   void init ();
    //   void finalize ();
    //   void registSendBuff (void * buff, unsigned size);
    //   void registRecvBuff (void * buff, unsigned size);
    //   void clearRigistedSend ();
    //   void clearRigistedRecv ();
    //   void buildSend ();
    //   void buildRecv ();
    //   void Isend (int dest, int tag);
    //   void Irecv (int src,  int tag);
    //   bool testSend ();
    //   bool testRecv ();
    //   void waitSend ();
    //   void waitRecv ();
    // }
  }
}


    
#endif
