#ifndef __Parallel_interface_h__
#define __Parallel_interface_h__

namespace Parallel{
  namespace Interface{
    void initEnvironment (int * argc, char *** argv);
    void finalizeEnvironment ();
    void initCart (const int & nx,
		   const int & ny,
		   const int & nz);
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

    unsigned numThreadsInCell  () ;
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
