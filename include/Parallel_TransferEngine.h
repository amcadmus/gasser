#ifndef __Parallel_TransferEngine_h__
#define __Parallel_TransferEngine_h__

#define MPI_CODE
#include "Parallel_Environment.h"
#include "Parallel_DataTransferBlock.h"
#include "Parallel_CellList.h"
#include "Parallel_TransferEngineCompatible.h"
#include "mpi.h"
#include "common.h"


namespace Parallel{
  
  class TransferEngine 
  {
    int * blockLength;
    MPI_Aint * shiftIndex;
    void * buff;
    MPI_Aint add_buff;
    MPI_Datatype dataType;
    unsigned count;
    unsigned memSize;
    MPI_Request request;
    MPI_Status status;
private:
    void resize (unsigned memSize);
public:
    TransferEngine ();
    ~TransferEngine ();
public:
    void clear ();
    void clearRegistered ();
    void registerBuff (void * buff, size_t size);
    void registerBuff (TransferEngineCompatible & data);
    void build ();
    void Isend (int dest, int tag);
    void Irecv (int src,  int tag);
    void buildPersistentSend (int dest, int tag);
    void buildPersistentRecv (int src,  int tag);
    void startPersistentRequest ();
    bool test ();
    void wait ();
  };

//   class SummationEngine 
//   {
// private:
//     ScalorType * sumScalorBuff;
//     IndexType    sumScalorBuffSize;
//     IndexType * sumIndexBuff;
//     IndexType   sumIndexBuffSize;
//     IntScalorType * sumIntScalorBuff;
//     IntScalorType   sumIntScalorBuffSize;
// public:
//     SummationEngine ();
//     ~SummationEngine ();
// public:
//     ScalorType sumScalor    (ScalorType * data, int num, ScalorType ** result);
//     ScalorType sumScalorAll (ScalorType * data, int num, ScalorType ** result);
//     IntScalorType sumIntScalor    (IntScalorType * data, int num, IntScalorType **result);
//     IntScalorType sumIntScalorAll (IntScalorType * data, int num, IntScalorType **result);
//     IndexType sumIndex    (IndexType * data, int num, IndexType ** result);
//     IndexType sumIndexAll (IndexType * data, int num, IndexType ** result);
//   };

  class SummationEngine
  {
public:
    void sum (ScalorType * data,
	      int num,
	      ScalorType * result);
    void sum (IntScalorType * data,
	      int num,
	      IntScalorType * result);
    void sum (IndexType * data,
	      int num,
	      IndexType * result);    
    void sum (double * data,
	      int num,
	      double * result);    
    void sumAll (ScalorType * data,
		 int num,
		 ScalorType * result);
    void sumAll (IntScalorType * data,
		 int num,
		 IntScalorType * result);
    void sumAll (IndexType * data,
		 int num,
		 IndexType * result);  
    void sumAll (double * data,
		 int num,
		 double * result);  
  };
}

#endif
