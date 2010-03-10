#ifndef __Parallel_TransferEngine_h__
#define __Parallel_TransferEngine_h__

#define MPI_CODE
#include "Parallel_Environment.h"
#include "Parallel_DataTransferBlock.h"
#include "mpi.h"
#include "common.h"


namespace Parallel{

  class TransferEngine 
  {
    int * blockLength;
    MPI_Aint * shiftIndex;
    void * buff;
    MPI_Aint add_buff;
    bool dataTypeBuilt;
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
    void registerBuff (const DataTransferBlock & block);
    void build ();
    void Isend (int dest, int tag);
    void Irecv (int src,  int tag);
    bool test ();
    void wait ();
  };

  class SummationEngine 
  {
private:
    ScalorType * sumScalorBuff;
    IndexType    sumScalorBuffSize;
    IndexType * sumIndexBuff;
    IndexType   sumIndexBuffSize;
    IntScalorType * sumIntScalorBuff;
    IntScalorType   sumIntScalorBuffSize;
public:
    SummationEngine ();
    ~SummationEngine ();
public:
    ScalorType sumScalor    (ScalorType * data, int num, ScalorType ** result);
    ScalorType sumScalorAll (ScalorType * data, int num, ScalorType ** result);
    IntScalorType sumIntScalor    (IntScalorType * data, int num, IntScalorType **result);
    IntScalorType sumIntScalorAll (IntScalorType * data, int num, IntScalorType **result);
    IndexType sumIndex    (IndexType * data, int num, IndexType ** result);
    IndexType sumIndexAll (IndexType * data, int num, IndexType ** result);
  };
  
}



#endif
