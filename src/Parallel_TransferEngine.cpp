#define MPI_CODE

#include "Parallel_TransferEngine.h"
#include "Parallel_Environment.h"

#include "common.h"
#include <stdlib.h>

#include "compile_error_mixcode.h"

Parallel::TransferEngine::
TransferEngine ()
    : buff (NULL),
      blockLength(NULL),
      shiftIndex(NULL),
      dataTypeBuilt (false),
      count (0),
      memSize (0)
{
}

void Parallel::TransferEngine::
resize (unsigned memSize_)
{
  memSize = memSize_;
  blockLength = (int *) realloc (blockLength, sizeof(int) * memSize);
  if (blockLength == NULL) throw (MDExcptFailedReallocOnHost("blockLength", sizeof(int)*memSize));  
  shiftIndex = (MPI_Aint *) realloc (shiftIndex, sizeof(MPI_Aint) * memSize);
  if (shiftIndex == NULL) throw (MDExcptFailedReallocOnHost("shiftIndex", sizeof(MPI_Aint)*memSize));
}

void Parallel::TransferEngine::
registerBuff (void * ptr, size_t size)
{
  if (count == memSize){
    memSize ++;
    memSize <<= 1;
  }
  resize (memSize);

  if (count == 0){
    buff = ptr;
    MPI_Address (buff, &add_buff);
  }
  blockLength[count] = size;
  MPI_Aint add_ptr;
  MPI_Address (ptr, &add_ptr);
  shiftIndex [count] = (add_ptr - add_buff);
  count ++;
}

void Parallel::TransferEngine::
build ()
{
  MPI_Type_hindexed (count, blockLength, shiftIndex, MPI_BYTE, &dataType);
  MPI_Type_commit (&dataType);
  dataTypeBuilt = true;
}

void Parallel::TransferEngine::
Isend (int dest, int tag)
{
  MPI_Isend (buff, 1, dataType, dest, tag, Parallel::Environment::communicator(), &request);
}

void Parallel::TransferEngine::
Irecv (int src,  int tag)
{
  MPI_Irecv (buff, 1, dataType, src,  tag, Parallel::Environment::communicator(), &request);
}

bool Parallel::TransferEngine::
test ()
{
  int flag;
  MPI_Test (&request, &flag, &status);
  return flag != 0;
}

void Parallel::TransferEngine::
wait ()
{
  MPI_Wait (&request, &status);
}

void Parallel::TransferEngine::
clear ()
{
  freeAPointer ((void**)&blockLength);
  freeAPointer ((void**)&shiftIndex);
  buff = NULL;
  if (dataTypeBuilt){
    MPI_Type_free (&dataType);
    dataTypeBuilt = false;
  }
  count = memSize = 0;
}

void Parallel::TransferEngine::
clearRegistered ()
{
  count = 0;
  buff = NULL;
  if (dataTypeBuilt){
    MPI_Type_free (&dataType);
    dataTypeBuilt = false;
  }
}

Parallel::TransferEngine::
~TransferEngine ()
{
  freeAPointer ((void**)&blockLength);
  freeAPointer ((void**)&shiftIndex);
}


Parallel::SummationEngine::
SummationEngine ()
    : sumScalorBuff (NULL),
      sumScalorBuffSize (0),
      sumIndexBuff (NULL),
      sumIndexBuffSize (0),
      sumIntScalorBuff (NULL),
      sumIntScalorBuffSize (0)
{
}

Parallel::SummationEngine::
~SummationEngine ()
{
  freeAPointer ((void**)&sumScalorBuff);
  freeAPointer ((void**)&sumIntScalorBuff);
  freeAPointer ((void**)&sumIndexBuff);
}

ScalorType Parallel::SummationEngine::
sumScalor (ScalorType * data, int num, ScalorType ** result)
{
  Environment env;
  if (num > sumScalorBuffSize && env.myRank() == 0){
    sumScalorBuff = (ScalorType *)realloc(sumScalorBuff, num * sizeof(ScalorType));
    if (sumScalorBuff == NULL) throw (MDExcptFailedReallocOnHost("SummationEngine::result", sizeof(ScalorType)*num));
    sumScalorBuffSize = num;
  }
  
  MPI_Reduce (data, sumScalorBuff, num, MPI_FLOAT, MPI_SUM, 0, Parallel::Environment::communicator());;
  *result = sumScalorBuff;
}


ScalorType Parallel::SummationEngine::
sumScalorAll (ScalorType * data, int num, ScalorType ** result)
{
  if (num > sumScalorBuffSize){
    sumScalorBuff = (ScalorType *)realloc(sumScalorBuff, num * sizeof(ScalorType));
    if (sumScalorBuff == NULL) throw (MDExcptFailedReallocOnHost("SummationEngine::result", sizeof(ScalorType)*num));
    sumScalorBuffSize = num;
  }
  
  MPI_Allreduce (data, sumScalorBuff, num, MPI_FLOAT, MPI_SUM, Parallel::Environment::communicator());;
  *result = sumScalorBuff;
}

IndexType Parallel::SummationEngine::
sumIndex (IndexType * data, int num, IndexType ** result)
{
  Environment env;
  if (num > sumIndexBuffSize && env.myRank() == 0){
    sumIndexBuff = (IndexType *)realloc(sumIndexBuff, num * sizeof(IndexType));
    if (sumIndexBuff == NULL) throw (MDExcptFailedReallocOnHost("SummationEngine::result", sizeof(IndexType)*num));
    sumIndexBuffSize = num;
  }
  
  MPI_Reduce (data, sumIndexBuff, num, MPI_UNSIGNED, MPI_SUM, 0, Parallel::Environment::communicator());;
  *result = sumIndexBuff;
}


IndexType Parallel::SummationEngine::
sumIndexAll (IndexType * data, int num, IndexType ** result)
{
  if (num > sumIndexBuffSize){
    sumIndexBuff = (IndexType *)realloc(sumIndexBuff, num * sizeof(IndexType));
    if (sumIndexBuff == NULL) throw (MDExcptFailedReallocOnHost("SummationEngine::result", sizeof(IndexType)*num));
    sumIndexBuffSize = num;
  }
  // for (int i = 0; i < num; ++i){
  //   sumIndexBuff[i] = 0;
  // }
  
  MPI_Allreduce (data, sumIndexBuff, num, MPI_UNSIGNED, MPI_SUM, Parallel::Environment::communicator());;
  *result = sumIndexBuff;
}

