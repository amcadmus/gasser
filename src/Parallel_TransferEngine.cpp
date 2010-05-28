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
      count (0),
      memSize (0),
      request (MPI_REQUEST_NULL),
      dataType (MPI_DATATYPE_NULL)
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
  if (dataType != MPI_DATATYPE_NULL) MPI_Type_free (&dataType);
  MPI_Type_hindexed (count, blockLength, shiftIndex, MPI_BYTE, &dataType);
  MPI_Type_commit (&dataType);
}

void Parallel::TransferEngine::
Isend (int dest, int tag)
{
  if (request != MPI_REQUEST_NULL) MPI_Request_free (&request);
  MPI_Isend (buff, 1, dataType, dest, tag, Parallel::Environment::communicator(), &request);
}

void Parallel::TransferEngine::
Irecv (int src,  int tag)
{
  if (request != MPI_REQUEST_NULL) MPI_Request_free (&request);
  MPI_Irecv (buff, 1, dataType, src,  tag, Parallel::Environment::communicator(), &request);
}

void Parallel::TransferEngine::
buildPersistentSend (int dest, int tag)
{
  if (dataType != MPI_DATATYPE_NULL) MPI_Type_free (&dataType);
  MPI_Type_hindexed (count, blockLength, shiftIndex, MPI_BYTE, &dataType);
  MPI_Type_commit (&dataType);
  if (request != MPI_REQUEST_NULL) MPI_Request_free (&request);
  MPI_Send_init (buff, 1, dataType, dest, tag, Parallel::Environment::communicator(),
		 &request);
}

void Parallel::TransferEngine::
buildPersistentRecv (int src,  int tag)
{
  if (dataType != MPI_DATATYPE_NULL) MPI_Type_free (&dataType);
  MPI_Type_hindexed (count, blockLength, shiftIndex, MPI_BYTE, &dataType);
  MPI_Type_commit (&dataType);
  if (request != MPI_REQUEST_NULL) MPI_Request_free (&request);
  MPI_Recv_init (buff, 1, dataType, src, tag, Parallel::Environment::communicator(),
		 &request);
}

void Parallel::TransferEngine::
startPersistentRequest ()
{
  MPI_Start (&request);
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
  if (request != MPI_REQUEST_NULL) MPI_Request_free (&request);
  if (dataType != MPI_DATATYPE_NULL) MPI_Type_free (&dataType);
  count = memSize = 0;
}

void Parallel::TransferEngine::
clearRegistered ()
{
  count = 0;
  buff = NULL;
  if (request != MPI_REQUEST_NULL) MPI_Request_free (&request);
  if (dataType != MPI_DATATYPE_NULL) MPI_Type_free (&dataType);
}

Parallel::TransferEngine::
~TransferEngine ()
{
  // clearRegistered ();
  // freeAPointer ((void**)&blockLength);
  // freeAPointer ((void**)&shiftIndex);
  clear ();
}



// // Parallel::SummationEngine::
// // SummationEngine ()
// //     : sumScalorBuff (NULL),
// //       sumScalorBuffSize (0),
// //       sumIndexBuff (NULL),
// //       sumIndexBuffSize (0),
// //       sumIntScalorBuff (NULL),
// //       sumIntScalorBuffSize (0)
// // {
// // }

// // Parallel::SummationEngine::
// // ~SummationEngine ()
// // {
// //   freeAPointer ((void**)&sumScalorBuff);
// //   freeAPointer ((void**)&sumIntScalorBuff);
// //   freeAPointer ((void**)&sumIndexBuff);
// // }

// // ScalorType Parallel::SummationEngine::
// // sumScalor (ScalorType * data, int num, ScalorType ** result)
// // {
// //   Environment env;
// //   if (num > sumScalorBuffSize && env.myRank() == 0){
// //     sumScalorBuff = (ScalorType *)realloc(sumScalorBuff, num * sizeof(ScalorType));
// //     if (sumScalorBuff == NULL) throw (MDExcptFailedReallocOnHost("SummationEngine::result", sizeof(ScalorType)*num));
// //     sumScalorBuffSize = num;
// //   }
  
// //   MPI_Reduce (data, sumScalorBuff, num, MPI_FLOAT, MPI_SUM, 0, Parallel::Environment::communicator());;
// //   *result = sumScalorBuff;
// // }

void Parallel::SummationEngine::
sum (ScalorType * data, int num, ScalorType * result)
{
  MPI_Reduce (data, result, num, MPI_FLOAT, MPI_SUM, 0,
	      Parallel::Environment::communicator());
}

void Parallel::SummationEngine::
sum (IntScalorType * data, int num, IntScalorType * result)
{
  MPI_Reduce (data, result, num, MPI_INT, MPI_SUM, 0,
	      Parallel::Environment::communicator());
}

void Parallel::SummationEngine::
sum (IndexType * data, int num, IndexType * result)
{
  MPI_Reduce (data, result, num, MPI_UNSIGNED, MPI_SUM, 0,
	      Parallel::Environment::communicator());
}

void Parallel::SummationEngine::
sum (double * data, int num, double * result)
{
  MPI_Reduce (data, result, num, MPI_DOUBLE, MPI_SUM, 0,
	      Parallel::Environment::communicator());
}

void Parallel::SummationEngine::
sumAll (ScalorType * data, int num, ScalorType * result)
{
  MPI_Allreduce (data, result, num, MPI_FLOAT, MPI_SUM, 
		 Parallel::Environment::communicator());
}

void Parallel::SummationEngine::
sumAll (IntScalorType * data, int num, IntScalorType * result)
{
  MPI_Allreduce (data, result, num, MPI_INT, MPI_SUM, 
		 Parallel::Environment::communicator());
}

void Parallel::SummationEngine::
sumAll (IndexType * data, int num, IndexType * result)
{
  MPI_Allreduce (data, result, num, MPI_UNSIGNED, MPI_SUM, 
		 Parallel::Environment::communicator());
}

void Parallel::SummationEngine::
sumAll (double * data, int num, double * result)
{
  MPI_Allreduce (data, result, num, MPI_DOUBLE, MPI_SUM, 
		 Parallel::Environment::communicator());
}



// // ScalorType Parallel::SummationEngine::
// // sumScalorAll (ScalorType * data, int num, ScalorType ** result)
// // {
// //   if (num > sumScalorBuffSize){
// //     sumScalorBuff = (ScalorType *)realloc(sumScalorBuff, num * sizeof(ScalorType));
// //     if (sumScalorBuff == NULL) throw (MDExcptFailedReallocOnHost("SummationEngine::result", sizeof(ScalorType)*num));
// //     sumScalorBuffSize = num;
// //   }
  
// //   MPI_Allreduce (data, sumScalorBuff, num, MPI_FLOAT, MPI_SUM, Parallel::Environment::communicator());;
// //   *result = sumScalorBuff;
// // }

// // IndexType Parallel::SummationEngine::
// // sumIndex (IndexType * data, int num, IndexType ** result)
// // {
// //   Environment env;
// //   if (num > sumIndexBuffSize && env.myRank() == 0){
// //     sumIndexBuff = (IndexType *)realloc(sumIndexBuff, num * sizeof(IndexType));
// //     if (sumIndexBuff == NULL) throw (MDExcptFailedReallocOnHost("SummationEngine::result", sizeof(IndexType)*num));
// //     sumIndexBuffSize = num;
// //   }
  
// //   MPI_Reduce (data, sumIndexBuff, num, MPI_UNSIGNED, MPI_SUM, 0, Parallel::Environment::communicator());;
// //   *result = sumIndexBuff;
// // }


// // IndexType Parallel::SummationEngine::
// // sumIndexAll (IndexType * data, int num, IndexType ** result)
// // {
// //   if (num > sumIndexBuffSize){
// //     sumIndexBuff = (IndexType *)realloc(sumIndexBuff, num * sizeof(IndexType));
// //     if (sumIndexBuff == NULL) throw (MDExcptFailedReallocOnHost("SummationEngine::result", sizeof(IndexType)*num));
// //     sumIndexBuffSize = num;
// //   }
// //   // for (int i = 0; i < num; ++i){
// //   //   sumIndexBuff[i] = 0;
// //   // }
  
// //   MPI_Allreduce (data, sumIndexBuff, num, MPI_UNSIGNED, MPI_SUM, Parallel::Environment::communicator());;
// //   *result = sumIndexBuff;
// // }


// void Parallel::TransferEngine::
// registerBuff (const DataTransferBlock & block)
// {
//   for (IndexType i = 0; i < NumDataItemsInMDData; ++i){
//     registerBuff (block.pointer[i], block.size[i]);
//   }
// }

// void Parallel::TransferEngine::
// registerBuff (HostSubCellList &hsubCell, const MDDataItemMask_t mask)
// {
//   IndexType num;
//   void ** buffs = NULL;
//   size_t * sizes = NULL;
//   hsubCell.collectBuffInfo (mask, &num, &buffs, &sizes);

//   for (IndexType i = 0; i < num; ++i){
//     registerBuff (buffs[i], sizes[i]);
//   }
// }

void Parallel::TransferEngine::
registerBuff (TransferEngineCompatible & data)
{
  IndexType num;
  void ** buffs;
  size_t * sizes;

  data.getTransBuffs (&num, &buffs, &sizes);

  for (unsigned i = 0; i < num; ++i){
    registerBuff (buffs[i], sizes[i]);
  }
}

