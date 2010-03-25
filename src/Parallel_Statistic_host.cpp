#define HOST_CODE

#include "Parallel_Statistic.h"
#include "Parallel_TransferEngine.h"
#include "compile_error_mixcode.h"

Parallel::HostStatistic::
HostStatistic ()
    : localData (NULL), globalData (NULL)
{
}

Parallel::HostStatistic::
~HostStatistic ()
{
  freeAPointer ((void**)&localData);
  freeAPointer ((void**)&globalData);
}

void Parallel::HostStatistic::
reinit (const HostCellListedMDData & sys)
{
  volume = sys.getGlobalBox().size.x *
      sys.getGlobalBox().size.y * sys.getGlobalBox().size.z;
  volumei = 1./volume;
  
  size = sizeof(ScalorType) * NumberOfStatisticItems;
  localData  = (ScalorType *) malloc (size);
  if (localData == NULL){
    throw MDExcptFailedMallocOnHost ("HostStatistic::reinit", "localData", size);
  }
  globalData = (ScalorType *) malloc (size);
  if (globalData == NULL){
    throw MDExcptFailedMallocOnHost ("MDStatistic::MDStatistic", "hdata", size);
  }

#pragma unroll NumberOfStatisticItems
  for (IndexType i = 0; i < NumberOfStatisticItems; ++i){
    localData[i] = 0.f;
    globalData[i] = 0.f;
  }
}


void Parallel::HostStatistic::
collectData ()
{
  Parallel::SummationEngine eng;
  eng.sum (localData, NumberOfStatisticItems, globalData);
}

void Parallel::HostStatistic::
collectDataAll ()
{
  Parallel::SummationEngine eng;
  eng.sumAll (localData, NumberOfStatisticItems, globalData);
}

void Parallel::HostStatistic::
collectData (const mdStatisticItem_t item)
{
  Parallel::SummationEngine eng;
  eng.sum (&localData[item], 1, &globalData[item]);
}

void Parallel::HostStatistic::
collectDataAll (const mdStatisticItem_t item)
{
  Parallel::SummationEngine eng;
  eng.sumAll (&localData[item], 1, &globalData[item]);
}

