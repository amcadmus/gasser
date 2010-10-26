#ifndef __Parallel_Timer_h_wanghan__
#define __Parallel_Timer_h_wanghan__

#include <stdio.h>
#include "common.h"
#include "Stopwatch.h"


namespace Parallel{
  namespace Timer{  
    typedef float TimeType;
    const IndexType ParallelItemShift		= 10000;
    const IndexType SizeOfTimeRecordArray	= 32;
    const IndexType MaxWordsLength		= 1024;
    const int       PrintStartPosition		= 64 - 3;
    
    extern TimeType  deviceRecord	[SizeOfTimeRecordArray];
    extern TimeType  hostRecordReal	[SizeOfTimeRecordArray];
    extern TimeType  hostRecordUser	[SizeOfTimeRecordArray];
    extern char      deviceWords	[SizeOfTimeRecordArray][MaxWordsLength];
    extern char      deviceSpace	[SizeOfTimeRecordArray][MaxWordsLength];
    extern char      hostWords		[SizeOfTimeRecordArray][MaxWordsLength];
    extern char      hostSpace		[SizeOfTimeRecordArray][MaxWordsLength];
    
    extern bool hostTimerInited;
    extern bool deviceTimerInited;
  
    enum timeItem {
      item_Total				= 0 + ParallelItemShift,
      
      item_ApplyBondaryCondition		= 1,
      item_BuildCellList			= 2,
      item_Integrate				= 5,
      item_RemoveTransFreedom			= 6,
      item_NonBondedInteraction			= 7,
      item_NonBondedInterStatistic		= 8,
      item_BondedInteraction			= 9,
      item_BondedInterStatistic			= 10,
      item_AngleInteraction			= 11,
      item_AngleInterStatistic			= 12,
      item_BuildBondList			= 13,
      item_ClearInteraction			= 14,
    
      item_DataTransfer				= 19,
      item_DataIO				= 20,

      item_Redistribute				= 1 + ParallelItemShift,
      item_Redistribute_SyncNum			= 3 + ParallelItemShift,
      item_Redistribute_BuildEngine		= 6 + ParallelItemShift,
      item_Redistribute_Transfer		= 4 + ParallelItemShift,
      item_Redistribute_DHCopy			= 5 + ParallelItemShift,
      item_Redistribute_Barrier			= 7 + ParallelItemShift,
      item_TransferGhost			= 2 + ParallelItemShift,
      item_TransferGhost_Tranfer		= 9 + ParallelItemShift,
      item_TransferGhost_DHCopy			= 8 + ParallelItemShift,
      item_TransferGhost_DHCopy_Pack		= 10 + ParallelItemShift,
      item_TransferGhost_DHCopy_Unpack		= 11 + ParallelItemShift,
      item_TransferGhost_DHCopy_Copy		= 12 + ParallelItemShift,
    };
    typedef enum timeItem timeItem_t;

    void  printRecord (FILE * fp);
    float printDeviceItem (FILE * fp, timeItem_t item);
    float printHostItem   (FILE * fp, timeItem_t item);
    void  printTotalTime  (FILE * fp);

    namespace HostTimer {
      extern Stopwatch watch [SizeOfTimeRecordArray];
      void init ();
      void finalize ();
      void reset ();
      void tic (timeItem_t item);
      void toc (timeItem_t item);
    }
    
      
#ifdef DEVICE_CODE
    namespace DeviceTimer{
      extern cudaEvent_t start [SizeOfTimeRecordArray];
      extern cudaEvent_t stop  [SizeOfTimeRecordArray];
      void init ();
      void finalize ();
      void reset ();
      void tic (timeItem_t timeItem,
		cudaStream_t stream = 0);
      TimeType toc (timeItem_t timeItem,
		    cudaStream_t stream = 0);
    }
#endif
  }
}





#ifdef DEVICE_CODE
inline void Parallel::Timer::DeviceTimer::
tic(Parallel::Timer::timeItem_t item,
    cudaStream_t stream)
{
  cudaEventRecord(start[item], stream);
}
#endif

#ifdef DEVICE_CODE
inline Parallel::Timer::TimeType Parallel::Timer::DeviceTimer::
toc(Parallel::Timer::timeItem_t item,
    cudaStream_t stream)
{
  TimeType tmptime;
  
  cudaEventRecord(stop[item], stream);
  cudaEventSynchronize (stop[item]);
  cudaEventElapsedTime (&tmptime, start[item], stop[item]);

  Parallel::Timer::deviceRecord[item] += tmptime;
  return tmptime;
}
#endif


#endif
