#define DEVICE_CODE
#include "Parallel_Timer.h"
#include "compile_error_mixcode.h"

using namespace Parallel::Timer;

bool     Parallel::Timer::deviceTimerInited = false;
TimeType Parallel::Timer::deviceRecord[SizeOfTimeRecordArray];
char     Parallel::Timer::deviceWords[SizeOfTimeRecordArray][MaxWordsLength] = {{'\0'}};
char	 Parallel::Timer::deviceSpace[SizeOfTimeRecordArray][MaxWordsLength];

cudaEvent_t	Parallel::Timer::DeviceTimer::start [SizeOfTimeRecordArray];
cudaEvent_t	Parallel::Timer::DeviceTimer::stop  [SizeOfTimeRecordArray];

void Parallel::Timer::DeviceTimer::
init ()
{
  for (IndexType i = 0; i < SizeOfTimeRecordArray; ++i){
    cudaEventCreate (&start[i]);
    cudaEventCreate (&stop [i]);
    deviceRecord[i] = 0.f;
  }
  strncpy (deviceWords[item_ApplyBondaryCondition] ,
	   "Apply bondary condition",
	   MaxWordsLength);
  strncpy (deviceWords[item_BuildCellList],
	   "Build cell list",
	   MaxWordsLength);
  strncpy (deviceWords[item_Integrate],
	   "Integrate",
	   MaxWordsLength);
  strncpy (deviceWords[item_RemoveTransFreedom],
	   "Remove translational freedom",
	   MaxWordsLength);
  strncpy (deviceWords[item_NonBondedInteraction],
	   "Non bonded interaction",
	   MaxWordsLength);
  strncpy (deviceWords[item_NonBondedInterStatistic],
	   "Non bonded interaction with st",
	   MaxWordsLength);
  strncpy (deviceWords[item_BondedInteraction],
	   "Bonded interaction",
	   MaxWordsLength);
  strncpy (deviceWords[item_BondedInterStatistic],
	   "Bonded interaction with st",
	   MaxWordsLength);
  strncpy (deviceWords[item_BuildBondList],
	   "Build bond list",
	   MaxWordsLength);
  strncpy (deviceWords[item_DataIO],
	   "Data IO",
	   MaxWordsLength);
  strncpy (deviceWords[item_ClearInteraction],
	   "Clear interaction",
	   MaxWordsLength);
  
  for (IndexType i = 0; i < SizeOfTimeRecordArray; ++i){
    int wordLength = strlen (deviceWords[i]);
    int nspace = PrintStartPosition - wordLength;
    if (nspace < 0){
      nspace = 0;
    }
    int j = 0;
    for (; j < nspace; ++j){
      deviceSpace[i][j] = ' ';
    }
    deviceSpace[i][j] ='\0';
  }
      
  deviceTimerInited = true;
}

void Parallel::Timer::DeviceTimer::
finalize ()
{
  for (IndexType i = 0; i < SizeOfTimeRecordArray; ++i){
    cudaEventDestroy (start[i]);
    cudaEventDestroy (stop [i]);
  }
  deviceTimerInited = false;
}


void Parallel::Timer::DeviceTimer::
reset ()
{
  for (IndexType i = 0; i < SizeOfTimeRecordArray; ++i){
    deviceRecord[i] = 0.f;
  }
}

