#define DEVICE_CODE
#include "Parallel_Timer.h"
#include "compile_error_mixcode.h"

bool
Parallel::Timer::deviceTimerInited = false;
Parallel::Timer::TimeType
Parallel::Timer::deviceRecord[SizeOfTimeRecordArray];
char
Parallel::Timer::deviceWords[SizeOfTimeRecordArray][MaxWordsLength] = {{'\0'}};
char		Parallel::Timer::deviceSpace[SizeOfTimeRecordArray][MaxWordsLength];
cudaEvent_t	Parallel::Timer::DeviceTimer::start [SizeOfTimeRecordArray];
cudaEvent_t	Parallel::Timer::DeviceTimer::stop  [SizeOfTimeRecordArray];

void Parallel::Timer::DeviceTimer::
init ()
{
  for (IndexType i = 0; i < SizeOfTimeRecordArray; ++i){
    cudaEventCreate (&start[i]);
    cudaEventCreate (&stop [i]);
    Parallel::Timer::deviceRecord[i] = 0.f;
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
  for (IndexType i = 0; i < SizeOfTimeRecordArray; ++i){
    int wordLength = strlen (deviceWords[i]);
    int nspace = PrintStartPosition - wordLength;
    if (nspace < 1){
      nspace = 1;
    }
    int j = 0;
    Parallel::Timer::deviceSpace[i][j++] = ':';
    for (; j < nspace; ++j){
      Parallel::Timer::deviceSpace[i][j] = ' ';
    }
    Parallel::Timer::deviceSpace[i][j] ='\0';
  }
      
  Parallel::Timer::deviceTimerInited = true;
}

void Parallel::Timer::DeviceTimer::
finalize ()
{
  for (IndexType i = 0; i < SizeOfTimeRecordArray; ++i){
    cudaEventDestroy (start[i]);
    cudaEventDestroy (stop [i]);
  }
  Parallel::Timer::deviceTimerInited = false;
}


void Parallel::Timer::DeviceTimer::
reset ()
{
  for (IndexType i = 0; i < SizeOfTimeRecordArray; ++i){
    Parallel::Timer::deviceRecord[i] = 0.f;
  }
}

