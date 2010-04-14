#include "Parallel_Timer.h"
#include <stdlib.h>
#include "compile_error_mixcode.h"

using namespace Parallel::Timer;

bool      Parallel::Timer::hostTimerInited = false;
TimeType  Parallel::Timer::hostRecordReal [SizeOfTimeRecordArray] = {0};
TimeType  Parallel::Timer::hostRecordUser [SizeOfTimeRecordArray] = {0};
Stopwatch HostTimer::watch [SizeOfTimeRecordArray];
char     Parallel::Timer::hostWords[SizeOfTimeRecordArray][MaxWordsLength] = {{'\0'}};
char	 Parallel::Timer::hostSpace[SizeOfTimeRecordArray][MaxWordsLength];


void Parallel::Timer::
printDeviceItem (FILE * fp,
		 timeItem_t item)
{
  fprintf(fp, "# %s:%s%1.3es   ( % 3.1f% )\n",
	  deviceWords[item],
	  deviceSpace[item],
	  deviceRecord[item] * 0.001,
	  deviceRecord[item] * 100 * 0.001 /
	  hostRecordUser[item_Total - ParallelItemShift]);
}

void Parallel::Timer::
printHostItem (FILE * fp,
	       timeItem_t item)
{
  fprintf(fp, "# %s:%s%1.3es (User), %1.3es (Real)  (% 3.1f%, % 3.1f% )\n",
	  hostWords[item - ParallelItemShift],
	  hostSpace[item - ParallelItemShift],
	  hostRecordUser[item - ParallelItemShift],
	  hostRecordReal[item - ParallelItemShift],
	  hostRecordUser[item - ParallelItemShift] * 100 /
	  hostRecordUser[item_Total - ParallelItemShift],
	  hostRecordReal[item - ParallelItemShift] * 100 /
	  hostRecordUser[item_Total - ParallelItemShift]);
}

void Parallel::Timer::
printTotalTime (FILE * fp)
{
  fprintf (fp, "# %s:%s%1.3es (User), %1.3es (Real)\n",
	   hostWords[item_Total - ParallelItemShift],
	   hostSpace[item_Total - ParallelItemShift],	   
	   hostRecordUser[item_Total - ParallelItemShift],
	   hostRecordReal[item_Total - ParallelItemShift]);
}


void Parallel::Timer::
printRecord (FILE * fp)
{
  if (!hostTimerInited && !deviceTimerInited) return;
  if (hostRecordUser[item_Total - ParallelItemShift] == 0) return;

  printDeviceItem (fp, item_BuildCellList);
  printDeviceItem (fp, item_ApplyBondaryCondition);
  printDeviceItem (fp, item_NonBondedInteraction);
  printDeviceItem (fp, item_NonBondedInterStatistic);
  printDeviceItem (fp, item_Integrate);

  printHostItem (fp, item_Redistribute);
  printHostItem (fp, item_Redistribute_Data);
  printHostItem (fp, item_Redistribute_Transfer);
  printHostItem (fp, item_TransferGhost);

  printTotalTime (fp);
}

void Parallel::Timer::HostTimer::
tic (timeItem_t item)
{
  watch[item - ParallelItemShift].start();
}

void Parallel::Timer::HostTimer::
toc (timeItem_t item)
{
  watch[item - ParallelItemShift].stop();
  hostRecordReal[item - ParallelItemShift] +=
      watch[item - ParallelItemShift].real();
  hostRecordUser[item - ParallelItemShift] +=
      watch[item - ParallelItemShift].user();
}

void Parallel::Timer::HostTimer::
init ()
{
  strncpy (hostWords[item_Total - ParallelItemShift] ,
	   "Total time comsuption is",
	   MaxWordsLength);
  strncpy (hostWords[item_Redistribute - ParallelItemShift] ,
	   "Redistribute atoms",
	   MaxWordsLength);
  strncpy (hostWords[item_Redistribute_Data - ParallelItemShift] ,
	   "Redistribute atoms build data struct",
	   MaxWordsLength);
  strncpy (hostWords[item_Redistribute_Transfer - ParallelItemShift] ,
	   "Redistribute atoms transfer",
	   MaxWordsLength);
  strncpy (hostWords[item_TransferGhost - ParallelItemShift] ,
	   "Transfer Ghosts",
	   MaxWordsLength);

  for (IndexType i = 0; i < SizeOfTimeRecordArray; ++i){
    int wordLength = strlen (hostWords[i]);
    int nspace = PrintStartPosition - wordLength;
    if (nspace < 0){
      nspace = 0;
    }
    int j = 0;
    for (; j < nspace; ++j){
      hostSpace[i][j] = ' ';
    }
    hostSpace[i][j] ='\0';
  }

  hostTimerInited = true;
}

void Parallel::Timer::HostTimer::
finalize()
{
  hostTimerInited = false;
}
