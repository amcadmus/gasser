#include "Parallel_Timer.h"
// #include "common.h"
#include <stdlib.h>
#include "compile_error_mixcode.h"

using namespace Parallel::Timer;

bool      Parallel::Timer::hostTimerInited = false;
TimeType  Parallel::Timer::hostRecordReal [SizeOfTimeRecordArray] = {0};
TimeType  Parallel::Timer::hostRecordUser [SizeOfTimeRecordArray] = {0};
Stopwatch HostTimer::watch [SizeOfTimeRecordArray];
char     Parallel::Timer::hostWords[SizeOfTimeRecordArray][MaxWordsLength] = {{'\0'}};
char	 Parallel::Timer::hostSpace[SizeOfTimeRecordArray][MaxWordsLength];


float Parallel::Timer::
printDeviceItem (FILE * fp,
		 timeItem_t item)
{
  float tmp;
  fprintf(fp, "# %s:%s%1.3es   ( % 3.1f% )\n",
	  deviceWords[item],
	  deviceSpace[item],
	  deviceRecord[item] * 0.001,
	  tmp = (deviceRecord[item] * 100 * 0.001 /
		 hostRecordUser[item_Total - ParallelItemShift]));
  return tmp;
}

float Parallel::Timer::
printHostItem (FILE * fp,
	       timeItem_t item)
{
  float tmp;
  fprintf(fp, "# %s:%s%1.3es (User), %1.3es (Real)  (% 3.1f%, % 3.1f% )\n",
	  hostWords[item - ParallelItemShift],
	  hostSpace[item - ParallelItemShift],
	  hostRecordUser[item - ParallelItemShift],
	  hostRecordReal[item - ParallelItemShift],
	  tmp = (hostRecordUser[item - ParallelItemShift] * 100 /
		 hostRecordUser[item_Total - ParallelItemShift]),
	  hostRecordReal[item - ParallelItemShift] * 100 /
	  hostRecordUser[item_Total - ParallelItemShift]);
  return tmp;
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

  float totalPercent (0);
  totalPercent += printDeviceItem (fp, item_BuildCellList);
  totalPercent += printDeviceItem (fp, item_BuildBondList);
  totalPercent += printDeviceItem (fp, item_ApplyBondaryCondition);
  totalPercent += printDeviceItem (fp, item_RemoveTransFreedom);
  totalPercent += printDeviceItem (fp, item_NonBondedInteraction);
  totalPercent += printDeviceItem (fp, item_NonBondedInterStatistic);
  totalPercent += printDeviceItem (fp, item_BondedInteraction);
  totalPercent += printDeviceItem (fp, item_BondedInterStatistic);
  totalPercent += printDeviceItem (fp, item_Integrate);
  totalPercent += printDeviceItem (fp, item_DataIO);

  totalPercent += printHostItem (fp, item_Redistribute);
  printHostItem (fp, item_Redistribute_Data);
  printHostItem (fp, item_Redistribute_Data0);
  printHostItem (fp, item_Redistribute_Transfer);
  printHostItem (fp, item_Redistribute_DHCopy);
  totalPercent += printHostItem (fp, item_TransferGhost);

  fprintf (fp, "# Total percentage:");
  for (int i = 0; i < PrintStartPosition - 16; ++i){
    fprintf (fp, " ");
  }
  fprintf (fp, "%3.1f \%\n", totalPercent);
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
  strncpy (hostWords[item_Redistribute_Data0 - ParallelItemShift] ,
	   "Redistribute atoms build data struct 0",
	   MaxWordsLength);
  strncpy (hostWords[item_Redistribute_Transfer - ParallelItemShift] ,
	   "Redistribute atoms transfer",
	   MaxWordsLength);
  strncpy (hostWords[item_Redistribute_DHCopy - ParallelItemShift] ,
	   "Redistribute atoms device-host copy",
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
