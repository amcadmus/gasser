#include "Parallel_Timer.h"
#include "Parallel_Interface.h"
// #include "common.h"
#include <stdlib.h>
#include <stdio.h>
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

  for (int i = 0; i < Parallel::Interface::numProc(); ++i){
    if (i == Parallel::Interface::myRank()){
      for (unsigned j = 0; j < PrintStartPosition + 32; ++j){
	// fprintf (fp, "##################################################\n");
	fputc ('#', fp);
      }
      fputc ('\n', fp);
      fprintf (fp, "# Summary for efficiency, rank %d\n", i);      
      float totalPercent (0);
      totalPercent += printDeviceItem (fp, item_BuildCellList);
      totalPercent += printDeviceItem (fp, item_BuildBondList);
      totalPercent += printDeviceItem (fp, item_ApplyBondaryCondition);
      totalPercent += printDeviceItem (fp, item_RemoveTransFreedom);
      totalPercent += printDeviceItem (fp, item_ClearInteraction);
      totalPercent += printDeviceItem (fp, item_NonBondedInteraction);
      totalPercent += printDeviceItem (fp, item_NonBondedInterStatistic);
      totalPercent += printDeviceItem (fp, item_BondedInteraction);
      totalPercent += printDeviceItem (fp, item_BondedInterStatistic);
      totalPercent += printDeviceItem (fp, item_Integrate);
      totalPercent += printDeviceItem (fp, item_DataIO);

      totalPercent += printHostItem (fp, item_Redistribute);
      printHostItem (fp, item_Redistribute_SyncNum);
      printHostItem (fp, item_Redistribute_BuildEngine);
      printHostItem (fp, item_Redistribute_Transfer);
      printHostItem (fp, item_Redistribute_DHCopy);
      printHostItem (fp, item_Redistribute_Barrier);
      totalPercent += printHostItem (fp, item_TransferGhost);
      printHostItem (fp, item_TransferGhost_Tranfer);
      printHostItem (fp, item_TransferGhost_DHCopy);
      printHostItem (fp, item_TransferGhost_DHCopy_Pack);
      printHostItem (fp, item_TransferGhost_DHCopy_Copy);
      printHostItem (fp, item_TransferGhost_DHCopy_Unpack);

      fprintf (fp, "# Total percentage:");
      for (int i = 0; i < PrintStartPosition - 16; ++i){
	fprintf (fp, " ");
      }
      fprintf (fp, "%3.1f \%\n", totalPercent);
      printTotalTime (fp);
    }
    Parallel::Interface::barrier();
  }
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
  strncpy (hostWords[item_Redistribute_SyncNum - ParallelItemShift] ,
	   "Redistribute atoms sync number to trans",
	   MaxWordsLength);
  strncpy (hostWords[item_Redistribute_BuildEngine - ParallelItemShift] ,
	   "Redistribute atoms build engine",
	   MaxWordsLength);
  strncpy (hostWords[item_Redistribute_Transfer - ParallelItemShift] ,
	   "Redistribute atoms transfer",
	   MaxWordsLength);
  strncpy (hostWords[item_Redistribute_DHCopy - ParallelItemShift] ,
	   "Redistribute atoms device-host copy",
	   MaxWordsLength);
  strncpy (hostWords[item_Redistribute_Barrier - ParallelItemShift] ,
	   "Redistribute barrier",
	   MaxWordsLength);
  strncpy (hostWords[item_TransferGhost - ParallelItemShift] ,
	   "Transfer Ghosts",
	   MaxWordsLength);
  strncpy (hostWords[item_TransferGhost_DHCopy - ParallelItemShift] ,
	   "Transfer Ghosts device-host copy",
	   MaxWordsLength);
  strncpy (hostWords[item_TransferGhost_DHCopy_Pack - ParallelItemShift] ,
	   "Transfer Ghosts device-host copy pack",
	   MaxWordsLength);
  strncpy (hostWords[item_TransferGhost_DHCopy_Copy - ParallelItemShift] ,
	   "Transfer Ghosts device-host copy copy",
	   MaxWordsLength);
  strncpy (hostWords[item_TransferGhost_DHCopy_Unpack - ParallelItemShift] ,
	   "Transfer Ghosts device-host copy unpack",
	   MaxWordsLength);
  strncpy (hostWords[item_TransferGhost_Tranfer - ParallelItemShift] ,
	   "Transfer Ghosts transfer",
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
