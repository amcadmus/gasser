#include "Parallel_Timer.h"
#include <stdlib.h>
#include "compile_error_mixcode.h"

using namespace Parallel::Timer;

TimeType
Parallel::Timer::hostRecordReal [SizeOfTimeRecordArray] = {0};
TimeType
Parallel::Timer::hostRecordUser [SizeOfTimeRecordArray] = {0};
Stopwatch HostTimer::watch [SizeOfTimeRecordArray];

void Parallel::Timer::
printDeviceItem (FILE * fp,
		 timeItem_t item)
{
  fprintf(fp, "%s%s%1.3es   (% 3.1f%)\n",
	  deviceWords[item],
	  deviceSpace[item],
	  deviceRecord[item],
	  deviceRecord[item] * 100 /
	  hostRecordUser[item_Total - ParallelItemShift]);
}

void Parallel::Timer::
printRecord (FILE * fp)
{
  if (hostRecordUser[item_Total - ParallelItemShift] == 0) return;

  printDeviceItem (fp, item_BuildCellList);
  printDeviceItem (fp, item_ApplyBondaryCondition);
  printDeviceItem (fp, item_NonBondedInterStatistic);
  printDeviceItem (fp, item_Integrate);
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

