#ifndef __MDTimer_interface_h_wanghan__
#define __MDTimer_interface_h_wanghan__

#include "MDTimer.h"
#include <stdio.h>

class MDTimer
{
  cudaEvent_t start[NumberOfMemberInMdTime];
  cudaEvent_t stop [NumberOfMemberInMdTime];
  TimeType timeRecord [NumberOfMemberInMdTime];
public:
  MDTimer();
  ~MDTimer();
  void tic (mdTimeItem_t timeItem,
	    cudaStream_t stream = 0);
  TimeType toc (mdTimeItem_t timeItem,
		cudaStream_t stream = 0);
  void printRecord (FILE * fstream);
};

inline void MDTimer::tic(mdTimeItem_t item, cudaStream_t stream)
{
  cudaEventRecord(start[item], stream);
}

inline TimeType MDTimer::toc(mdTimeItem_t item, cudaStream_t stream)
{
  float tmptime;
  
  cudaEventRecord(stop[item], stream);
  cudaEventSynchronize (stop[item]);
  cudaEventElapsedTime (&tmptime, start[item], stop[item]);

  timeRecord[item] += TimeType(tmptime);
  return TimeType(tmptime);
}

#endif
