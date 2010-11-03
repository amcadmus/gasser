#ifndef __TwinRangeCorrectionRecorder_h__
#define __TwinRangeCorrectionRecorder_h__

#include "MDSystem_interface.h"
#include "Statistic_interface.h"
#include "SumVector.h"

class TwinRangeCorrectionRecorder
{
  ScalorType		energyCorr;
  ScalorType		pressureCorr;
  bool			malloced;
  dim3			atomGridDim;
  dim3			myBlockDim;
public:
  MDStatistic		st;
  ScalorType *		forcx;
  ScalorType *		forcy;
  ScalorType *		forcz;
private:
  void clear ();
public:
  TwinRangeCorrectionRecorder ();
  TwinRangeCorrectionRecorder (const MDSystem & sys,
			       const IndexType & NThread);
  ~TwinRangeCorrectionRecorder ();
  void reinit (const MDSystem & sys,
	       const IndexType & NThread);
public:
  const ScalorType & energyCorrection   () const {return energyCorr;}
  const ScalorType & pressureCorrection () const {return pressureCorr;}
  ScalorType & energyCorrection   () {return energyCorr;}
  ScalorType & pressureCorrection () {return pressureCorr;}
public:
  void correct (MDSystem & sys,
		MDTimer * timer = NULL) const;
  void correct (MDSystem & sys,
		MDStatistic & st,
		MDTimer * timer = NULL) const;
}
    ;


#endif
