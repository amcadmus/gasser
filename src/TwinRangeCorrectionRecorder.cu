#define DEVICE_CODE

#include "TwinRangeCorrectionRecorder.h"
#include "Auxiliary.h"

void TwinRangeCorrectionRecorder::
clear ()
{
  if (malloced) {
    cudaFree (forcx);
    cudaFree (forcy);
    cudaFree (forcz);
    checkCUDAError ("TwinRangeCorrectionRecorder::clear");
    malloced = false;
  }
}

void TwinRangeCorrectionRecorder::
reinit (const MDSystem & sys,
	const IndexType & NThread)
{
  clear();

  myBlockDim.x = NThread;
  IndexType nob = (sys.ddata.numAtom + myBlockDim.x - 1) / (myBlockDim.x);
  atomGridDim = toGridDim (nob);

  st.reinit (sys);
  
  energyCorr = ScalorType(0);
  pressureCorr = ScalorType(0);
  cudaMalloc ((void**)&forcx, sizeof(ScalorType) * sys.ddata.numAtom);
  cudaMalloc ((void**)&forcy, sizeof(ScalorType) * sys.ddata.numAtom);
  cudaMalloc ((void**)&forcz, sizeof(ScalorType) * sys.ddata.numAtom);
  checkCUDAError ("TwinRangeCorrectionRecorder::reinit");
  malloced = true;
}

TwinRangeCorrectionRecorder::
TwinRangeCorrectionRecorder ()
    : energyCorr(ScalorType(0)),
      pressureCorr(ScalorType(0)),
      malloced (false)
{
}

TwinRangeCorrectionRecorder::
~TwinRangeCorrectionRecorder ()
{
  clear();
}

TwinRangeCorrectionRecorder::
TwinRangeCorrectionRecorder (const MDSystem & sys,
			     const IndexType & NThread)
    : energyCorr(ScalorType(0)),
      pressureCorr(ScalorType(0)),
      malloced (false)
{
  reinit (sys, NThread);
}

// both blockDim and gridDim are 1
static __global__ void
applyTwinRangeEnergyPressureCorrection (ScalorType * ddata,
					ScalorType energyCorr,
					ScalorType pressureCorr)
{
  ddata[mdStatisticTwinRangeEnergyCorrection] = energyCorr;
  ddata[mdStatisticTwinRangePressureCorrection] = pressureCorr;
}


void TwinRangeCorrectionRecorder::
correct (MDSystem & sys,
	 MDTimer * timer) const
{
  if (timer != NULL) timer->tic(mdTimeNonBondedInteraction);
  addProperty
      <<<atomGridDim, myBlockDim>>> (
	  sys.ddata.forcx, forcx, sys.ddata.numAtom);
  addProperty
      <<<atomGridDim, myBlockDim>>> (
	  sys.ddata.forcy, forcy, sys.ddata.numAtom);
  addProperty
      <<<atomGridDim, myBlockDim>>> (
	  sys.ddata.forcz, forcz, sys.ddata.numAtom);
  checkCUDAError ("TwinRangeCorrectionRecorder::correct");
  if (timer != NULL) timer->toc(mdTimeNonBondedInteraction);
}

void TwinRangeCorrectionRecorder::
correct (MDSystem & sys,
	 MDStatistic & st,
	 MDTimer * timer) const
{
  if (timer != NULL) timer->tic(mdTimeNBInterStatistic);
  correct (sys);
  applyTwinRangeEnergyPressureCorrection
      <<<1, 1>>> (
	  st.ddata,
	  energyCorr,
	  pressureCorr);
  checkCUDAError ("TwinRangeCorrectionRecorder::correct");
  if (timer != NULL) timer->toc(mdTimeNBInterStatistic);
}



