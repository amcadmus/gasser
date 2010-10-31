#define DEVICE_CODE

#include "Displacement_interface.h"

void Displacement_max::
reinit (const MDSystem & sys,
	const IndexType & NTread)
{
  myBlockDim.y = 1;
  myBlockDim.z = 1;
  myBlockDim.x = NTread;

  IndexType nob;
  if (sys.ddata.numAtom % myBlockDim.x == 0){
    nob = sys.ddata.numAtom / myBlockDim.x;
  } else {
    nob = sys.ddata.numAtom / myBlockDim.x + 1;
  }
  atomGridDim = toGridDim (nob);

  clearDisplacement ();
  mallocDisplacemant (sys);

  max.reinit (nob, NThreadForSum);
  displacement_sbuffSize = myBlockDim.x * sizeof(IndexType);
}

Displacement_max::
Displacement_max (const MDSystem & sys,
		  const IndexType & NThread)
    : malloced (false)
{
  reinit (sys, NThread);
}

Displacement_max::
~Displacement_max ()
{
  clearDisplacement();
}	

void Displacement_max ::
mallocDisplacemant (const MDSystem & sys)
{
  cudaMalloc ((void **)& backupCoord,  sizeof(CoordType) *sys.ddata.numAtom);
  // reshuffle backup
  cudaMalloc ((void **)& bkbackupCoord,sizeof(CoordType) *sys.ddata.numAtom);
  cudaMalloc ((void **)& dresult,  sizeof(ScalorType));
  
  checkCUDAError ("Displacement_max::init allocations");
  malloced = true;
}

void Displacement_max::
clearDisplacement ()
{
  if ( malloced ){
    cudaFree (backupCoord);
    cudaFree (bkbackupCoord);
    cudaFree (dresult);
    malloced = false;
    checkCUDAError ("NeighborList::clearJudgeStuff");
  }
}

void Displacement_max::
recordCoord (const MDSystem & sys)
{
  cpyProperty <<<atomGridDim, myBlockDim>>> (
      backupCoord,
      sys.ddata.coord,
      sys.ddata.numAtom);
  checkCUDAError ("NeighborList::init backup coords");
}

static __global__ void
displacement_max_block (const IndexType numAtom,
			const RectangularBox box,
			const CoordType * coord,
			const CoordType * backupCoord,
			ScalorType * buff)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;

  extern __shared__ volatile ScalorType mydiff[];

  ScalorType dx=0.f, dy=0.f, dz=0.f;
  if (ii < numAtom){
    dx = coord[ii].x - backupCoord[ii].x;
    dy = coord[ii].y - backupCoord[ii].y;
    dz = coord[ii].z - backupCoord[ii].z;
    shortestImage (box, &dx, &dy, &dz);
  }
  mydiff[threadIdx.x] = sqrtf(dx*dx + dy*dy + dz*dz);

  maxVectorBlockBuffer_2 (mydiff);
  if (threadIdx.x == 0){
    buff[bid] = mydiff[0];
  }
}


ScalorType Displacement_max::
calMaxDisplacemant (const MDSystem & sys,
		    MDTimer *timer)
{
  if (timer != NULL) timer->tic(mdTimeJudgeRebuild);
  
  displacement_max_block 
      <<<atomGridDim, myBlockDim,
      displacement_sbuffSize>>> (
  	  sys.ddata.numAtom,
	  sys.box,
	  sys.ddata.coord,
	  backupCoord,
  	  max.getBuff());
  max.maxBuff (dresult, 0);
  
  cudaMemcpy (&hresult, dresult, sizeof(ScalorType), cudaMemcpyDeviceToHost);

  if (timer != NULL) timer->toc(mdTimeJudgeRebuild);
  return hresult;
}











void Displacement_mean::
reinit (const MDSystem & sys,
	const IndexType & NTread)
{
  myBlockDim.y = 1;
  myBlockDim.z = 1;
  myBlockDim.x = NTread;

  IndexType nob;
  if (sys.ddata.numAtom % myBlockDim.x == 0){
    nob = sys.ddata.numAtom / myBlockDim.x;
  } else {
    nob = sys.ddata.numAtom / myBlockDim.x + 1;
  }
  atomGridDim = toGridDim (nob);

  clearDisplacement ();
  mallocDisplacemant (sys);

  sum.reinit (nob, NThreadForSum);
  displacement_sbuffSize = myBlockDim.x * sizeof(IndexType);
}

Displacement_mean::
Displacement_mean (const MDSystem & sys,
		   const IndexType & NThread)
    : malloced (false)
{
  reinit (sys, NThread);
}

Displacement_mean::
~Displacement_mean ()
{
  clearDisplacement();
}	

void Displacement_mean ::
mallocDisplacemant (const MDSystem & sys)
{
  cudaMalloc ((void **)& backupCoord,  sizeof(CoordType) *sys.ddata.numAtom);
  // reshuffle backup
  cudaMalloc ((void **)& bkbackupCoord,sizeof(CoordType) *sys.ddata.numAtom);
  cudaMalloc ((void **)& dresult,  sizeof(ScalorType));
  
  checkCUDAError ("Displacement_mean::init allocations");
  malloced = true;
}

void Displacement_mean::
clearDisplacement ()
{
  if ( malloced ){
    cudaFree (backupCoord);
    cudaFree (bkbackupCoord);
    cudaFree (dresult);
    malloced = false;
    checkCUDAError ("NeighborList::clearJudgeStuff");
  }
}

void Displacement_mean::
recordCoord (const MDSystem & sys)
{
  cpyProperty <<<atomGridDim, myBlockDim>>> (
      backupCoord,
      sys.ddata.coord,
      sys.ddata.numAtom);
  checkCUDAError ("NeighborList::init backup coords");
}

static __global__ void
displacement_mean_block (const IndexType numAtom,
			 const RectangularBox box,
			 const CoordType * coord,
			 const CoordType * backupCoord,
			 ScalorType * buff)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;

  extern __shared__ volatile ScalorType mydiff[];

  ScalorType dx=0.f, dy=0.f, dz=0.f;
  if (ii < numAtom){
    dx = coord[ii].x - backupCoord[ii].x;
    dy = coord[ii].y - backupCoord[ii].y;
    dz = coord[ii].z - backupCoord[ii].z;
    shortestImage (box, &dx, &dy, &dz);
  }
  mydiff[threadIdx.x] = sqrtf(dx*dx + dy*dy + dz*dz);

  sumVectorBlockBuffer_2 (mydiff);
  if (threadIdx.x == 0){
    buff[bid] = mydiff[0];
  }
}


ScalorType Displacement_mean::
calMeanDisplacemant (const MDSystem & sys,
		     MDTimer *timer)
{
  if (timer != NULL) timer->tic(mdTimeJudgeRebuild);
  
  displacement_mean_block 
      <<<atomGridDim, myBlockDim,
      displacement_sbuffSize>>> (
  	  sys.ddata.numAtom,
	  sys.box,
	  sys.ddata.coord,
	  backupCoord,
  	  sum.getBuff());
  sum.sumBuff (dresult, 0);
  
  cudaMemcpy (&hresult, dresult, sizeof(ScalorType), cudaMemcpyDeviceToHost);

  if (timer != NULL) timer->toc(mdTimeJudgeRebuild);
  return hresult / sys.ddata.numAtom;
}




