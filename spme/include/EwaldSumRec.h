#ifndef __EwaldSumRec_h_wanghan__
#define __EwaldSumRec_h_wanghan__

#define DEVICE_CODE

#include "systemDefines.h"
#include "Statistic_interface.h"
#include "MDSystem_interface.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>

class EwaldSumRec
{
  IntVectorType K;
  IntVectorType N;
  IndexType nele;
  ScalorType volume;
  ScalorType beta;
  MatrixType vecA;
  MatrixType vecAStar;
  cufftComplex * Sm;
  ScalorType * fm;
  dim3 meshGridDim;
  dim3 meshBlockDim;
  dim3 atomGridDim;
  dim3 atomBlockDim;
  bool malloced;
private:
  void calV();
  void calAStar();
  void freeAll ();
public:
  EwaldSumRec ();
  ~EwaldSumRec ();
  void reinit (const MatrixType & vecA,
	       const IntVectorType & K,
	       const ScalorType & beta,
	       const IndexType & natom,
	       const IndexType & meshNThread = 128,
	       const IndexType & atomNThread = 128);
public:
  void applyInteraction (MDSystem & sys,
			 MDStatistic * st = NULL,
			 MDTimer *timer = NULL);
};



// mesh grid and block
void __global__
cal_fm (const IntVectorType K,
	const MatrixType vecAStar,
	const ScalorType beta,
	const ScalorType volume,
	ScalorType * fm,
	const IndexType nele);
// mesh grid and block
// sbuffsize: sizeof(CoordType) * blockDim.x
void __global__
cal_Sm (const IntVectorType K,
	const MatrixType vecAStar,
	const CoordType * coord,
	const ScalorType * charge,
	const IndexType natom,
	cufftComplex * Sm,
	const IndexType nele);
// atom grid and atom block	    
void __global__     
applyForce (const IntVectorType K,
	    const MatrixType vecAStar,
	    const ScalorType * fm,
	    const cufftComplex * Sm,
	    const IndexType nele,
	    const CoordType * coord,
	    const ScalorType * charge,
	    ScalorType * forcx,
	    ScalorType * forcy,
	    ScalorType * forcz,
	    const IndexType natom);



#endif
