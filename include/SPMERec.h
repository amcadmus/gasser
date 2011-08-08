#ifndef __SPMERec_h_wanghan__
#define __SPMERec_h_wanghan__

#define DEVICE_CODE

#include "systemDefines.h"
#include "Statistic_interface.h"
#include "MDSystem_interface.h"
#include "SumVector.h"
#include "MDError_interface.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>

class SPMERecIk
{
  IndexType order;
  IntVectorType K;
  IndexType nele;
  ScalorType volume;
  ScalorType beta;
  MatrixType vecA;
  MatrixType vecAStar;
  dim3 meshGridDim;
  dim3 meshGridDim_half;
  dim3 meshBlockDim;
  dim3 atomGridDim;
  dim3 atomBlockDim;
  bool malloced;
private:
  ScalorType * vecbx;
  ScalorType * vecby;
  ScalorType * vecbz;
  void calB ();
  cufftReal * Q;
  cufftComplex * psiF;
  cufftComplex * phiF0;
  cufftComplex * phiF1;
  cufftComplex * phiF2;
  cufftComplex * QF;
  cufftComplex * QFxPsiF;
  cufftComplex * QFxPhiF0;
  cufftComplex * QFxPhiF1;
  cufftComplex * QFxPhiF2;
  cufftReal * QconvPsi;
  cufftReal * QconvPhi0;
  cufftReal * QconvPhi1;
  cufftReal * QconvPhi2;
private:
  MDError err;
  IndexType * nlist_n;
  IndexType * nlist_list;
  IndexType nlist_stride;
  IndexType nlist_length;
  void buildNeighborList (const MDSystem & sys);
private:
  void calV();
  void calAStar();
  void freeAll ();
public:
  SPMERecIk ();
  ~SPMERecIk ();
  void reinit (const MatrixType & vecA,
	       const IntVectorType & K,
	       const IndexType & order,
	       const ScalorType & beta,
	       const IndexType & natom,
	       const IndexType & meshNThread = 128,
	       const IndexType & atomNThread = 128);
}
    ;


__global__ void
cal_Bx (const IntVectorType K,
	const IndexType order,
	ScalorType * b);
__global__ void
cal_By (const IntVectorType K,
	const IndexType order,
	ScalorType * b);
__global__ void
cal_Bz (const IntVectorType K,
	const IndexType order,
	ScalorType * b);
__global__ void
cal_PsiFPhiF (const IntVectorType K,
	      const MatrixType vecAStar,
	      const ScalorType beta,
	      const ScalorType volume,
	      const ScalorType * bx,
	      const ScalorType * by,
	      const ScalorType * bz,
	      cufftComplex * psiF,
	      cufftComplex * phiF0,
	      cufftComplex * phiF1,
	      cufftComplex * phiF2);






#endif



