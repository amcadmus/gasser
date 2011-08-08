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

//
// consider NPT ensemble!!!!
//

class SPMERecIk
{
  IndexType order;
  IntVectorType K;
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
  cufftHandle planForward;
  cufftHandle planBackward;
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
  void calQ (const MDSystem & sys);
  void applyInteraction (MDSystem & sys,
			 MDStatistic * pst,
			 MDTimer * timer);
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


// mesh gridDim and blockDim
__global__ void 
initMeshNeighborList (const IntVectorType K,
		      const MatrixType vecAStar,
		      IndexType * nlist_n,
		      IndexType * nlist_list,
		      const IndexType nlist_stride,
		      const IndexType nlist_length);
// atom grid and block      
__global__ void
buildMeshNeighborList (const IntVectorType K,
		       const MatrixType vecAStar,
		       const IndexType order,
		       const CoordType * coord,
		       const IndexType natom,
		       IndexType * nlist_n,
		       IndexType * nlist_list,
		       const IndexType nlist_stride,
		       const IndexType nlist_length,
		       mdError_t * ptr_de );
// mesh gridDim and blockDim
__global__ void
cal_Q (const IntVectorType K,
       const MatrixType vecAStar,
       const IndexType order,
       const IndexType * nlist_n,
       const IndexType * nlist_list,
       const IndexType nlist_stride,
       cufftReal * Q);
// half mesh grid and block
__global__ void
timeQFPsiF (const cufftComplex * QF,
	    const cufftComplex * PsiF,
	    cufftComplex * QFxPsiF,
	    const IndexType nelehalf,
	    const ScalorType sizei);
// half mesh grid and block
__global__ void
timeQFPhiF (const cufftComplex * QF,
	    const cufftComplex * PhiF0,
	    const cufftComplex * PhiF1,
	    const cufftComplex * PhiF2,
	    cufftComplex * QFxPhiF0,
	    cufftComplex * QFxPhiF1,
	    cufftComplex * QFxPhiF2,
	    const IndexType nelehalf,
	    const ScalorType sizei);
// atom grid and block
__global__ void
calForce (const IntVectorType K,
	  const MatrixType vecAStar,
	  const IndexType order,
	  const CoordType * coord,
	  const ScalorType * charge,
	  const IndexType natom,
	  const cufftReal * QconvPhi0,
	  const cufftReal * QconvPhi1,
	  const cufftReal * QconvPhi2,
	  ScalorType * forcx,
	  ScalorType * forcy,
	  ScalorType * forcz,
	  mdError_t * ptr_de );


#endif



