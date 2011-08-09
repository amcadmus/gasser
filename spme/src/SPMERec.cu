#define DEVICE_CODE

#include "common.h"
#include "Ewald.h"
#include "SPMERec.h"
#include "CardinalBspline.h"

__global__ void
calBx (const IntVectorType K,
       const IndexType order,
       ScalorType * b)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;

  if (ii < int(K.x)){
    cufftComplex fenzi, fenmu;
    ScalorType tmp = ScalorType(2.f * M_PIF) * ScalorType((order-1) * ii) / ScalorType (K.x);
    fenzi.x = cosf(tmp);
    fenzi.y = sinf(tmp);

    fenmu.x = fenmu.y = 0.f;
    for (IndexType k = 0; k < order-1; k ++){
      ScalorType scale = BSplineValue (order, double(k+1));
      tmp = 2 * M_PI * ScalorType(ii * k) / ScalorType (K.x);
      fenmu.x += scale * cosf(tmp);
      fenmu.y += scale * sinf(tmp);
    }
    cufftComplex btmp;
    ScalorType scale = 1./ (fenmu.x*fenmu.x + fenmu.y*fenmu.y);
    btmp.x = scale * (fenzi.x * fenmu.x + fenzi.y * fenmu.y);
    btmp.y = scale * (fenzi.y * fenmu.x - fenzi.x * fenmu.y);
    b[ii] = btmp.x*btmp.x + btmp.y*btmp.y;
  }
}

__global__ void
calBy (const IntVectorType K,
       const IndexType order,
       ScalorType * b)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;

  if (ii < int(K.y)){
    cufftComplex fenzi, fenmu;
    ScalorType tmp = ScalorType(2.f * M_PIF) * ScalorType((order-1) * ii) / ScalorType (K.y);
    fenzi.x = cosf(tmp);
    fenzi.y = sinf(tmp);

    fenmu.x = fenmu.y = 0.f;
    for (IndexType k = 0; k < order-1; k ++){
      ScalorType scale = BSplineValue (order, double(k+1));
      tmp = 2.f * M_PIF * ScalorType(ii * k) / ScalorType (K.y);
      fenmu.x += scale * cosf(tmp);
      fenmu.y += scale * sinf(tmp);
    }
    cufftComplex btmp;
    ScalorType scale = 1./ (fenmu.x*fenmu.x + fenmu.y*fenmu.y);
    btmp.x = scale * (fenzi.x * fenmu.x + fenzi.y * fenmu.y);
    btmp.y = scale * (fenzi.y * fenmu.x - fenzi.x * fenmu.y);
    b[ii] = btmp.x*btmp.x + btmp.y*btmp.y;
  }
}

__global__ void
calBz (const IntVectorType K,
       const IndexType order,
       ScalorType * b)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;

  if (ii < int(K.z)){
    cufftComplex fenzi, fenmu;
    ScalorType tmp = ScalorType(2.f * M_PIF) * ScalorType((order-1) * ii) / ScalorType (K.z);
    // fenzi.x = cosf(tmp);
    // fenzi.y = sinf(tmp);
    fenzi.x = cos(tmp);
    fenzi.y = sin(tmp);

    fenmu.x = fenmu.y = 0.f;
    for (IndexType k = 0; k < order-1; k ++){
      ScalorType scale = BSplineValue (order, double(k+1));
      tmp = ScalorType(2.f * M_PIF) * ScalorType(ii * k) / ScalorType (K.z);
      // fenmu.x += scale * cosf(tmp);
      // fenmu.y += scale * sinf(tmp);
      fenmu.x += scale * cos(tmp);
      fenmu.y += scale * sin(tmp);
    }
    cufftComplex btmp;
    ScalorType scale = ScalorType(1.)/ (fenmu.x*fenmu.x + fenmu.y*fenmu.y);
    btmp.x = scale * (fenzi.x * fenmu.x + fenzi.y * fenmu.y);
    btmp.y = scale * (fenzi.y * fenmu.x - fenzi.x * fenmu.y);
    b[ii] = btmp.x*btmp.x + btmp.y*btmp.y;
  }
}


__global__ void
calPsiFPhiF (const IntVectorType K,
	     const MatrixType vecAStar,
	     const ScalorType beta,
	     const ScalorType volume,
	     const ScalorType * bx,
	     const ScalorType * by,
	     const ScalorType * bz,
	     cufftComplex * psiF,
	     cufftComplex * phiF0,
	     cufftComplex * phiF1,
	     cufftComplex * phiF2)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;

  IntVectorType N(K);
  N.z = (N.z >> 1) + 1;
  IndexType nele = K.x * K.y * K.z;
  IndexType nelehalf = N.x * N.y * N.z;
  ScalorType scalor0 =  0.5f / (M_PIF * volume) * nele;
  ScalorType scalor1 = -2.0f / volume * nele;
  
  if (ii < nelehalf) {
    IntVectorType idx;
    index1to3 (ii, N, &idx);
    ScalorType valueB =  bx[idx.x] * by[idx.y] * bz[idx.z];
    if (idx.x > (K.x >> 1)) idx.x -= K.x;
    if (idx.y > (K.y >> 1)) idx.y -= K.y;
    // if (idx.z > (K.z >> 1)) idx.z -= K.z;
    VectorType mm;
    mm.x = idx.x * vecAStar.xx + idx.y * vecAStar.yx + idx.z * vecAStar.zx;
    mm.y = idx.x * vecAStar.xy + idx.y * vecAStar.yy + idx.z * vecAStar.zy;
    mm.z = idx.x * vecAStar.xz + idx.y * vecAStar.yz + idx.z * vecAStar.zz;
    ScalorType m2 = (mm.x*mm.x + mm.y*mm.y + mm.z*mm.z);
    if (m2 == 0){
      psiF[ii].x  = psiF[ii].y  = 0.f;
      phiF0[ii].x = phiF0[ii].y = 0.f;
      phiF1[ii].x = phiF1[ii].y = 0.f;
      phiF2[ii].x = phiF2[ii].y = 0.f;
    }
    else {
      ScalorType expp = kernel_rm1_rec_f (m2, beta) * valueB;
      psiF[ii].x = scalor0 * expp;
      psiF[ii].y = 0.f;
      phiF0[ii].x = 0.f;
      phiF1[ii].x = 0.f;
      phiF2[ii].x = 0.f;
      phiF0[ii].y = scalor1 * expp * mm.x;
      phiF1[ii].y = scalor1 * expp * mm.y;
      phiF2[ii].y = scalor1 * expp * mm.z;
    }
  }
}

// mesh gridDim and blockDim
__global__ void 
initMeshNeighborList (const IntVectorType K,
		      const MatrixType vecAStar,
		      IndexType * nlist_n,
		      IndexType * nlist_list,
		      const IndexType nlist_stride,
		      const IndexType nlist_length)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;

  if (ii < K.x * K.y * K.z){
    nlist_n[ii] = 0;
  }
}


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
		       mdError_t * ptr_de )
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;

  if (ii < natom){
    CoordType my_coord (coord[ii]);
    VectorType uu;
    uu.x = K.x * (vecAStar.xx * my_coord.x + vecAStar.xy * my_coord.y + vecAStar.xz *my_coord.z);
    uu.y = K.y * (vecAStar.yx * my_coord.x + vecAStar.yy * my_coord.y + vecAStar.yz *my_coord.z);
    uu.z = K.z * (vecAStar.zx * my_coord.x + vecAStar.zy * my_coord.y + vecAStar.zz *my_coord.z);
    if      (uu.x <  0  ) uu.x += K.x;
    else if (uu.x >= K.x) uu.x -= K.x;
    if      (uu.y <  0  ) uu.y += K.y;
    else if (uu.y >= K.y) uu.y -= K.y;
    if      (uu.z <  0  ) uu.z += K.z;
    else if (uu.z >= K.z) uu.z -= K.z;
    IntVectorType meshIdx_base;
    meshIdx_base.x = int(uu.x);
    meshIdx_base.y = int(uu.y);
    meshIdx_base.z = int(uu.z);
    // if      (meshIdx_base.x < 0   ) meshIdx_base.x += K.x;
    // else if (meshIdx_base.x >= K.x) meshIdx_base.x -= K.x;
    // if      (meshIdx_base.y < 0   ) meshIdx_base.y += K.y;
    // else if (meshIdx_base.y >= K.y) meshIdx_base.y -= K.y;
    // if      (meshIdx_base.z < 0   ) meshIdx_base.z += K.z;
    // else if (meshIdx_base.z >= K.z) meshIdx_base.z -= K.z;
    if (meshIdx_base.x >= K.x || meshIdx_base.x < 0 ||
	meshIdx_base.y >= K.y || meshIdx_base.y < 0 ||
	meshIdx_base.z >= K.z || meshIdx_base.z < 0 ){
      *ptr_de = mdErrorOverFlowMeshIdx;
    }
    else {
      IntVectorType vd;
      IntVectorType meshIdx;
      for (vd.x = 0; vd.x < int(order); ++vd.x){
	meshIdx.x = meshIdx_base.x - vd.x;
	if (meshIdx.x < 0) meshIdx.x += K.x;
	for (vd.y = 0; vd.y < int(order); ++vd.y){
	  meshIdx.y = meshIdx_base.y - vd.y;
	  if (meshIdx.y < 0) meshIdx.y += K.y;
	  for (vd.z = 0; vd.z < int(order); ++vd.z){
	    meshIdx.z = meshIdx_base.z - vd.z;
	    if (meshIdx.z < 0) meshIdx.z += K.z;
	    //
	    // possible improvement!!!!
	    //
	    IndexType posi = index3to1 (meshIdx, K);
	    IndexType natomInList = atomicInc (&(nlist_n[posi]), (nlist_length));
	    if (natomInList >= nlist_length){
	      *ptr_de = mdErrorShortMeshNeighborList;
	    }
	    else {
	      nlist_list[nlist_stride * natomInList + posi] = ii;
	    }
	  }
	}
      }
    }
  }
}


__global__ void
timeQFPsiF (const cufftComplex * QF,
	    const cufftComplex * PsiF,
	    cufftComplex * QFxPsiF,
	    const IndexType nelehalf,
	    const ScalorType sizei)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  
  if (ii < nelehalf){
    QFxPsiF[ii].x = (QF[ii].x * PsiF[ii].x - QF[ii].y * PsiF[ii].y) * sizei;
    QFxPsiF[ii].y = (QF[ii].x * PsiF[ii].y + QF[ii].y * PsiF[ii].x) * sizei;
  }
}

__global__ void
timeQFPhiF (const cufftComplex * QF,
	    const cufftComplex * PhiF0,
	    const cufftComplex * PhiF1,
	    const cufftComplex * PhiF2,
	    cufftComplex * QFxPhiF0,
	    cufftComplex * QFxPhiF1,
	    cufftComplex * QFxPhiF2,
	    const IndexType nelehalf,
	    const ScalorType sizei)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  
  if (ii < nelehalf){
    QFxPhiF0[ii].x = (QF[ii].x * PhiF0[ii].x - QF[ii].y * PhiF0[ii].y) * sizei;
    QFxPhiF0[ii].y = (QF[ii].x * PhiF0[ii].y + QF[ii].y * PhiF0[ii].x) * sizei;
    QFxPhiF1[ii].x = (QF[ii].x * PhiF1[ii].x - QF[ii].y * PhiF1[ii].y) * sizei;
    QFxPhiF1[ii].y = (QF[ii].x * PhiF1[ii].y + QF[ii].y * PhiF1[ii].x) * sizei;
    QFxPhiF2[ii].x = (QF[ii].x * PhiF2[ii].x - QF[ii].y * PhiF2[ii].y) * sizei;
    QFxPhiF2[ii].y = (QF[ii].x * PhiF2[ii].y + QF[ii].y * PhiF2[ii].x) * sizei;
  }
}

// mesh grid and block
__global__ void
calEnergy (const cufftReal * Q,
	   const cufftReal * QConvPsi,
	   ScalorType * buff_e,
	   const IndexType nele)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;

  if (ii < nele){
    buff_e[ii] = Q[ii] * QConvPsi[ii];
  }
}

// atom grid and block
__global__ void
calForce (const IntVectorType K,
	  const MatrixType vecAStar,
	  const IndexType order,
	  const CoordType * coord,
	  const ScalorType * charge,
	  const IndexType natom,
	  const cufftReal * QConvPhi0,
	  const cufftReal * QConvPhi1,
	  const cufftReal * QConvPhi2,
	  ScalorType * forcx,
	  ScalorType * forcy,
	  ScalorType * forcz,
	  mdError_t * ptr_de )
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;

  if (ii < natom){
    CoordType my_coord (coord[ii]);
    VectorType uu;
    uu.x = K.x * (vecAStar.xx * my_coord.x + vecAStar.xy * my_coord.y + vecAStar.xz *my_coord.z);
    uu.y = K.y * (vecAStar.yx * my_coord.x + vecAStar.yy * my_coord.y + vecAStar.yz *my_coord.z);
    uu.z = K.z * (vecAStar.zx * my_coord.x + vecAStar.zy * my_coord.y + vecAStar.zz *my_coord.z);
    if      (uu.x <  0  ) uu.x += K.x;
    else if (uu.x >= K.x) uu.x -= K.x;
    if      (uu.y <  0  ) uu.y += K.y;
    else if (uu.y >= K.y) uu.y -= K.y;
    if      (uu.z <  0  ) uu.z += K.z;
    else if (uu.z >= K.z) uu.z -= K.z;
    IntVectorType meshIdx;
    meshIdx.x = int(uu.x);
    meshIdx.y = int(uu.y);
    meshIdx.z = int(uu.z);
    if (meshIdx.x >= K.x || meshIdx.x < 0 ||
	meshIdx.y >= K.y || meshIdx.y < 0 ||
	meshIdx.z >= K.z || meshIdx.z < 0 ){
      *ptr_de = mdErrorOverFlowMeshIdx;
    }
    else {
      ScalorType Mnx[MaxSPMEOrder];
      ScalorType Mny[MaxSPMEOrder];
      ScalorType Mnz[MaxSPMEOrder];
      {
	VectorType value;
	value.x = uu.x - meshIdx.x;
	value.y = uu.y - meshIdx.y;
	value.z = uu.z - meshIdx.z;
	for (IndexType jj = 0; jj < order; ++jj){
	  Mnx[jj] = BSplineValue (order, value.x);
	  Mny[jj] = BSplineValue (order, value.y);
	  Mnz[jj] = BSplineValue (order, value.z);
	  value.x += 1.;
	  value.y += 1.;
	  value.z += 1.;
	}
      }
      ScalorType mycharge = charge[ii];
      VectorType fsum;
      fsum.x = fsum.y = fsum.z = ScalorType(0.f);
      IntVectorType myMeshIdx;
      IntVectorType dk;
      for (dk.x = 0; dk.x < order; ++dk.x){
	myMeshIdx.x = meshIdx.x - dk.x;
	if (myMeshIdx.x < 0) myMeshIdx.x += K.x;
	for (dk.y = 0; dk.y < order; ++dk.y){
	  myMeshIdx.y = meshIdx.y - dk.y;
	  if (myMeshIdx.y < 0) myMeshIdx.y += K.y;
	  for (dk.z = 0; dk.z < order; ++dk.z){
	    myMeshIdx.z = meshIdx.z - dk.z;
	    if (myMeshIdx.z < 0) myMeshIdx.z += K.z;
	    //
	    // possible improvement!!!!
	    //
	    IndexType index = index3to1 (myMeshIdx, K);
	    ScalorType myP = Mnx[dk.x] * Mny[dk.y] * Mnz[dk.z];
	    //
	    // possible improvement of reading mem!!!!
	    //
	    fsum.x += myP * QConvPhi0[index];
	    fsum.y += myP * QConvPhi1[index];
	    fsum.z += myP * QConvPhi2[index];
	  }
	}
      }
      forcx[ii] += mycharge * fsum.x;
      forcy[ii] += mycharge * fsum.y;
      forcz[ii] += mycharge * fsum.z;
    }
  }
}


