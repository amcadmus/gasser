#define DEVICE_CODE

#include "common.h"
#include "Ewald.h"
#include "SPMERec.h"
#include "CardinalBspline.h"

__global__ void
cal_Bx (const IntVectorType K,
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
cal_By (const IntVectorType K,
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
cal_Bz (const IntVectorType K,
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
	      cufftComplex * phiF2)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;

  IntVectorType N(K);
  N.z = (N.z >> 1) + 1;
  IndexType nele = K.x * K.y * K.z;
  IndexType nelehalf = N.x * N.y * N.z;
  ScalorType scalor0 =  1.f / (M_PIF * volume) * nele;
  ScalorType scalor1 = -2.f / volume * nele;
  
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
    IntVectorType meshIdx_base;
    meshIdx_base.x = int(uu.x);
    meshIdx_base.y = int(uu.y);
    meshIdx_base.z = int(uu.z);
    if      (meshIdx_base.x < 0   ) meshIdx_base.x += K.x;
    else if (meshIdx_base.x >= K.x) meshIdx_base.x -= K.x;
    if      (meshIdx_base.y < 0   ) meshIdx_base.y += K.y;
    else if (meshIdx_base.y >= K.y) meshIdx_base.y -= K.y;
    if      (meshIdx_base.z < 0   ) meshIdx_base.z += K.z;
    else if (meshIdx_base.z >= K.z) meshIdx_base.z -= K.z;
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
	    IndexType posi = index3to1 (meshIdx, K);
	    IndexType natomInList = atomicInc (&(nlist_n[posi]), (nlist_length + 1));
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

	

