#include "Ewald.h"
#include "EwaldSumRec.h"


// mesh grid and block
void __global__
cal_fm (const IntVectorType K,
	const MatrixType vecAStar,
	const ScalorType beta,
	const ScalorType volume,
	ScalorType * fm,
	const IndexType nele)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;

  IntVectorType N;
  N.x = ((K.x>>1) << 1) + 1;
  N.y = ((K.y>>1) << 1) + 1;
  N.z = ((K.z>>1) << 1) + 1;

  if (ii < nele){
    IntVectorType idx;
    index1to3 (ii, N, &idx);
    IntVectorType im;
    im.x = idx.x - (K.x >> 1);
    im.y = idx.y - (K.y >> 1);
    im.z = idx.z - (K.z >> 1);
    VectorType mm;
    mm.x = im.x * vecAStar.xx + im.y * vecAStar.yx + im.z * vecAStar.zx;
    mm.y = im.x * vecAStar.xy + im.y * vecAStar.yy + im.z * vecAStar.zy;
    mm.z = im.x * vecAStar.xz + im.y * vecAStar.yz + im.z * vecAStar.zz;
    ScalorType m2 = (mm.x*mm.x + mm.y*mm.y + mm.z*mm.z);
    if (m2 != 0.f){
      fm[ii] = kernel_rm1_rec_f (m2, beta) / (2.f * M_PIF * volume);
    }
    else {
      fm[ii] = 0.f;
    }
  }
}

// mesh grid and block
// sbuffsize: sizeof(CoordType) * blockDim.x
void __global__
cal_Sm (const IntVectorType K,
	const MatrixType vecAStar,
	const CoordType * coord,
	const ScalorType * charge,
	const IndexType natom,
	cufftComplex * Sm,
	const IndexType nele)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  
  extern __shared__ char my_sbuff[];
  CoordType * sbuff_coord = (CoordType *) my_sbuff;
  sbuff_coord[tid].x = sbuff_coord[tid].y = sbuff_coord[tid].z = ScalorType(0.f);
  __syncthreads();

  VectorType mm;
  if (ii < nele){
    IntVectorType N;
    N.x = ((K.x>>1) << 1) + 1;
    N.y = ((K.y>>1) << 1) + 1;
    N.z = ((K.z>>1) << 1) + 1;
    // N.x = ((K.x>>1) << 1) + 1;
    // N.y = (K.y/2) * 2 + 1;
    // N.z = (K.z/2) * 2 + 1;
    IntVectorType idx;
    index1to3 (ii, N, &idx);
    IntVectorType im;
    im.x = idx.x - (K.x >> 1);
    im.y = idx.y - (K.y >> 1);
    im.z = idx.z - (K.z >> 1);
    mm.x = im.x * vecAStar.xx + im.y * vecAStar.yx + im.z * vecAStar.zx;
    mm.y = im.x * vecAStar.xy + im.y * vecAStar.yy + im.z * vecAStar.zy;
    mm.z = im.x * vecAStar.xz + im.y * vecAStar.yz + im.z * vecAStar.zz;
  }

  cufftComplex sum;
  sum.x = sum.y = 0.f;
  for (IndexType start = 0; start < natom; start += blockDim.x) {
    __syncthreads ();
    if (start + tid < natom){
      sbuff_coord[tid] = coord[start + tid];
      sbuff_coord[tid].w = charge[start + tid];
    }
    __syncthreads ();
    for (IndexType jj = 0; jj < blockDim.x; ++jj){
      if (jj + start == natom) break;
      if (ii < nele){
	ScalorType mr = 2.f * M_PIF * (
	    mm.x * sbuff_coord[jj].x +
	    mm.y * sbuff_coord[jj].y +
	    mm.z * sbuff_coord[jj].z );
	sum.x += sbuff_coord[jj].w * cos(mr);
	sum.y += sbuff_coord[jj].w * sin(mr);
      }
    }
  }

  if (ii < nele){
    Sm[ii] = sum;
  }
}


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
	    const IndexType natom)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  
  extern __shared__ char my_sbuff[];
  ScalorType * sbuff_fm = (ScalorType *) my_sbuff;
  cufftComplex * sbuff_Sm = (cufftComplex *) &(sbuff_fm[((blockDim.x >> 1) << 1) + 2]);

  VectorType fsum;
  fsum.x = fsum.y = fsum.z = ScalorType(0.f);
  VectorType fsum_small;
  fsum_small.x = fsum_small.y = fsum_small.z = ScalorType(0.f);
  ScalorType my_charge (0.f);
  CoordType  my_coord;
  my_coord.x = my_coord.y = my_coord.z = 0.f;
  if (ii < natom){
    my_charge = charge[ii];
    my_coord  = coord[ii];
  }
  IntVectorType N;
  N.x = ((K.x>>1) << 1) + 1;
  N.y = ((K.y>>1) << 1) + 1;
  N.z = ((K.z>>1) << 1) + 1;
      
  for (IndexType start = 0; start < nele; start += blockDim.x){
    __syncthreads();
    if (start + tid < nele){
      sbuff_fm[tid] = fm[start + tid];
      sbuff_Sm[tid] = Sm[start + tid];
    }
    __syncthreads();
    for (IndexType jj = 0; jj < blockDim.x; ++jj){
      IndexType meshIdx = start + jj;
      if (meshIdx == nele) break;
      IntVectorType idx;
      index1to3 (meshIdx, N, &idx);
      IntVectorType im;
      im.x = idx.x - (K.x >> 1);
      im.y = idx.y - (K.y >> 1);
      im.z = idx.z - (K.z >> 1);
      VectorType mm;
      mm.x = im.x * vecAStar.xx + im.y * vecAStar.yx + im.z * vecAStar.zx;
      mm.y = im.x * vecAStar.xy + im.y * vecAStar.yy + im.z * vecAStar.zy;
      mm.z = im.x * vecAStar.xz + im.y * vecAStar.yz + im.z * vecAStar.zz;

      if (ii < natom){
	ScalorType mr = 2.f * M_PIF * (
	    mm.x * my_coord.x +
	    mm.y * my_coord.y +
	    mm.z * my_coord.z );
	ScalorType scalor =
	    my_charge *
	    4.f * M_PIF *
	    sbuff_fm[jj] *
	    (sbuff_Sm[jj].x * sinf(mr) - sbuff_Sm[jj].y * cosf(mr));
	// ScalorType tmp;
	// tmp = scalor * mm.x;
	// if (tmp < 1e-3f) fsum_small.x += tmp;
	// else fsum.x += tmp;
	// tmp = scalor * mm.y;
	// if (tmp < 1e-3f) fsum_small.y += tmp;
	// else fsum.y += tmp;
	// tmp = scalor * mm.z;
	// if (tmp < 1e-3f) fsum_small.z += tmp;
	// else fsum.z += tmp;	  
	fsum.x += scalor * mm.x;
	fsum.y += scalor * mm.y;
	fsum.z += scalor * mm.z;
	// if (ii == 0){
	//   printf ("%f\n", fsum.y);
	// }
      }
    }
  }
  
  if (ii < natom){
    // forcx[ii] = fsum.x + fsum_small.x;
    // forcy[ii] = fsum.y + fsum_small.y;
    // forcz[ii] = fsum.z + fsum_small.z;
    forcx[ii] += fsum.x;
    forcy[ii] += fsum.y;
    forcz[ii] += fsum.z;
    // if (ii == 0){
    //   printf ("^^^^^^^^^^^^^^ %f\n", forcy[ii]);
    // }
  }  
}

void __global__
cal_energyPressure (const IntVectorType K,
		    const MatrixType vecAStar,
		    const ScalorType beta,
		    const ScalorType volume,
		    const ScalorType * fm,
		    const cufftComplex * Sm,
		    ScalorType * buff_e,
		    ScalorType * buff_pxx,
		    ScalorType * buff_pyy,
		    ScalorType * buff_pzz,
		    const IndexType nele)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  
  IntVectorType N;
  N.x = ((K.x>>1) << 1) + 1;
  N.y = ((K.y>>1) << 1) + 1;
  N.z = ((K.z>>1) << 1) + 1;

  if (ii < nele){
    IntVectorType idx;
    index1to3 (ii, N, &idx);
    IntVectorType im;
    im.x = idx.x - (K.x >> 1);
    im.y = idx.y - (K.y >> 1);
    im.z = idx.z - (K.z >> 1);
    VectorType mm;
    mm.x = im.x * vecAStar.xx + im.y * vecAStar.yx + im.z * vecAStar.zx;
    mm.y = im.x * vecAStar.xy + im.y * vecAStar.yy + im.z * vecAStar.zy;
    mm.z = im.x * vecAStar.xz + im.y * vecAStar.yz + im.z * vecAStar.zz;
    ScalorType m2 = (mm.x*mm.x + mm.y*mm.y + mm.z*mm.z);

    if (m2 == 0){
      buff_e[ii] = buff_pxx[ii] = buff_pyy[ii] = buff_pzz[ii] = 0.f;
    }
    else {
      ScalorType energy = fm[ii] * (Sm[ii].x * Sm[ii].x + Sm[ii].y * Sm[ii].y);
      buff_e[ii] = energy;
      ScalorType tmp = M_PIF / beta;
      tmp = (ScalorType(1.) + tmp * tmp * m2) / m2 * 2.f;
      buff_pxx[ii] = energy * (ScalorType(1.f) - tmp * mm.x * mm.x);
      buff_pyy[ii] = energy * (ScalorType(1.f) - tmp * mm.y * mm.y);
      buff_pzz[ii] = energy * (ScalorType(1.f) - tmp * mm.z * mm.z);    
    }
  }
}

