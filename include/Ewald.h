#ifndef __Ewald_h_wanghan__
#define __Ewald_h_wanghan__

#define DEVICE_CODE

#include "common.h"

__device__ ScalorType
kernel_rm1_rec_f (const ScalorType m2,
		  const ScalorType beta)
{
  ScalorType tmp = M_PIF / beta;
  return expf (- tmp*tmp*m2) / m2;
}

__device__ IndexType 
index3to1 (const IndexType ix,
	   const IndexType iy,
	   const IndexType iz,
	   const IndexType nx,
	   const IndexType ny,
	   const IndexType nz )
{
  return iz + nz * (iy + ny * ix);
}

__device__ IndexType 
index3to1 (const IntVectorType i,
	   const IntVectorType N)
{
  return i.z + N.z * (i.y + N.y * i.x);
}


__device__ void
index1to3 (const IndexType input,
	   const IndexType nx,
	   const IndexType ny,
	   const IndexType nz,
	   IndexType * ix,
	   IndexType * iy,
	   IndexType * iz) 
{
  IndexType tmp = input;
  *iz = tmp % (nz);
  tmp = (tmp - *iz) / nz;
  *iy = tmp % (ny);
  *ix =  (tmp - *iy) / ny;
}


__device__ void
index1to3 (const IndexType input,
	   const IntVectorType N,
	   IntVectorType * i)
{
  IndexType tmp = input;
  i->z = tmp % (N.z);
  tmp = (tmp - i->z) / N.z;
  i->y = tmp % (N.y);
  i->x =  (tmp - i->y) / N.y;
}


#endif
