#ifndef __Ewald_h_wanghan__
#define __Ewald_h_wanghan__

#define DEVICE_CODE

#include "common.h"

__device__ ScalorType static
kernel_rm1_rec_f (const ScalorType m2,
		  const ScalorType beta)
{
  ScalorType tmp = M_PIF / beta;
  return expf (- tmp*tmp*m2) / m2;
}

const double global_twoSqrtPIi (1.12837916709551);

__device__ ScalorType static
kernel_rm1_dir_forceScale (const ScalorType r,
			   const ScalorType beta)
{
  ScalorType ri (ScalorType(1.)/r);
  ScalorType betar = beta * r;
  return - (erfcf(betar) * ri +
	    beta * global_twoSqrtPIi *
	    exp(-betar*betar)) * ri * ri;
}

__device__ ScalorType static
kernel_rm1_dir_energy_forceScale (const ScalorType r,
				  const ScalorType beta,
				  ScalorType * fscale)
{
  double ri = double(1.f)/r;
  double betar = beta * r;
  double erfcBetaRRi = erfc (betar) * ri;
  * fscale= - (erfcBetaRRi +
	       beta * global_twoSqrtPIi *
	       exp (- betar * betar)) * ri * ri;
  return erfcBetaRRi;
}


__device__ IndexType static
index3to1 (const IndexType ix,
	   const IndexType iy,
	   const IndexType iz,
	   const IndexType nx,
	   const IndexType ny,
	   const IndexType nz )
{
  return iz + nz * (iy + ny * ix);
}

__device__ IndexType static
index3to1 (const IntVectorType i,
	   const IntVectorType N)
{
  return i.z + N.z * (i.y + N.y * i.x);
}


__device__ void static
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


__device__ void static
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
