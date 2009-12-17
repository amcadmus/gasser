#ifndef __AngleInteraction_h_wanghan__
#define __AngleInteraction_h_wanghan__

#include "common.h"

enum mdAngleInteraction{
  mdForceAngleHarmonic		= 0
};

typedef int mdAngleInteraction_t;


namespace AngleHarmonic {
    typedef enum paramIndex{
      k		= 0,
      theta0	= 1
    } paramIndex_t;

    __host__ void initParameter (ScalorType * param,
				 ScalorType k,
				 ScalorType theta0);
    __device__ void force0 (const ScalorType * param,
			    const ScalorType diff0x,
			    const ScalorType diff0y,
			    const ScalorType diff0z,
			    const ScalorType diff1x,
			    const ScalorType diff1y,
			    const ScalorType diff1z,
			    ScalorType * f0x,
			    ScalorType * f0y,
			    ScalorType * f0z);
    __device__ ScalorType forcePoten0 (const ScalorType * param,
				       const ScalorType diff0x,
				       const ScalorType diff0y,
				       const ScalorType diff0z,
				       const ScalorType diff1x,
				       const ScalorType diff1y,
				       const ScalorType diff1z,
				       ScalorType * f0x,
				       ScalorType * f0y,
				       ScalorType * f0z);
    __device__ void force1 (const ScalorType * param,
			    const ScalorType diff0x,
			    const ScalorType diff0y,
			    const ScalorType diff0z,
			    const ScalorType diff1x,
			    const ScalorType diff1y,
			    const ScalorType diff1z,
			    ScalorType * f0x,
			    ScalorType * f0y,
			    ScalorType * f0z);
    __device__ ScalorType forcePoten1 (const ScalorType * param,
				       const ScalorType diff0x,
				       const ScalorType diff0y,
				       const ScalorType diff0z,
				       const ScalorType diff1x,
				       const ScalorType diff1y,
				       const ScalorType diff1z,
				       ScalorType * f0x,
				       ScalorType * f0y,
				       ScalorType * f0z);
    
};


__device__ void AngleHarmonic::force0 (const ScalorType * param,
				      const ScalorType diff0x,
				      const ScalorType diff0y,
				      const ScalorType diff0z,
				      const ScalorType diff1x,
				      const ScalorType diff1y,
				      const ScalorType diff1z,
				      ScalorType * f0x,
				      ScalorType * f0y,
				      ScalorType * f0z)
{
  ScalorType c00 = (diff0x * diff0x + diff0y * diff0y + diff0z * diff0z);
  ScalorType c01 = (diff0x * diff1x + diff0y * diff1y + diff0z * diff1z);
  ScalorType c11 = (diff1x * diff1x + diff1y * diff1y + diff1z * diff1z);

  ScalorType sqrtc00i = 1.f / sqrtf(c00);
  ScalorType sqrtc11i = 1.f / sqrtf(c11);

  ScalorType tmp = sqrtc00i * sqrtc11i;
  ScalorType costheta = c01 * tmp;
  ScalorType theta = acosf (costheta);

  ScalorType scalor = 2 * param[k] * (theta - param[theta0]);
  scalor /= sinf(theta);
  scalor *= -tmp;
  
  *f0x = scalor * (diff1x - c01 * sqrtc00i * sqrtc00i * diff0x);
  *f0y = scalor * (diff1y - c01 * sqrtc00i * sqrtc00i * diff0y);
  *f0z = scalor * (diff1z - c01 * sqrtc00i * sqrtc00i * diff0z);
}


__device__ ScalorType AngleHarmonic::forcePoten0 (const ScalorType * param,
						  const ScalorType diff0x,
						  const ScalorType diff0y,
						  const ScalorType diff0z,
						  const ScalorType diff1x,
						  const ScalorType diff1y,
						  const ScalorType diff1z,
						  ScalorType * f0x,
						  ScalorType * f0y,
						  ScalorType * f0z)
{
  ScalorType c00 = (diff0x * diff0x + diff0y * diff0y + diff0z * diff0z);
  ScalorType c01 = (diff0x * diff1x + diff0y * diff1y + diff0z * diff1z);
  ScalorType c11 = (diff1x * diff1x + diff1y * diff1y + diff1z * diff1z);

  ScalorType sqrtc00i = 1.f / sqrtf(c00);
  ScalorType sqrtc11i = 1.f / sqrtf(c11);

  ScalorType tmp = sqrtc00i * sqrtc11i;
  ScalorType costheta = c01 * tmp;
  ScalorType theta = acosf (costheta);

  ScalorType scalor = 2 * param[k] * (theta - param[theta0]);
  scalor /= sinf(theta);
  scalor *= -tmp;
  
  *f0x = scalor * (diff1x - c01 * sqrtc00i * sqrtc00i * diff0x);
  *f0y = scalor * (diff1y - c01 * sqrtc00i * sqrtc00i * diff0y);
  *f0z = scalor * (diff1z - c01 * sqrtc00i * sqrtc00i * diff0z);

  return param[k] * (theta - param[theta0]) * (theta - param[theta0]);
}


__device__ void AngleHarmonic::force1 (const ScalorType * param,
				      const ScalorType diff0x,
				      const ScalorType diff0y,
				      const ScalorType diff0z,
				      const ScalorType diff1x,
				      const ScalorType diff1y,
				      const ScalorType diff1z,
				      ScalorType * f1x,
				      ScalorType * f1y,
				      ScalorType * f1z)
{
  ScalorType c00 = (diff0x * diff0x + diff0y * diff0y + diff0z * diff0z);
  ScalorType c01 = (diff0x * diff1x + diff0y * diff1y + diff0z * diff1z);
  ScalorType c11 = (diff1x * diff1x + diff1y * diff1y + diff1z * diff1z);

  ScalorType sqrtc00i = 1.f / sqrtf(c00);
  ScalorType sqrtc11i = 1.f / sqrtf(c11);

  ScalorType tmp = sqrtc00i * sqrtc11i;
  ScalorType costheta = c01 * tmp;
  ScalorType theta = acosf (costheta);

  ScalorType scalor = 2 * param[k] * (theta - param[theta0]);
  scalor /= sinf(theta);
  scalor *= tmp;
  
  *f1x = scalor * (diff1x - c01 * sqrtc00i * sqrtc00i * diff0x - diff0x + c01 * sqrtc11i * sqrtc11i * diff1x);
  *f1y = scalor * (diff1y - c01 * sqrtc00i * sqrtc00i * diff0y - diff0y + c01 * sqrtc11i * sqrtc11i * diff1y);
  *f1z = scalor * (diff1z - c01 * sqrtc00i * sqrtc00i * diff0z - diff0z + c01 * sqrtc11i * sqrtc11i * diff1z);
}


__device__ ScalorType AngleHarmonic::forcePoten1 (const ScalorType * param,
						  const ScalorType diff0x,
						  const ScalorType diff0y,
						  const ScalorType diff0z,
						  const ScalorType diff1x,
						  const ScalorType diff1y,
						  const ScalorType diff1z,
						  ScalorType * f1x,
						  ScalorType * f1y,
						  ScalorType * f1z)
{
  ScalorType c00 = (diff0x * diff0x + diff0y * diff0y + diff0z * diff0z);
  ScalorType c01 = (diff0x * diff1x + diff0y * diff1y + diff0z * diff1z);
  ScalorType c11 = (diff1x * diff1x + diff1y * diff1y + diff1z * diff1z);

  ScalorType sqrtc00i = 1.f / sqrtf(c00);
  ScalorType sqrtc11i = 1.f / sqrtf(c11);

  ScalorType tmp = sqrtc00i * sqrtc11i;
  ScalorType costheta = c01 * tmp;
  ScalorType theta = acosf (costheta);

  ScalorType scalor = 2 * param[k] * (theta - param[theta0]);
  scalor /= sinf(theta);
  scalor *= tmp;
  
  *f1x = scalor * (diff1x - c01 * sqrtc00i * sqrtc00i * diff0x - diff0x + c01 * sqrtc11i * sqrtc11i * diff1x);
  *f1y = scalor * (diff1y - c01 * sqrtc00i * sqrtc00i * diff0y - diff0y + c01 * sqrtc11i * sqrtc11i * diff1y);
  *f1z = scalor * (diff1z - c01 * sqrtc00i * sqrtc00i * diff0z - diff0z + c01 * sqrtc11i * sqrtc11i * diff1z);

  return param[k] * (theta - param[theta0]) * (theta - param[theta0]);  
}



__device__ void angleForce0 (const mdAngleInteraction_t ftype,
			     ScalorType * param,
			     const ScalorType diff0x,
			     const ScalorType diff0y,
			     const ScalorType diff0z,
			     const ScalorType diff1x,
			     const ScalorType diff1y,
			     const ScalorType diff1z,
			     ScalorType * f0x,
			     ScalorType * f0y,
			     ScalorType * f0z)
{
  if (ftype == mdForceAngleHarmonic){
    AngleHarmonic::force0 (param,
			   diff0x, diff0y, diff0z, diff1x, diff1y, diff1z,
			   f0x, f0y, f0z);
  }
}

__device__ void angleForcePoten0 (const mdAngleInteraction_t ftype,
				  ScalorType * param,
				  const ScalorType diff0x,
				  const ScalorType diff0y,
				  const ScalorType diff0z,
				  const ScalorType diff1x,
				  const ScalorType diff1y,
				  const ScalorType diff1z,
				  ScalorType * f0x,
				  ScalorType * f0y,
				  ScalorType * f0z,
				  ScalorType * dp)
{
  if (ftype == mdForceAngleHarmonic){
    *dp = AngleHarmonic::forcePoten0 (param,
				      diff0x, diff0y, diff0z,
				      diff1x, diff1y, diff1z,
				      f0x, f0y, f0z);
  }
}


__device__ void angleForce1 (const mdAngleInteraction_t ftype,
			     ScalorType * param,
			     const ScalorType diff0x,
			     const ScalorType diff0y,
			     const ScalorType diff0z,
			     const ScalorType diff1x,
			     const ScalorType diff1y,
			     const ScalorType diff1z,
			     ScalorType * f1x,
			     ScalorType * f1y,
			     ScalorType * f1z)
{
  if (ftype == mdForceAngleHarmonic){
    AngleHarmonic::force1 (param,
			   diff0x, diff0y, diff0z, diff1x, diff1y, diff1z,
			   f1x, f1y, f1z);
  }
}

__device__ void angleForcePoten1 (const mdAngleInteraction_t ftype,
				  ScalorType * param,
				  const ScalorType diff0x,
				  const ScalorType diff0y,
				  const ScalorType diff0z,
				  const ScalorType diff1x,
				  const ScalorType diff1y,
				  const ScalorType diff1z,
				  ScalorType * f1x,
				  ScalorType * f1y,
				  ScalorType * f1z,
				  ScalorType * dp)
{
  if (ftype == mdForceAngleHarmonic){
    *dp = AngleHarmonic::forcePoten1 (param,
				      diff0x, diff0y, diff0z,
				      diff1x, diff1y, diff1z,
				      f1x, f1y, f1z);
  }
}






#endif
