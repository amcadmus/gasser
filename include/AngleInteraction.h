#ifndef __AngleInteraction_h_wanghan__
#define __AngleInteraction_h_wanghan__

#include "common.h"
#include "Interaction.h"

enum mdAngleInteractionNParam{
  mdForceNParamAngleHarmonic	= 2
};

class AngleInteractionParameter : public InteractionParamter
{
public:
  bool same (const AngleInteractionParameter & f1) const ;
  bool operator == (const AngleInteractionParameter & f1) const;
};

class AngleHarmonicParameter : public AngleInteractionParameter
{
  ScalorType param [mdForceNParamAngleHarmonic];
public:
  virtual InteractionType type () const;
  virtual unsigned numParam () const ;
  virtual const ScalorType * c_ptr () const ;
  void init (ScalorType k,
	     ScalorType theta0);
};

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



__device__ void angleForce (const bool center,
			    const InteractionType ftype,
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
    if (center){
      AngleHarmonic::force1 (param,
			     diff0x, diff0y, diff0z, diff1x, diff1y, diff1z,
			     f0x, f0y, f0z);
    } else {
      AngleHarmonic::force0 (param,
			     diff0x, diff0y, diff0z, diff1x, diff1y, diff1z,
			     f0x, f0y, f0z);
    }
  }
}

__device__ void angleForcePoten (const bool center,
				 const InteractionType ftype,
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
    if (center){
      *dp = AngleHarmonic::forcePoten1 (param,
					diff0x, diff0y, diff0z,
					diff1x, diff1y, diff1z,
					f0x, f0y, f0z);
    } else {
      *dp = AngleHarmonic::forcePoten0 (param,
					diff0x, diff0y, diff0z,
					diff1x, diff1y, diff1z,
					f0x, f0y, f0z);
    }
  }
}


inline __host__ void AngleHarmonic::initParameter (ScalorType * param,
						   ScalorType k_,
						   ScalorType theta0_)
{
  param[k] = k_;
  param[theta0] = M_PIF - theta0_;
}

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

  ScalorType ex = (diff1x - c01 * sqrtc00i * sqrtc00i * diff0x);
  ScalorType ey = (diff1y - c01 * sqrtc00i * sqrtc00i * diff0y);
  ScalorType ez = (diff1z - c01 * sqrtc00i * sqrtc00i * diff0z);

  ScalorType theta = acosf (c01 * sqrtc11i * sqrtc00i);

  ScalorType scalor = 2 * param[k] * (theta - param[theta0]);
  scalor *= sqrtc00i;
  scalor *= -1. / sqrtf (ex*ex + ey*ey + ez*ez);

  *f0x = ex * scalor;
  *f0y = ey * scalor;
  *f0z = ez * scalor;
  
  // ScalorType tmp = sqrtc00i * sqrtc11i;
  // ScalorType costheta = c01 * tmp;
  // ScalorType theta = acosf (costheta);

  // ScalorType scalor = 2 * param[k] * (theta - param[theta0]);
  // scalor /= sinf(theta);
  // scalor *= -tmp;
  
  // *f0x = scalor * (diff1x - c01 * sqrtc00i * sqrtc00i * diff0x);
  // *f0y = scalor * (diff1y - c01 * sqrtc00i * sqrtc00i * diff0y);
  // *f0z = scalor * (diff1z - c01 * sqrtc00i * sqrtc00i * diff0z);
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

  ScalorType ex = (diff1x - c01 * sqrtc00i * sqrtc00i * diff0x);
  ScalorType ey = (diff1y - c01 * sqrtc00i * sqrtc00i * diff0y);
  ScalorType ez = (diff1z - c01 * sqrtc00i * sqrtc00i * diff0z);

  ScalorType theta = acosf (c01 * sqrtc11i * sqrtc00i);

  ScalorType scalor = 2 * param[k] * (theta - param[theta0]);
  scalor *= sqrtc00i;
  scalor *= -1. / sqrtf (ex*ex + ey*ey + ez*ez);

  *f0x = ex * scalor;
  *f0y = ey * scalor;
  *f0z = ez * scalor;

  // ScalorType c00 = (diff0x * diff0x + diff0y * diff0y + diff0z * diff0z);
  // ScalorType c01 = (diff0x * diff1x + diff0y * diff1y + diff0z * diff1z);
  // ScalorType c11 = (diff1x * diff1x + diff1y * diff1y + diff1z * diff1z);

  // ScalorType sqrtc00i = 1.f / sqrtf(c00);
  // ScalorType sqrtc11i = 1.f / sqrtf(c11);

  // ScalorType tmp = sqrtc00i * sqrtc11i;
  // ScalorType costheta = c01 * tmp;
  // ScalorType theta = acosf (costheta);

  // ScalorType scalor = 2 * param[k] * (theta - param[theta0]);
  // scalor /= sinf(theta);
  // scalor *= -tmp;
  
  // *f0x = scalor * (diff1x - c01 * sqrtc00i * sqrtc00i * diff0x);
  // *f0y = scalor * (diff1y - c01 * sqrtc00i * sqrtc00i * diff0y);
  // *f0z = scalor * (diff1z - c01 * sqrtc00i * sqrtc00i * diff0z);

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

  ScalorType ex = (diff1x - c01 * sqrtc00i * sqrtc00i * diff0x);
  ScalorType ey = (diff1y - c01 * sqrtc00i * sqrtc00i * diff0y);
  ScalorType ez = (diff1z - c01 * sqrtc00i * sqrtc00i * diff0z);

  ScalorType gx = (diff0x - c01 * sqrtc11i * sqrtc11i * diff1x);
  ScalorType gy = (diff0y - c01 * sqrtc11i * sqrtc11i * diff1y);
  ScalorType gz = (diff0z - c01 * sqrtc11i * sqrtc11i * diff1z);
  
  ScalorType theta = acosf (c01 * sqrtc11i * sqrtc00i);

  ScalorType scalor = 2 * param[k] * (theta - param[theta0]);
  ScalorType scalor0 = (1.f / sqrtf (ex*ex + ey*ey + ez*ez)) * sqrtc00i * scalor;
  ScalorType scalor1 = (1.f / sqrtf (gx*gx + gy*gy + gz*gz)) * sqrtc11i * scalor;

  *f1x = scalor0 * ex - scalor1 * gx;
  *f1y = scalor0 * ey - scalor1 * gy;
  *f1z = scalor0 * ez - scalor1 * gz;

  // ScalorType c00 = (diff0x * diff0x + diff0y * diff0y + diff0z * diff0z);
  // ScalorType c01 = (diff0x * diff1x + diff0y * diff1y + diff0z * diff1z);
  // ScalorType c11 = (diff1x * diff1x + diff1y * diff1y + diff1z * diff1z);

  // ScalorType sqrtc00i = 1.f / sqrtf(c00);
  // ScalorType sqrtc11i = 1.f / sqrtf(c11);

  // ScalorType tmp = sqrtc00i * sqrtc11i;
  // ScalorType costheta = c01 * tmp;
  // ScalorType theta = acosf (costheta);

  // ScalorType scalor = 2 * param[k] * (theta - param[theta0]);
  // scalor /= sinf(theta);
  // scalor *= tmp;
  
  // *f1x = scalor * (diff1x - c01 * sqrtc00i * sqrtc00i * diff0x - diff0x + c01 * sqrtc11i * sqrtc11i * diff1x);
  // *f1y = scalor * (diff1y - c01 * sqrtc00i * sqrtc00i * diff0y - diff0y + c01 * sqrtc11i * sqrtc11i * diff1y);
  // *f1z = scalor * (diff1z - c01 * sqrtc00i * sqrtc00i * diff0z - diff0z + c01 * sqrtc11i * sqrtc11i * diff1z);
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

  ScalorType ex = (diff1x - c01 * sqrtc00i * sqrtc00i * diff0x);
  ScalorType ey = (diff1y - c01 * sqrtc00i * sqrtc00i * diff0y);
  ScalorType ez = (diff1z - c01 * sqrtc00i * sqrtc00i * diff0z);

  ScalorType gx = (diff0x - c01 * sqrtc11i * sqrtc11i * diff1x);
  ScalorType gy = (diff0y - c01 * sqrtc11i * sqrtc11i * diff1y);
  ScalorType gz = (diff0z - c01 * sqrtc11i * sqrtc11i * diff1z);
  
  ScalorType theta = acosf (c01 * sqrtc11i * sqrtc00i);

  ScalorType scalor = 2 * param[k] * (theta - param[theta0]);
  ScalorType scalor0 = (1.f / sqrtf (ex*ex + ey*ey + ez*ez)) * sqrtc00i * scalor;
  ScalorType scalor1 = (1.f / sqrtf (gx*gx + gy*gy + gz*gz)) * sqrtc11i * scalor;

  *f1x = scalor0 * ex - scalor1 * gx;
  *f1y = scalor0 * ey - scalor1 * gy;
  *f1z = scalor0 * ez - scalor1 * gz;

  // ScalorType c00 = (diff0x * diff0x + diff0y * diff0y + diff0z * diff0z);
  // ScalorType c01 = (diff0x * diff1x + diff0y * diff1y + diff0z * diff1z);
  // ScalorType c11 = (diff1x * diff1x + diff1y * diff1y + diff1z * diff1z);

  // ScalorType sqrtc00i = 1.f / sqrtf(c00);
  // ScalorType sqrtc11i = 1.f / sqrtf(c11);

  // ScalorType tmp = sqrtc00i * sqrtc11i;
  // ScalorType costheta = c01 * tmp;
  // ScalorType theta = acosf (costheta);

  // ScalorType scalor = 2 * param[k] * (theta - param[theta0]);
  // scalor /= sinf(theta);
  // scalor *= tmp;
  
  // *f1x = scalor * (diff1x - c01 * sqrtc00i * sqrtc00i * diff0x - diff0x + c01 * sqrtc11i * sqrtc11i * diff1x);
  // *f1y = scalor * (diff1y - c01 * sqrtc00i * sqrtc00i * diff0y - diff0y + c01 * sqrtc11i * sqrtc11i * diff1y);
  // *f1z = scalor * (diff1z - c01 * sqrtc00i * sqrtc00i * diff0z - diff0z + c01 * sqrtc11i * sqrtc11i * diff1z);

  return param[k] * (theta - param[theta0]) * (theta - param[theta0]);  
}



__device__ void angleForce0 (const InteractionType ftype,
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

__device__ void angleForcePoten0 (const InteractionType ftype,
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


__device__ void angleForce1 (const InteractionType ftype,
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

__device__ void angleForcePoten1 (const InteractionType ftype,
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
