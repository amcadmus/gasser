#ifndef __AngleInteraction_h_wanghan__
#define __AngleInteraction_h_wanghan__

#include "common.h"
#include "Interaction.h"

enum mdAngleInteractionNParam{
  mdForceNParamAngleHarmonic	= 2,
  mdForceNParamCosAngle0	= 1
};

class AngleHarmonicParameter : public AngleInteractionParameter
{
  ScalorType param [mdForceNParamAngleHarmonic];
public:
  AngleHarmonicParameter () {}
  AngleHarmonicParameter (ScalorType k,
			  ScalorType theta0);
  void reinit (ScalorType k,
	       ScalorType theta0);
  virtual InteractionType type () const;
  virtual unsigned numParam () const ;
  virtual ScalorType * c_ptr () ;
  virtual const ScalorType * c_ptr () const ;
};

class CosAngle0Parameter : public AngleInteractionParameter
{
  ScalorType param [mdForceNParamCosAngle0];
public:
  CosAngle0Parameter () {}
  CosAngle0Parameter (ScalorType k);
  void reinit (ScalorType k);
  virtual InteractionType type () const;
  virtual unsigned numParam () const ;
  virtual ScalorType * c_ptr () ;
  virtual const ScalorType * c_ptr () const ;
};


namespace AngleHarmonic {
  typedef enum paramIndex{
    k		= 0,
    theta0	= 1
  } paramIndex_t;

  __host__ void initParameter (ScalorType * param,
			       ScalorType k,
			       ScalorType theta0);
  static __device__ void
  force0 (const ScalorType * param,
	  const ScalorType diff0x,
	  const ScalorType diff0y,
	  const ScalorType diff0z,
	  const ScalorType diff1x,
	  const ScalorType diff1y,
	  const ScalorType diff1z,
	  ScalorType * f0x,
	  ScalorType * f0y,
	  ScalorType * f0z);
  static __device__ ScalorType
  forcePoten0 (const ScalorType * param,
	       const ScalorType diff0x,
	       const ScalorType diff0y,
	       const ScalorType diff0z,
	       const ScalorType diff1x,
	       const ScalorType diff1y,
	       const ScalorType diff1z,
	       ScalorType * f0x,
	       ScalorType * f0y,
	       ScalorType * f0z);
  static __device__ void
  force1 (const ScalorType * param,
	  const ScalorType diff0x,
	  const ScalorType diff0y,
	  const ScalorType diff0z,
	  const ScalorType diff1x,
	  const ScalorType diff1y,
	  const ScalorType diff1z,
	  ScalorType * f0x,
	  ScalorType * f0y,
	  ScalorType * f0z);
  static __device__ ScalorType
  forcePoten1 (const ScalorType * param,
	       const ScalorType diff0x,
	       const ScalorType diff0y,
	       const ScalorType diff0z,
	       const ScalorType diff1x,
	       const ScalorType diff1y,
	       const ScalorType diff1z,
	       ScalorType * f0x,
	       ScalorType * f0y,
	       ScalorType * f0z);
  static __device__ void
  force01 (const ScalorType * param,
	   const ScalorType diff0x,
	   const ScalorType diff0y,
	   const ScalorType diff0z,
	   const ScalorType diff1x,
	   const ScalorType diff1y,
	   const ScalorType diff1z,
	   ScalorType * f0x,
	   ScalorType * f0y,
	   ScalorType * f0z,
	   ScalorType * f1x,
	   ScalorType * f1y,
	   ScalorType * f1z);    
  static __device__ ScalorType
  forcePoten01 (const ScalorType * param,
		const ScalorType diff0x,
		const ScalorType diff0y,
		const ScalorType diff0z,
		const ScalorType diff1x,
		const ScalorType diff1y,
		const ScalorType diff1z,
		ScalorType * f0x,
		ScalorType * f0y,
		ScalorType * f0z,
		ScalorType * f1x,
		ScalorType * f1y,
		ScalorType * f1z);
};


namespace CosAngle0 {
  typedef enum paramIndex{
    k		= 0
  } paramIndex_t;

  __host__ void initParameter (ScalorType * param,
			       ScalorType k);
  static __device__ void
  force0 (const ScalorType * param,
	  const ScalorType diff0x,
	  const ScalorType diff0y,
	  const ScalorType diff0z,
	  const ScalorType diff1x,
	  const ScalorType diff1y,
	  const ScalorType diff1z,
	  ScalorType * f0x,
	  ScalorType * f0y,
	  ScalorType * f0z);
  static __device__ ScalorType
  forcePoten0 (const ScalorType * param,
	       const ScalorType diff0x,
	       const ScalorType diff0y,
	       const ScalorType diff0z,
	       const ScalorType diff1x,
	       const ScalorType diff1y,
	       const ScalorType diff1z,
	       ScalorType * f0x,
	       ScalorType * f0y,
	       ScalorType * f0z);
  static __device__ void
  force1 (const ScalorType * param,
	  const ScalorType diff0x,
	  const ScalorType diff0y,
	  const ScalorType diff0z,
	  const ScalorType diff1x,
	  const ScalorType diff1y,
	  const ScalorType diff1z,
	  ScalorType * f0x,
	  ScalorType * f0y,
	  ScalorType * f0z);
  static __device__ ScalorType
  forcePoten1 (const ScalorType * param,
	       const ScalorType diff0x,
	       const ScalorType diff0y,
	       const ScalorType diff0z,
	       const ScalorType diff1x,
	       const ScalorType diff1y,
	       const ScalorType diff1z,
	       ScalorType * f0x,
	       ScalorType * f0y,
	       ScalorType * f0z);
    
  static __device__ void
  force01 (const ScalorType * param,
	   const ScalorType diff0x,
	   const ScalorType diff0y,
	   const ScalorType diff0z,
	   const ScalorType diff1x,
	   const ScalorType diff1y,
	   const ScalorType diff1z,
	   ScalorType * f0x,
	   ScalorType * f0y,
	   ScalorType * f0z,
	   ScalorType * f1x,
	   ScalorType * f1y,
	   ScalorType * f1z);    
  static __device__ ScalorType
  forcePoten01 (const ScalorType * param,
		const ScalorType diff0x,
		const ScalorType diff0y,
		const ScalorType diff0z,
		const ScalorType diff1x,
		const ScalorType diff1y,
		const ScalorType diff1z,
		ScalorType * f0x,
		ScalorType * f0y,
		ScalorType * f0z,
		ScalorType * f1x,
		ScalorType * f1y,
		ScalorType * f1z);
};



static __device__ void
angleForce (const bool center,
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
  if (ftype == mdForceCosAngle0){
    if (center){
      CosAngle0::force1 (param,
			 diff0x, diff0y, diff0z, diff1x, diff1y, diff1z,
			 f0x, f0y, f0z);
    } else {
      CosAngle0::force0 (param,
			 diff0x, diff0y, diff0z, diff1x, diff1y, diff1z,
			 f0x, f0y, f0z);
    }
  }
}

static __device__ void
angleForcePoten (const bool center,
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
  if (ftype == mdForceCosAngle0){
    if (center){
      *dp = CosAngle0::forcePoten1 (param,
				    diff0x, diff0y, diff0z,
				    diff1x, diff1y, diff1z,
				    f0x, f0y, f0z);
    } else {
      *dp = CosAngle0::forcePoten0 (param,
				    diff0x, diff0y, diff0z,
				    diff1x, diff1y, diff1z,
				    f0x, f0y, f0z);
    }
  }
}


static __device__ void
angleForce (const bool center,
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
	    ScalorType * f1x,
	    ScalorType * f1y,
	    ScalorType * f1z)
{
  if (ftype == mdForceAngleHarmonic){
    AngleHarmonic::force01 (param,
			    diff0x, diff0y, diff0z, diff1x, diff1y, diff1z,
			    f0x, f0y, f0z, f1x, f1y, f1z);
  }
  if (ftype == mdForceCosAngle0){
    CosAngle0::force01 (param,
			diff0x, diff0y, diff0z, diff1x, diff1y, diff1z,
			f0x, f0y, f0z, f1x, f1y, f1z);
  }
}

static __device__ void
angleForcePoten (const bool center,
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
		 ScalorType * f1x,
		 ScalorType * f1y,
		 ScalorType * f1z,
		 ScalorType * dp)
{
  if (ftype == mdForceAngleHarmonic){
    *dp = AngleHarmonic::forcePoten01 (param,
				       diff0x, diff0y, diff0z,
				       diff1x, diff1y, diff1z,
				       f0x, f0y, f0z,
				       f1x, f1y, f1z);
  }
  if (ftype == mdForceCosAngle0){
    *dp = CosAngle0::forcePoten01 (param,
				   diff0x, diff0y, diff0z,
				   diff1x, diff1y, diff1z,
				   f0x, f0y, f0z,
				   f1x, f1y, f1z);
  }
}




inline __host__ void AngleHarmonic::initParameter (ScalorType * param,
						   ScalorType k_,
						   ScalorType theta0_)
{
  param[k] = k_;
  param[theta0] = M_PIF - theta0_;
}

__device__ void AngleHarmonic::
force0 (const ScalorType * param,
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

  ScalorType theta;
  ScalorType thetatmp = c01 * sqrtc11i * sqrtc00i;
  if (thetatmp > 1.f) thetatmp = 1.f;
  theta = acosf (thetatmp);

  ScalorType scalor = 2 * param[k] * (theta - param[theta0]);
  scalor *= sqrtc00i;
  ex *= 1e6;
  ey *= 1e6;
  ez *= 1e6;

  ScalorType tmp = sqrtf (ex*ex + ey*ey + ez*ez);
  if (tmp != 0.f) {
    scalor *= -1. / tmp;
    *f0x = ex * scalor;
    *f0y = ey * scalor;
    *f0z = ez * scalor;
  }
  else {
    *f0x = -scalor;
    *f0y = 0.f;
    *f0z = 0.f;
  }
  
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


__device__ ScalorType AngleHarmonic::
forcePoten0 (const ScalorType * param,
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

  ScalorType theta;
  ScalorType thetatmp = c01 * sqrtc11i * sqrtc00i;
  if (thetatmp > 1.f) thetatmp = 1.f;
  theta = acosf (thetatmp);

  ScalorType scalor = 2 * param[k] * (theta - param[theta0]);
  scalor *= sqrtc00i;
  ex *= 1e6;
  ey *= 1e6;
  ez *= 1e6;
  
  ScalorType tmp = sqrtf ((ex*ex + ey*ey + ez*ez));
  if (tmp != 0.f){
    scalor *= -1.f / tmp;
    *f0x = ex * scalor;
    *f0y = ey * scalor;
    *f0z = ez * scalor;
  }
  else{
    *f0x = -scalor;
    *f0y = 0.f;
    *f0z = 0.f;
  }

  // ScalorType scalor = 2.f * param[k] * (theta - param[theta0]);
  // scalor *= sqrtc00i;
  // ex *= 1e6;
  // ey *= 1e6;
  // ez *= 1e6;
  // ScalorType tmp = sqrtf (ex*ex + ey*ey + ez*ez);
  // scalor *= -1.f / tmp;
  // *f0x = ex * scalor;
  // *f0y = ey * scalor;
  // *f0z = ez * scalor;

  
  ////////////////////////////////////////////////////////////////
  
  // if (tmp != 0){
  //   scalor *= -1.f / tmp;
  //   *f0x = ex * scalor;
  //   *f0y = ey * scalor;
  //   *f0z = ez * scalor;
  // }
  // else {
  //   *f0x = -1.f * scalor;
  //   *f0y = 0.f;
  //   *f0z = 0.f;
  // }

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


__device__ void AngleHarmonic::
force1 (const ScalorType * param,
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
  
  ScalorType theta;
  ScalorType thetatmp = c01 * sqrtc11i * sqrtc00i;
  if (thetatmp > 1.f) thetatmp = 1.f;
  theta = acosf (thetatmp);

  ScalorType scalor = 2 * param[k] * (theta - param[theta0]);
  ex *= 1e6f;
  ey *= 1e6f;
  ez *= 1e6f;
  gx *= 1e6f;
  gy *= 1e6f;
  gz *= 1e6f;

  ScalorType tmp0 = sqrtf (ex*ex + ey*ey + ez*ez);
  ScalorType tmp1 = sqrtf (gx*gx + gy*gy + gz*gz);
  ScalorType scalor0 = sqrtc00i * scalor;
  ScalorType scalor1 = sqrtc11i * scalor;
  
  *f1x = *f1y = *f1z = 0.f;
  
  if (tmp0 != 0){
    scalor0 *= ( 1.f / tmp0 );
    *f1x += scalor0 * ex;
    *f1y += scalor0 * ey;
    *f1z += scalor0 * ez;
  }
  else {
    *f1x += scalor0;
  }
  if (tmp1 != 0){
    scalor1 *= ( 1.f / tmp1 );
    *f1x -= scalor1 * gx;
    *f1y -= scalor1 * gy;
    *f1z -= scalor1 * gz;
  }
  else {
    *f1x -= scalor1;
  }  

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


__device__ ScalorType AngleHarmonic::
forcePoten1 (const ScalorType * param,
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
  
  ScalorType theta;
  ScalorType thetatmp = c01 * sqrtc11i * sqrtc00i;
  if (thetatmp > 1.f) thetatmp = 1.f;
  theta = acosf (thetatmp);
  
  ScalorType scalor = 2 * param[k] * (theta - param[theta0]);
  ex *= 1e6f;
  ey *= 1e6f;
  ez *= 1e6f;
  gx *= 1e6f;
  gy *= 1e6f;
  gz *= 1e6f;

  ScalorType tmp0 = sqrtf (ex*ex + ey*ey + ez*ez);
  ScalorType tmp1 = sqrtf (gx*gx + gy*gy + gz*gz);
  ScalorType scalor0 = sqrtc00i * scalor;
  ScalorType scalor1 = sqrtc11i * scalor;
  
  *f1x = *f1y = *f1z = 0.f;
  
  if (tmp0 != 0){
    scalor0 *= ( 1.f / tmp0 );
    *f1x += scalor0 * ex;
    *f1y += scalor0 * ey;
    *f1z += scalor0 * ez;
  }
  else {
    *f1x += scalor0;
  }
  if (tmp1 != 0){
    scalor1 *= ( 1.f / tmp1 );
    *f1x -= scalor1 * gx;
    *f1y -= scalor1 * gy;
    *f1z -= scalor1 * gz;
  }
  else {
    *f1x -= scalor1;
  }  

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



__device__ void AngleHarmonic::
force01 (const ScalorType * param,
	 const ScalorType diff0x,
	 const ScalorType diff0y,
	 const ScalorType diff0z,
	 const ScalorType diff1x,
	 const ScalorType diff1y,
	 const ScalorType diff1z,
	 ScalorType * f0x,
	 ScalorType * f0y,
	 ScalorType * f0z,
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
  
  ScalorType theta;
  ScalorType thetatmp = c01 * sqrtc11i * sqrtc00i;
  if (thetatmp > 1.f) thetatmp = 1.f;
  theta = acosf (thetatmp);

  ScalorType scalor = 2 * param[k] * (theta - param[theta0]);
  ex *= 1e6f;
  ey *= 1e6f;
  ez *= 1e6f;
  gx *= 1e6f;
  gy *= 1e6f;
  gz *= 1e6f;

  ScalorType tmp0 = sqrtf (ex*ex + ey*ey + ez*ez);
  ScalorType tmp1 = sqrtf (gx*gx + gy*gy + gz*gz);
  ScalorType scalor0 = sqrtc00i * scalor;
  ScalorType scalor1 = sqrtc11i * scalor;
  
  *f0x = *f0y = *f0z = 0.f;
  *f1x = *f1y = *f1z = 0.f;
  
  if (tmp0 != 0){
    scalor0 *= ( 1.f / tmp0 );
    *f0x += scalor0 * ex;
    *f0y += scalor0 * ey;
    *f0z += scalor0 * ez;
  }
  else {
    *f0x += scalor0;
  }
  if (tmp1 != 0){
    scalor1 *= ( 1.f / tmp1 );
    *f1x -= scalor1 * gx;
    *f1y -= scalor1 * gy;
    *f1z -= scalor1 * gz;
  }
  else {
    *f1x -= scalor1;
  }  

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


__device__ ScalorType AngleHarmonic::
forcePoten01 (const ScalorType * param,
	      const ScalorType diff0x,
	      const ScalorType diff0y,
	      const ScalorType diff0z,
	      const ScalorType diff1x,
	      const ScalorType diff1y,
	      const ScalorType diff1z,
	      ScalorType * f0x,
	      ScalorType * f0y,
	      ScalorType * f0z,
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
  
  ScalorType theta;
  ScalorType thetatmp = c01 * sqrtc11i * sqrtc00i;
  if (thetatmp > 1.f) thetatmp = 1.f;
  theta = acosf (thetatmp);
  
  ScalorType scalor = 2 * param[k] * (theta - param[theta0]);
  ex *= 1e6f;
  ey *= 1e6f;
  ez *= 1e6f;
  gx *= 1e6f;
  gy *= 1e6f;
  gz *= 1e6f;

  ScalorType tmp0 = sqrtf (ex*ex + ey*ey + ez*ez);
  ScalorType tmp1 = sqrtf (gx*gx + gy*gy + gz*gz);
  ScalorType scalor0 = sqrtc00i * scalor;
  ScalorType scalor1 = sqrtc11i * scalor;
  
  *f0x = *f0y = *f0z = 0.f;
  *f1x = *f1y = *f1z = 0.f;
  
  if (tmp0 != 0){
    scalor0 *= ( 1.f / tmp0 );
    *f0x += scalor0 * ex;
    *f0y += scalor0 * ey;
    *f0z += scalor0 * ez;
  }
  else {
    *f0x += scalor0;
  }
  if (tmp1 != 0){
    scalor1 *= ( 1.f / tmp1 );
    *f1x -= scalor1 * gx;
    *f1y -= scalor1 * gy;
    *f1z -= scalor1 * gz;
  }
  else {
    *f1x -= scalor1;
  }  

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





















// __device__ void angleForce0 (const InteractionType ftype,
// 			     ScalorType * param,
// 			     const ScalorType diff0x,
// 			     const ScalorType diff0y,
// 			     const ScalorType diff0z,
// 			     const ScalorType diff1x,
// 			     const ScalorType diff1y,
// 			     const ScalorType diff1z,
// 			     ScalorType * f0x,
// 			     ScalorType * f0y,
// 			     ScalorType * f0z)
// {
//   if (ftype == mdForceAngleHarmonic){
//     AngleHarmonic::force0 (param,
// 			   diff0x, diff0y, diff0z, diff1x, diff1y, diff1z,
// 			   f0x, f0y, f0z);
//   }
//   if (ftype == mdForceCosAngle0){
//     CosAngle0::force0 (param,
// 		       diff0x, diff0y, diff0z, diff1x, diff1y, diff1z,
// 		       f0x, f0y, f0z);
//   }
// }

// __device__ void angleForcePoten0 (const InteractionType ftype,
// 				  ScalorType * param,
// 				  const ScalorType diff0x,
// 				  const ScalorType diff0y,
// 				  const ScalorType diff0z,
// 				  const ScalorType diff1x,
// 				  const ScalorType diff1y,
// 				  const ScalorType diff1z,
// 				  ScalorType * f0x,
// 				  ScalorType * f0y,
// 				  ScalorType * f0z,
// 				  ScalorType * dp)
// {
//   if (ftype == mdForceAngleHarmonic){
//     *dp = AngleHarmonic::forcePoten0 (param,
// 				      diff0x, diff0y, diff0z,
// 				      diff1x, diff1y, diff1z,
// 				      f0x, f0y, f0z);
//   }
//   if (ftype == mdForceCosAngle0){
//     *dp = CosAngle0::forcePoten0 (param,
// 				  diff0x, diff0y, diff0z,
// 				  diff1x, diff1y, diff1z,
// 				  f0x, f0y, f0z);
//   }
// }


// __device__ void angleForce1 (const InteractionType ftype,
// 			     ScalorType * param,
// 			     const ScalorType diff0x,
// 			     const ScalorType diff0y,
// 			     const ScalorType diff0z,
// 			     const ScalorType diff1x,
// 			     const ScalorType diff1y,
// 			     const ScalorType diff1z,
// 			     ScalorType * f1x,
// 			     ScalorType * f1y,
// 			     ScalorType * f1z)
// {
//   if (ftype == mdForceAngleHarmonic){
//     AngleHarmonic::force1 (param,
// 			   diff0x, diff0y, diff0z, diff1x, diff1y, diff1z,
// 			   f1x, f1y, f1z);
//   }
//   if (ftype == mdForceCosAngle0){
//     CosAngle0::force1 (param,
// 		       diff0x, diff0y, diff0z, diff1x, diff1y, diff1z,
// 		       f1x, f1y, f1z);
//   }
// }

// __device__ void angleForcePoten1 (const InteractionType ftype,
// 				  ScalorType * param,
// 				  const ScalorType diff0x,
// 				  const ScalorType diff0y,
// 				  const ScalorType diff0z,
// 				  const ScalorType diff1x,
// 				  const ScalorType diff1y,
// 				  const ScalorType diff1z,
// 				  ScalorType * f1x,
// 				  ScalorType * f1y,
// 				  ScalorType * f1z,
// 				  ScalorType * dp)
// {
//   if (ftype == mdForceAngleHarmonic){
//     *dp = AngleHarmonic::forcePoten1 (param,
// 				      diff0x, diff0y, diff0z,
// 				      diff1x, diff1y, diff1z,
// 				      f1x, f1y, f1z);
//   }
//   if (ftype == mdForceCosAngle0){
//     *dp = CosAngle0::forcePoten1 (param,
// 				  diff0x, diff0y, diff0z,
// 				  diff1x, diff1y, diff1z,
// 				  f1x, f1y, f1z);
//   }
// }















inline __host__ void CosAngle0::
initParameter (ScalorType * param,
	       ScalorType k_)
{
  param[k] = k_;
}

__device__ void CosAngle0::
force0 (const ScalorType * param,
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
  ScalorType scalor = param[k] * (-tmp);
  
  *f0x = scalor * (diff1x - c01 * sqrtc00i * sqrtc00i * diff0x);
  *f0y = scalor * (diff1y - c01 * sqrtc00i * sqrtc00i * diff0y);
  *f0z = scalor * (diff1z - c01 * sqrtc00i * sqrtc00i * diff0z);
}


__device__ ScalorType CosAngle0::forcePoten0 (const ScalorType * param,
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
  ScalorType scalor = param[k] * (-tmp);
  
  *f0x = scalor * (diff1x - c01 * sqrtc00i * sqrtc00i * diff0x);
  *f0y = scalor * (diff1y - c01 * sqrtc00i * sqrtc00i * diff0y);
  *f0z = scalor * (diff1z - c01 * sqrtc00i * sqrtc00i * diff0z);

  return param[k] * (1.f - costheta);
}


__device__ void CosAngle0::
force1 (const ScalorType * param,
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
  ScalorType scalor = param[k] * tmp;
  
  *f1x = scalor * (diff1x - c01 * sqrtc00i * sqrtc00i * diff0x - diff0x + c01 * sqrtc11i * sqrtc11i * diff1x);
  *f1y = scalor * (diff1y - c01 * sqrtc00i * sqrtc00i * diff0y - diff0y + c01 * sqrtc11i * sqrtc11i * diff1y);
  *f1z = scalor * (diff1z - c01 * sqrtc00i * sqrtc00i * diff0z - diff0z + c01 * sqrtc11i * sqrtc11i * diff1z);
}


__device__ ScalorType CosAngle0::
forcePoten1 (const ScalorType * param,
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
  ScalorType scalor = param[k] * tmp;
  
  *f1x = scalor * (diff1x - c01 * sqrtc00i * sqrtc00i * diff0x - diff0x + c01 * sqrtc11i * sqrtc11i * diff1x);
  *f1y = scalor * (diff1y - c01 * sqrtc00i * sqrtc00i * diff0y - diff0y + c01 * sqrtc11i * sqrtc11i * diff1y);
  *f1z = scalor * (diff1z - c01 * sqrtc00i * sqrtc00i * diff0z - diff0z + c01 * sqrtc11i * sqrtc11i * diff1z);

  return param[k] * (1.f - costheta) ;  
}


__device__ void CosAngle0::
force01 (const ScalorType * param,
	 const ScalorType diff0x,
	 const ScalorType diff0y,
	 const ScalorType diff0z,
	 const ScalorType diff1x,
	 const ScalorType diff1y,
	 const ScalorType diff1z,
	 ScalorType * f0x,
	 ScalorType * f0y,
	 ScalorType * f0z,
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
  ScalorType scalor = param[k] * tmp;

  *f0x = scalor * (diff1x - c01 * sqrtc00i * sqrtc00i * diff0x);
  *f0y = scalor * (diff1y - c01 * sqrtc00i * sqrtc00i * diff0y);
  *f0z = scalor * (diff1z - c01 * sqrtc00i * sqrtc00i * diff0z);
  
  *f1x = scalor * ( - diff0x + c01 * sqrtc11i * sqrtc11i * diff1x);
  *f1y = scalor * ( - diff0y + c01 * sqrtc11i * sqrtc11i * diff1y);
  *f1z = scalor * ( - diff0z + c01 * sqrtc11i * sqrtc11i * diff1z);
}


__device__ ScalorType CosAngle0::
forcePoten01 (const ScalorType * param,
	      const ScalorType diff0x,
	      const ScalorType diff0y,
	      const ScalorType diff0z,
	      const ScalorType diff1x,
	      const ScalorType diff1y,
	      const ScalorType diff1z,
	      ScalorType * f0x,
	      ScalorType * f0y,
	      ScalorType * f0z,
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
  ScalorType scalor = param[k] * tmp;
  
  *f0x = scalor * (diff1x - c01 * sqrtc00i * sqrtc00i * diff0x);
  *f0y = scalor * (diff1y - c01 * sqrtc00i * sqrtc00i * diff0y);
  *f0z = scalor * (diff1z - c01 * sqrtc00i * sqrtc00i * diff0z);

  *f1x = scalor * (- diff0x + c01 * sqrtc11i * sqrtc11i * diff1x);
  *f1y = scalor * (- diff0y + c01 * sqrtc11i * sqrtc11i * diff1y);
  *f1z = scalor * (- diff0z + c01 * sqrtc11i * sqrtc11i * diff1z);

  return param[k] * (1.f - costheta) ;  
}



#endif
