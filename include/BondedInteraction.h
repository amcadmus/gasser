#ifndef __BondedInteraction_h_wanghan__
#define __BondedInteraction_h_wanghan__

#include "common.h"

enum mdBondInteraction {
  mdForceHarmonicSpring		= 0,
  mdForceFENE			= 1
} ;
typedef enum mdBondInteraction mdBondInteraction_t;

enum mdBondInteractionNParam{
  mdForceNParamHarmonicSpring	= 2,
  mdForceNParamFENE		= 2
};

inline __host__ IndexType calNumBondParameter (mdBondInteraction_t type)
{
  switch (type){
  case mdForceHarmonicSpring:
      return mdForceNParamHarmonicSpring;
  case mdForceFENE:
      return mdForceNParamFENE;
  default:
      return 0;
  }
}

namespace HarmonicSpring {
    typedef enum paramIndex{
      epsilon		= 0,
      r0		= 1
    } paramIndex_t;
    
    __host__ void initParameter (ScalorType * param,
				 ScalorType k,
				 ScalorType r0);
    __device__ void force (const ScalorType * param,
			   ScalorType diffx,
			   ScalorType diffy,
			   ScalorType diffz,
			   ScalorType *fx, 
			   ScalorType *fy,
			   ScalorType *fz);
    __device__ ScalorType forcePoten (const ScalorType * param,
				      ScalorType diffx,
				      ScalorType diffy,
				      ScalorType diffz,
				      ScalorType *fx, 
				      ScalorType *fy,
				      ScalorType *fz);
};

namespace FENE {
    typedef enum paramIndex{
      epsilon		= 0,
      rinf2		= 1
    } paramIndex_t;

    __host__ void initParameter (ScalorType * param,
				 ScalorType k,
				 ScalorType rinf);
    __device__ void force (const ScalorType * param,
			   ScalorType diffx,
			   ScalorType diffy,
			   ScalorType diffz,
			   ScalorType *fx, 
			   ScalorType *fy,
			   ScalorType *fz);
    __device__ ScalorType forcePoten (const ScalorType * param,
				      ScalorType diffx,
				      ScalorType diffy,
				      ScalorType diffz,
				      ScalorType *fx, 
				      ScalorType *fy,
				      ScalorType *fz);
};


__device__ void 
bondForce (const mdBondInteraction_t ftype,
	   const ScalorType * param,
	   ScalorType diffx, ScalorType diffy, ScalorType diffz,
	   ScalorType *fx,   ScalorType *fy,   ScalorType *fz)
{
  if (ftype == mdForceHarmonicSpring){
      HarmonicSpring::force (param, diffx, diffy, diffz, fx, fy, fz);
  }
  if (ftype == mdForceFENE){
      FENE::force (param, diffx, diffy, diffz, fx, fy, fz);
  }
  // switch (ftype){
  // case mdForceHarmonicSpring:
  //     HarmonicSpring::force (param, diffx, diffy, diffz, fx, fy, fz);
  //     break;
  // case mdForceFENE:
  //     FENE::force (param, diffx, diffy, diffz, fx, fy, fz);
  //     break;
  // default:
  //     *fx = *fy = *fz = 0.f;
  //     break;
  // }
}


__device__ void
bondForcePoten (const mdBondInteraction_t ftype,
		const ScalorType * param,
		ScalorType diffx, ScalorType diffy, ScalorType diffz,
		ScalorType *fx,   ScalorType *fy,   ScalorType *fz,
		ScalorType * dp)
{
  if (ftype == mdForceHarmonicSpring){
    *dp = HarmonicSpring::forcePoten (param, diffx, diffy, diffz, fx, fy, fz);
  }
  if (ftype == mdForceFENE){
    *dp = FENE::forcePoten (param,  diffx, diffy, diffz, fx, fy, fz);
  }  
  // switch (ftype){
  // case mdForceHarmonicSpring:
  //     return HarmonicSpring::forcePoten (param, diffx, diffy, diffz, fx, fy, fz);
  // case mdForceFENE:
  //     return FENE::forcePoten (param,  diffx, diffy, diffz, fx, fy, fz);
  // default:
  //     *fx = *fy = *fz = 0.f;
  //     return 0.f;
  // }
}





















inline __host__ void HarmonicSpring::initParameter (ScalorType *param,
						    ScalorType ep_,
						    ScalorType r0_)
{
  param[epsilon] = ep_;
  param[r0] = r0_;
}

__device__ void HarmonicSpring::force (const ScalorType * param,
				       ScalorType diffx,
				       ScalorType diffy,
				       ScalorType diffz,
				       ScalorType *fx, 
				       ScalorType *fy,
				       ScalorType *fz)
{
  ScalorType dr = __fsqrt_rn(diffx*diffx + diffy*diffy + diffz*diffz);
  ScalorType dri = __frcp_rn(dr);

  ScalorType scalor = - param[epsilon] * (1.f - param[r0] * dri) ;
  *fx = diffx * scalor;
  *fy = diffy * scalor;
  *fz = diffz * scalor;
}


__device__ ScalorType HarmonicSpring::forcePoten (const ScalorType * param,
						  ScalorType diffx,
						  ScalorType diffy,
						  ScalorType diffz,
						  ScalorType *fx, 
						  ScalorType *fy,
						  ScalorType *fz)
{
  ScalorType dr = __fsqrt_rn(diffx*diffx + diffy*diffy + diffz*diffz);
  ScalorType dri = __frcp_rn(dr);

  ScalorType scalor = - param[epsilon] * (1.f - param[r0] * dri) ;
  *fx = diffx * scalor;
  *fy = diffy * scalor;
  *fz = diffz * scalor;

  return 0.5f * param[epsilon] * (dr - param[r0]) * (dr - param[r0]);
}


inline __host__ void FENE::initParameter (ScalorType * param,
					  ScalorType ep,
					  ScalorType rinf)
{
  param[epsilon] = ep;
  param[rinf2] = rinf * rinf;
}

__device__  void FENE::force (const ScalorType * param,
				       ScalorType diffx,
				       ScalorType diffy,
				       ScalorType diffz,
				       ScalorType *fx, 
				       ScalorType *fy,
			      ScalorType *fz)
{
  ScalorType dr2 = diffx*diffx + diffy*diffy + diffz*diffz;

  ScalorType scalor = - param[epsilon] * __fdiv_rn(param[rinf2], param[rinf2] - dr2);
  *fx = diffx * scalor;
  *fy = diffy * scalor;
  *fz = diffz * scalor;
}
  
__device__ ScalorType FENE::forcePoten (const ScalorType * param,
					ScalorType diffx,
					ScalorType diffy,
					ScalorType diffz,
					ScalorType *fx, 
					ScalorType *fy,
					ScalorType *fz)
{
  ScalorType dr2 = diffx*diffx + diffy*diffy + diffz*diffz;

  ScalorType core =  __fdiv_rn(param[rinf2], param[rinf2] - dr2);
  ScalorType scalor = - core * param[epsilon];
  *fx = diffx * scalor;
  *fy = diffy * scalor;
  *fz = diffz * scalor;
  ScalorType tmp = __frcp_rn(core);
  if (tmp <= 0.f){
    return 1e10;
  }
  else{
    return -0.5f * param[epsilon] * param[rinf2] * __logf (tmp);
  }
}
  





    

struct HarmonicSpringParameters
{
  ScalorType k;
  ScalorType r0;
};

__device__ void initHarmonicSpringParameters (HarmonicSpringParameters * param,
					      ScalorType k_,
					      ScalorType r0_);
__device__ void harmonicSpringForce (HarmonicSpringParameters param,
				     ScalorType diffx,
				     ScalorType diffy,
				     ScalorType diffz,
				     ScalorType *fx, ScalorType *fy, ScalorType *fz);
__device__ ScalorType harmonicSpringForcePoten (HarmonicSpringParameters param,
						ScalorType diffx,
						ScalorType diffy,
						ScalorType diffz,
						ScalorType *fx, 
						ScalorType *fy, 
						ScalorType *fz);

struct FENEBondParameters
{
  ScalorType k;
  ScalorType rinf2;
};
__device__ void initFENEBondParameters (FENEBondParameters * param,
					ScalorType k_,
					ScalorType rinf_);
__device__ void FENEBondForce (FENEBondParameters  param,
			       ScalorType diffx,
			       ScalorType diffy,
			       ScalorType diffz,
			       ScalorType *fx, 
			       ScalorType *fy, 
			       ScalorType *fz);
__device__ ScalorType ScalorTypeFENEBondForcePoten (FENEBondParameters  param,
						    ScalorType diffx,
						    ScalorType diffy,
						    ScalorType diffz,
						    ScalorType *fx, 
						    ScalorType *fy, 
						    ScalorType *fz);



////////////////////////////////////////////////////////////
// implementation Harmonic Spring
// /////////////////////////////////////////////////////////
__device__ void initHarmonicSpringParameters (HarmonicSpringParameters * param,
					      ScalorType k_,
					      ScalorType r0_)
{
  param->k = k_;
  param->r0 = r0_;
}

__device__ void harmonicSpringForce (HarmonicSpringParameters  param,
				     ScalorType diffx,
				     ScalorType diffy,
				     ScalorType diffz,
				     ScalorType *fx, ScalorType *fy, ScalorType *fz)
{
  ScalorType dr = __fsqrt_rn(diffx*diffx + diffy*diffy + diffz*diffz);
  ScalorType dri = __frcp_rn(dr);

  ScalorType scalor = - param.k * (1.f - param.r0 * dri) ;
  *fx = diffx * scalor;
  *fy = diffy * scalor;
  *fz = diffz * scalor;
}

__device__ ScalorType harmonicSpringForcePoten (HarmonicSpringParameters  param,
						ScalorType diffx,
						ScalorType diffy,
						ScalorType diffz,
						ScalorType *fx, 
						ScalorType *fy, 
						ScalorType *fz)
{
  ScalorType dr = __fsqrt_rn(diffx*diffx + diffy*diffy + diffz*diffz);
  ScalorType dri = __frcp_rn(dr);

  ScalorType scalor = - param.k * (1.f - param.r0 * dri) ;
  *fx = diffx * scalor;
  *fy = diffy * scalor;
  *fz = diffz * scalor;

  return 0.5f * param.k * (dr - param.r0) * (dr - param.r0);
}


////////////////////////////////////////////////////////////
// implementation  fene bond
// /////////////////////////////////////////////////////////
__device__ void initFENEBondParameters (FENEBondParameters * param,
					ScalorType k_,
					ScalorType rinf_)
{
  param->k = k_;
  param->rinf2 = rinf_ * rinf_;
}

__device__ void FENEBondForce (FENEBondParameters  param,
			       ScalorType diffx,
			       ScalorType diffy,
			       ScalorType diffz,
			       ScalorType *fx, 
			       ScalorType *fy, 
			       ScalorType *fz)
{
  ScalorType dr2 = diffx*diffx + diffy*diffy + diffz*diffz;

  ScalorType scalor = - param.k * __fdiv_rn(param.rinf2, param.rinf2 - dr2);
  *fx = diffx * scalor;
  *fy = diffy * scalor;
  *fz = diffz * scalor;
}

__device__ ScalorType FENEBondForcePoten (FENEBondParameters  param,
					  ScalorType diffx,
					  ScalorType diffy,
					  ScalorType diffz,
					  ScalorType *fx, 
					  ScalorType *fy, 
					  ScalorType *fz)
{
  ScalorType dr2 = diffx*diffx + diffy*diffy + diffz*diffz;

  ScalorType core =  __fdiv_rn(param.rinf2, param.rinf2 - dr2);
  ScalorType scalor = - core * param.k;
  *fx = diffx * scalor;
  *fy = diffy * scalor;
  *fz = diffz * scalor;

  return -0.5f * param.k * param.rinf2 * logf (__fdiv_rn(param.rinf2 - dr2, param.rinf2));
}





#endif
