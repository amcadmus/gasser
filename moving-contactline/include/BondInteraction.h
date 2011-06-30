#ifndef __BondInteraction_h_wanghan__
#define __BondInteraction_h_wanghan__

#include "common.h"
#include "Interaction.h"

enum mdBondInteractionNParam{
  mdForceNParamHarmonicSpring	= 2,
  mdForceNParamFENE		= 2,
  mdForceNParamFENE2		= 3
};


class HarmonicSpringParameter : public BondInteractionParameter
{
  ScalorType param [mdForceNParamHarmonicSpring];
public:
  HarmonicSpringParameter () {}
  HarmonicSpringParameter (ScalorType k,
			   ScalorType r0);
  void reinit (ScalorType k,
	       ScalorType r0);
  virtual InteractionType type () const;
  virtual unsigned numParam () const ;
  virtual ScalorType * c_ptr () ;
  virtual const ScalorType * c_ptr () const ;
};

class FENEParameter : public BondInteractionParameter
{
  ScalorType param [mdForceNParamFENE];
public:
  FENEParameter () {}
  FENEParameter (ScalorType k,
		 ScalorType rinf);
  void reinit (ScalorType k,
	       ScalorType rinf);
  virtual InteractionType type () const;
  virtual unsigned numParam () const ;
  virtual ScalorType * c_ptr () ;
  virtual const ScalorType * c_ptr () const ;
};

class FENE2Parameter : public BondInteractionParameter
{
  ScalorType param [mdForceNParamFENE2];
public:
  FENE2Parameter () {}
  FENE2Parameter (ScalorType k,
		  ScalorType rs,
		  ScalorType r0);
  void reinit (ScalorType k,
	       ScalorType rs,
	       ScalorType r0);
  virtual InteractionType type () const;
  virtual unsigned numParam () const ;
  virtual ScalorType * c_ptr () ;
  virtual const ScalorType * c_ptr () const ;
};


namespace HarmonicSpring {
  typedef enum paramIndex{
    epsilon		= 0,
    r0		= 1
  } paramIndex_t;
    
  __host__ void
  initParameter (ScalorType * param,
		 ScalorType k,
		 ScalorType r0);
  static __device__ void
  force (const ScalorType * param,
	 ScalorType diffx,
	 ScalorType diffy,
	 ScalorType diffz,
	 ScalorType *fx, 
	 ScalorType *fy,
	 ScalorType *fz);
  static __device__ ScalorType
  forcePoten (const ScalorType * param,
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

  __host__ void
  initParameter (ScalorType * param,
		 ScalorType k,
		 ScalorType rinf);
  static __device__ void
  force (const ScalorType * param,
	 ScalorType diffx,
	 ScalorType diffy,
	 ScalorType diffz,
	 ScalorType *fx, 
	 ScalorType *fy,
	 ScalorType *fz);
  static __device__ ScalorType
  forcePoten (const ScalorType * param,
	      ScalorType diffx,
	      ScalorType diffy,
	      ScalorType diffz,
	      ScalorType *fx, 
	      ScalorType *fy,
	      ScalorType *fz);
};

namespace FENE2 {
  typedef enum paramIndex{
    epsilon		= 0,
    rs2		= 1,
    r0p		= 2
  } paramIndex_t;

  __host__ void
  initParameter (ScalorType * param,
		 ScalorType k,
		 ScalorType rs,
		 ScalorType r0);
  static __device__ void
  force (const ScalorType * param,
	 ScalorType diffx,
	 ScalorType diffy,
	 ScalorType diffz,
	 ScalorType *fx, 
	 ScalorType *fy,
	 ScalorType *fz);
  static __device__ ScalorType
  forcePoten (const ScalorType * param,
	      ScalorType diffx,
	      ScalorType diffy,
	      ScalorType diffz,
	      ScalorType *fx, 
	      ScalorType *fy,
	      ScalorType *fz);
};



static __device__ void 
bondForce (const InteractionType ftype,
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
  if (ftype == mdForceFENE2){
    FENE2::force (param, diffx, diffy, diffz, fx, fy, fz);
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


static __device__ void
bondForcePoten (const InteractionType ftype,
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
  if (ftype == mdForceFENE2){
    *dp = FENE2::forcePoten (param,  diffx, diffy, diffz, fx, fy, fz);
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



inline __host__ void HarmonicSpring::
initParameter (ScalorType *param,
	       ScalorType ep_,
	       ScalorType r0_)
{
  param[epsilon] = ep_;
  param[r0] = r0_;
}

__device__ void HarmonicSpring::
force (const ScalorType * param,
       ScalorType diffx,
       ScalorType diffy,
       ScalorType diffz,
       ScalorType *fx, 
       ScalorType *fy,
       ScalorType *fz)
{
//   ScalorType dr = __fsqrt_rn(diffx*diffx + diffy*diffy + diffz*diffz);
//   ScalorType dri = __frcp_rn(dr);
  ScalorType dri = 1./sqrtf(diffx*diffx + diffy*diffy + diffz*diffz);
  
  ScalorType scalor = param[epsilon] * (1.f - param[r0] * dri) ;
  *fx = diffx * scalor;
  *fy = diffy * scalor;
  *fz = diffz * scalor;
}


__device__ ScalorType
HarmonicSpring::forcePoten (const ScalorType * param,
			    ScalorType diffx,
			    ScalorType diffy,
			    ScalorType diffz,
			    ScalorType *fx, 
			    ScalorType *fy,
			    ScalorType *fz)
{
//   ScalorType dr = __fsqrt_rn(diffx*diffx + diffy*diffy + diffz*diffz);
//   ScalorType dri = __frcp_rn(dr);
  ScalorType dr = diffx*diffx + diffy*diffy + diffz*diffz;
  ScalorType dri = 1./sqrtf(dr);
  dr *= dri;

  ScalorType scalor = param[epsilon] * (1.f - param[r0] * dri) ;
  *fx = diffx * scalor;
  *fy = diffy * scalor;
  *fz = diffz * scalor;

  return 0.5f * param[epsilon] * (dr - param[r0]) * (dr - param[r0]);
}


inline __host__ void FENE::
initParameter (ScalorType * param,
	       ScalorType ep,
	       ScalorType rinf)
{
  param[epsilon] = ep;
  param[rinf2] = rinf * rinf;
}

__device__  void FENE::
force (const ScalorType * param,
       ScalorType diffx,
       ScalorType diffy,
       ScalorType diffz,
       ScalorType *fx, 
       ScalorType *fy,
       ScalorType *fz)
{
  ScalorType dr2 = diffx*diffx + diffy*diffy + diffz*diffz;

  ScalorType scalor = param[epsilon] * __fdiv_rn(param[rinf2], param[rinf2] - dr2);
  *fx = diffx * scalor;
  *fy = diffy * scalor;
  *fz = diffz * scalor;
}
  
__device__ ScalorType FENE::
forcePoten (const ScalorType * param,
	    ScalorType diffx,
	    ScalorType diffy,
	    ScalorType diffz,
	    ScalorType *fx, 
	    ScalorType *fy,
	    ScalorType *fz)
{
  ScalorType dr2 = diffx*diffx + diffy*diffy + diffz*diffz;

  ScalorType core =  __fdiv_rn(param[rinf2], param[rinf2] - dr2);
  ScalorType scalor = core * param[epsilon];
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



inline __host__ void FENE2::
initParameter (ScalorType * param,
	       ScalorType ep,
	       ScalorType rs,
	       ScalorType r0_)
{
  param[epsilon] = ep;
  param[rs2] = rs * rs;
  param[r0p] = r0_;
}

__device__  void FENE2::
force (const ScalorType * param,
       ScalorType diffx,
       ScalorType diffy,
       ScalorType diffz,
       ScalorType *fx, 
       ScalorType *fy,
       ScalorType *fz)
{
  ScalorType dr = sqrtf(diffx*diffx + diffy*diffy + diffz*diffz);
  ScalorType r_r0 = dr - param[r0p];
  
  ScalorType scalor = param[epsilon] * param[rs2] * r_r0 / (dr * (param[rs2] - r_r0*r_r0));
  *fx = diffx * scalor;
  *fy = diffy * scalor;
  *fz = diffz * scalor;
}
  
__device__ ScalorType FENE2::
forcePoten (const ScalorType * param,
	    ScalorType diffx,
	    ScalorType diffy,
	    ScalorType diffz,
	    ScalorType *fx, 
	    ScalorType *fy,
	    ScalorType *fz)
{
  ScalorType dr = sqrtf(diffx*diffx + diffy*diffy + diffz*diffz);
  ScalorType r_r0 = dr - param[r0p];
  ScalorType core = param[rs2] - r_r0*r_r0;

  ScalorType scalor = param[epsilon] * param[rs2] * r_r0 / (core * dr);
  
  *fx = diffx * scalor;
  *fy = diffy * scalor;
  *fz = diffz * scalor;
  // printf ("dr is %f, scalor is %f\n", dr, scalor);
  
  ScalorType tmp = core / param[rs2];
  if (tmp <= 0.f){
    return 1e10;
  }
  else{
    return -0.5f* param[epsilon] * param[rs2] * __logf (tmp);
  }
}




#endif
