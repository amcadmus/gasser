#ifndef __tmp_h_wanghan__
#define __tmp_h_wanghan__

#include "common.h"
#include <math.h>

namespace tmptest{
    
    struct LJ6_12_Parameters
    {
      ScalorType epsilon;
      ScalorType sigma;
      ScalorType shift;
      ScalorType rcut;
    }
	;

    void initLJ6_12_Parameters (LJ6_12_Parameters * param,
				ScalorType rcut = 2.5,
				ScalorType epsilon = 1.,
				ScalorType sigma = 1.,
				ScalorType shift = 0.);
//  ScalorType LJ6_12_Potential (LJ6_12_Parameters param,
// 					ScalorType diffx,
// 					ScalorType diffy, 
// 					ScalorType diffz);
    void LJ6_12_Force (LJ6_12_Parameters param,
		       ScalorType diffx,
		       ScalorType diffy,
		       ScalorType diffz,
		       ScalorType *fx, ScalorType *fy, ScalorType *fz);
    ScalorType LJ6_12_Force_Poten (LJ6_12_Parameters param,
				   ScalorType diffx,
				   ScalorType diffy,
				   ScalorType diffz,
				   ScalorType *fx, ScalorType *fy, ScalorType *fz);

    struct LJ6_12_cap_Parameters
    {
      ScalorType epsilon;
      ScalorType sigma;
      ScalorType shift;
      ScalorType rcut;
      ScalorType cap;
      ScalorType capR;
      ScalorType pcapR;
    }
	;

    void initLJ6_12_cap_Parameters (LJ6_12_cap_Parameters * param,
				    ScalorType rcut = 2.5,
				    ScalorType epsilon = 1.,
				    ScalorType sigma = 1.,
				    ScalorType shift = 0.,
				    ScalorType cap = 1000);
    void LJ6_12_cap_Force (LJ6_12_cap_Parameters param,
			   ScalorType diffx,
			   ScalorType diffy,
			   ScalorType diffz,
			   ScalorType *fx, ScalorType *fy, ScalorType *fz);
    ScalorType LJ6_12_cap_Force_Poten (LJ6_12_cap_Parameters param,
				       ScalorType diffx,
				       ScalorType diffy,
				       ScalorType diffz,
				       ScalorType *fx, ScalorType *fy, ScalorType *fz);


    ScalorType lj6_12_originForce (LJ6_12_cap_Parameters param, ScalorType r);
    
    ScalorType lj6_12_originPoten (LJ6_12_cap_Parameters param, ScalorType r);
    
    ScalorType findCapR (LJ6_12_cap_Parameters param);
    

}
;



////////////////////////////////////////////////////////////
// implementation lj6-12
// /////////////////////////////////////////////////////////
void tmptest::initLJ6_12_Parameters (LJ6_12_Parameters * param,
				     ScalorType rcut,
				     ScalorType epsilon,
				     ScalorType sigma,
				     ScalorType shift)
{
  param->epsilon = epsilon;
  param->sigma = sigma;
  param->shift = shift;
  param->rcut = rcut;
}

//  ScalorType LJ6_12_Potential (LJ6_12_Parameters param,
// 					ScalorType diffx,
// 					ScalorType diffy, 
// 					ScalorType diffz)
// {
//   ScalorType dr2 = diffx*diffx + diffy*diffy + diffz*diffz;
//   if (dr2 > param.rcut*param.rcut) return 0;
//   ScalorType sri2 = param.sigma*param.sigma / dr2;
//   ScalorType sri6 = sri2*sri2*sri2;
//   return 4 * param.epsilon * (sri6*sri6 - sri6 + param.shift);
// }

void tmptest::LJ6_12_Force (LJ6_12_Parameters param,
			    ScalorType diffx,
			    ScalorType diffy,
			    ScalorType diffz,
			    ScalorType *fx, ScalorType *fy, ScalorType *fz)
{
  ScalorType dr2 = diffx*diffx + diffy*diffy + diffz*diffz;
  if (dr2 > param.rcut*param.rcut) {
    *fx = *fy = *fz = 0.;
    return;
  }
  ScalorType ri2 = 1./dr2;
  ScalorType sri2 = param.sigma*param.sigma * ri2;
  ScalorType sri6 = sri2*sri2*sri2;
  ScalorType scalor = 4 * param.epsilon * (12 * sri6*sri6 - 6 * sri6) * ri2;
  *fx = diffx * scalor;
  *fy = diffy * scalor;
  *fz = diffz * scalor;
}


ScalorType tmptest::LJ6_12_Force_Poten (LJ6_12_Parameters param,
					ScalorType diffx,
					ScalorType diffy,
					ScalorType diffz,
					ScalorType *fx, ScalorType *fy, ScalorType *fz)
{
  ScalorType dr2 = diffx*diffx + diffy*diffy + diffz*diffz;
  if (dr2 > param.rcut*param.rcut) {
    *fx = *fy = *fz = 0.;
    return 0.;
  }
  ScalorType ri2 = 1./dr2;
  ScalorType sri2 = param.sigma*param.sigma * ri2;
  ScalorType sri6 = sri2*sri2*sri2;
  ScalorType scalor = 4 * param.epsilon * (12 * sri6*sri6 - 6 * sri6) * ri2;
  *fx = diffx * scalor;
  *fy = diffy * scalor;
  *fz = diffz * scalor;
  return 4 * param.epsilon * (sri6*sri6 - sri6 + param.shift);
}


////////////////////////////////////////////////////////////
// implementation lj6-12-cap
// /////////////////////////////////////////////////////////

 ScalorType tmptest::lj6_12_originForce (LJ6_12_cap_Parameters param, ScalorType r)
{
  ScalorType dr2 = r*r;
  if (dr2 > param.rcut*param.rcut) {
    return 0.;
  }
  ScalorType ri = 1./r;
  ScalorType ri2 = ri*ri;
  ScalorType sri2 = param.sigma*param.sigma * ri2;
  ScalorType sri6 = sri2*sri2*sri2;
  return 4 * param.epsilon * (12 * sri6*sri6 - 6 * sri6) * ri;
}

 ScalorType tmptest::lj6_12_originPoten (LJ6_12_cap_Parameters param, ScalorType r)
{
  ScalorType dr2 = r*r;
  if (dr2 > param.rcut*param.rcut) {
    return 0.;
  }
  ScalorType ri = 1./r;
  ScalorType ri2 = ri*ri;
  ScalorType sri2 = param.sigma*param.sigma * ri2;
  ScalorType sri6 = sri2*sri2*sri2;
  return 4 * param.epsilon * (sri6*sri6 - sri6 + param.shift);
}

 ScalorType tmptest::findCapR (LJ6_12_cap_Parameters param)
{
  ScalorType point = 0.95 * param.rcut;
  while (lj6_12_originForce(param, point) < param.cap){
    point *= 0.95;
  }
  ScalorType lower = point;
  ScalorType upper = param.rcut;
  while (fabs(upper - lower) > 1e-4){
    ScalorType mid = 0.5 * (upper + lower);
    if (lj6_12_originForce(param, mid) > param.cap) {
      lower = mid;
    }
    else {
      upper = mid;
    }
  }
  return 0.5 * (upper + lower);
}

void tmptest::initLJ6_12_cap_Parameters (LJ6_12_cap_Parameters * param,
					 ScalorType rcut,
					 ScalorType epsilon,
					 ScalorType sigma,
					 ScalorType shift,
					 ScalorType cap)
{
  param->epsilon = epsilon;
  param->sigma = sigma;
  param->shift = shift;
  param->rcut = rcut;
  param->cap = cap;
  param->capR = findCapR (*param);
  param->pcapR = lj6_12_originPoten(*param, param->capR);
}


void tmptest::LJ6_12_cap_Force (LJ6_12_cap_Parameters param,
				ScalorType diffx,
				ScalorType diffy,
				ScalorType diffz,
				ScalorType *fx, ScalorType *fy, ScalorType *fz)
{
  ScalorType dr2 = diffx*diffx + diffy*diffy + diffz*diffz;
  if (dr2 > param.rcut*param.rcut) {
    *fx = *fy = *fz = 0.;
    return;
  }
  if (dr2 == 0){
    *fx = param.cap;
    *fy = 0;
    *fz = 0;
    return;
  }
  if (dr2 < param.capR*param.capR){
    ScalorType s = 1. / sqrtf(diffx*diffx + diffy*diffy + diffz*diffz);
    *fx = s * param.cap * diffx;
    *fy = s * param.cap * diffy;
    *fz = s * param.cap * diffz;
    return;
  }
  ScalorType ri2 = 1./dr2;
  ScalorType sri2 = param.sigma*param.sigma * ri2;
  ScalorType sri6 = sri2*sri2*sri2;
  ScalorType scalor = 4 * param.epsilon * (12 * sri6*sri6 - 6 * sri6) * ri2;
  *fx = diffx * scalor;
  *fy = diffy * scalor;
  *fz = diffz * scalor;
}


ScalorType tmptest::LJ6_12_cap_Force_Poten (LJ6_12_cap_Parameters param,
				   ScalorType diffx,
				   ScalorType diffy,
				   ScalorType diffz,
				   ScalorType *fx, ScalorType *fy, ScalorType *fz)
{
  ScalorType dr2 = diffx*diffx + diffy*diffy + diffz*diffz;
  if (dr2 > param.rcut*param.rcut) {
    *fx = *fy = *fz = 0.;
    return 0.;
  }
  if (dr2 < param.capR*param.capR){
    ScalorType s = 1. / sqrtf(diffx*diffx + diffy*diffy + diffz*diffz);
    *fx = s * param.cap * diffx;
    *fy = s * param.cap * diffy;
    *fz = s * param.cap * diffz;
    return -param.cap * (sqrtf(dr2) - param.capR) + param.pcapR ;
  }
  ScalorType ri2 = 1./dr2;
  ScalorType sri2 = param.sigma*param.sigma * ri2;
  ScalorType sri6 = sri2*sri2*sri2;
  ScalorType scalor = 4 * param.epsilon * (12 * sri6*sri6 - 6 * sri6) * ri2;
  *fx = diffx * scalor;
  *fy = diffy * scalor;
  *fz = diffz * scalor;
  return 4 * param.epsilon * (sri6*sri6 - sri6 + param.shift);
}






namespace tmpBoxGeometry{
    struct tmpBox 
    {
      VectorType size;
      VectorType sizei;
    };
    
    void setBoxSize (ScalorType x, ScalorType y, ScalorType z,
		     tmpBox *box);
    void moveParticleToBox (
	tmpBox box,
	ScalorType * x, ScalorType * y, ScalorType * z,
	IntScalorType * noix, 
	IntScalorType * noiy, 
	IntScalorType * noiz);
    void shortestImage (tmpBox box,
			ScalorType * x, ScalorType * y, ScalorType * z);
}
;

void tmpBoxGeometry::setBoxSize (
    ScalorType x, ScalorType y, ScalorType z,
    tmpBoxGeometry::tmpBox *box)
{
  box->size.x = x;
  box->size.y = y;
  box->size.z = z;
  box->sizei.x = 1/x;
  box->sizei.y = 1/y;
  box->sizei.z = 1/z;
}

void tmpBoxGeometry::moveParticleToBox(
    tmpBoxGeometry::tmpBox rectBox,
    ScalorType *x,  
    ScalorType *y, 
    ScalorType *z,
    IntScalorType *noix, 
    IntScalorType *noiy, 
    IntScalorType *noiz)
{
  IntScalorType tmp;
  tmp = floorf(*x * rectBox.sizei.x);
  *noix += tmp;
  *x -= tmp * rectBox.size.x;

  tmp = floorf(*y * rectBox.sizei.y);
  *noiy += tmp;
  *y -= tmp * rectBox.size.y;

  tmp = floorf(*z * rectBox.sizei.z);
  *noiz += tmp;
  *z -= tmp * rectBox.size.z;
}

void tmpBoxGeometry::shortestImage (
    tmpBoxGeometry::tmpBox box,
    ScalorType * x, ScalorType * y, ScalorType * z)
{
  *x -= floorf(*x * box.sizei.x + 0.5) * box.size.x;
  *y -= floorf(*y * box.sizei.y + 0.5) * box.size.y;
  *z -= floorf(*z * box.sizei.z + 0.5) * box.size.z;
}






#endif
