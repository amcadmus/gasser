#ifndef __NonBondedInteraction_h_wanghan__
#define __NonBondedInteraction_h_wanghan__

#include "common.h"
#include <math.h>

#pragma option -b
enum mdNBInteraction {
  mdForceNBNull				= 0,
  mdForceLennardJones6_12		= 1,
  mdForceLennardJones6_12_cap		= 2,
  mdForceCosTail			= 3,
};
typedef  int mdNBInteraction_t;

enum mdNBInteractionNParam {
  mdForceNParamNBNull			= 1,
  mdForceNParamLennardJones6_12		= 4,
  mdForceNParamLennardJones6_12_cap	= 7,
  mdForceNParamCosTail			= 6
};

inline __host__ IndexType calNumNBParameter (mdNBInteraction_t type) 
{
  switch (type){
  case mdForceNBNull:
      return mdForceNParamNBNull;
  case mdForceLennardJones6_12:
      return mdForceNParamLennardJones6_12;
  case mdForceLennardJones6_12_cap:
      return mdForceNParamLennardJones6_12_cap;
  case mdForceCosTail:
      return mdForceNParamCosTail;
  default:
      return 0;
  }
}
  

/**
 * Records all infromation of non-bonded interaction. Including the
 * type and parameters of the interactions.
 * 
 */
struct HostNBForceSetting
{
  mdNBInteraction_t * type;	/**< vector keeps types of non-bonded
				 * interactions */
  IndexType NNBForce;		/**< number of non-bonded interactions */
  ScalorType * param;		/**< vector keeps parameters of interactions */
  IndexType * paramPosi;	/**< i-th value is the start position
				 * of i-th interaction in the vector
				 * param */
  IndexType paramLength;	/**< length of the vector param */
};

/**
 * Maps the type of atoms to the index of non-bonded interaction.
 * 
 */
struct HostAtomNBForceTable
{
  ForceIndexType * data;	/**< data of the mapping table */
  IndexType dataLength;		/**< length of the data */
  IndexType NAtomType;		/**< number of atom type */
};

/**
 * The description of non-bonded interaction in a system.
 * 
 */
class SystemNBForce
{
  IndexType memNBForce;
  IndexType memNBParam;
public:
  HostAtomNBForceTable indexTable;
  HostNBForceSetting   setting;
public:
  SystemNBForce();
  ~SystemNBForce();
  /** 
   * Initialize the non-bonded interaction information in the system
   * 
   * @param NatomType Possible maximum number of atom types in the
   * system.  This value should be larger than or equal to the maximum
   * type (int value) will appear in the following function
   * addNBForce.
   */
  void init (const IndexType NAtomType);
  /** 
   * Add a non-bonded interaction to the system.
   * 
   * @param i One atom type.
   * @param j Another atom type.
   * @param forceType Non-bonded force type.
   * @param param Parameters of the force.
   */
  void addNBForce (const TypeType &i,
		   const TypeType &j, 
		   const mdNBInteraction_t & forceType,
		   const ScalorType * param);
};


  


namespace AtomNBForceTable{
    __host__ IndexType calDataLength (const IndexType Ntype);
    __host__ int deviceInitTable (const ForceIndexType * htableData,
				  ForceIndexType ** dtableData,
				  const IndexType Ntype);
    __host__ void setTableItem (ForceIndexType * tableData,
				const IndexType Ntype,
				const TypeType atom0, 
				const TypeType atom1,
				const ForceIndexType forceType);
    __device__ ForceIndexType calForceIndex (ForceIndexType * tableData, 
					     const IndexType Ntype,
					     const TypeType atom0, 
					     const TypeType atom1);
    __device__ ForceIndexType calForceIndex (volatile ForceIndexType * tableData, 
					     const IndexType Ntype,
					     const TypeType atom0, 
					     const TypeType atom1);
    __device__ IndexType dCalDataLength (const IndexType Ntype);
};

namespace NBForceSetting {
    __device__ void getParam (ForceIndexType forceIdx,
			      ScalorType * param, 
			      IndexType  * paramPosi,
			      ScalorType ** ptr_param);
    __device__ mdNBInteraction_t getType (ForceIndexType forceIdx,
					  const mdNBInteraction_t * typeTable);
};


namespace LennardJones6_12{
    typedef enum paramIndex {
      epsilon		= 0,
      sigma		= 1,
      shift		= 2,
      rcut		= 3
    } paramIndex_t;
    
    __host__ void initParameter (ScalorType * param, 
				 ScalorType epsilon,
				 ScalorType sigma,
				 ScalorType shift,
				 ScalorType rcut);
    __host__ ScalorType calRcut (const ScalorType * param);
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


namespace LennardJones6_12_cap {
    ScalorType lj6_12_originForce (ScalorType * param, ScalorType r);
    ScalorType lj6_12_originPoten (ScalorType * param, ScalorType r);
    ScalorType lj6_12_findCapR (ScalorType * param);

    typedef enum paramIndex {
      epsilon		= 0,
      sigma		= 1,
      shift		= 2,
      rcut		= 3,
      cap		= 4,
      capR		= 5,
      pcapR		= 6
    } paramIndex_t;

    __host__ void initParameter (ScalorType * param,
				 ScalorType epsilon,
				 ScalorType sigma,
				 ScalorType shift,
				 ScalorType rcut,
				 ScalorType cap);
    __host__ ScalorType calRcut (const ScalorType * param);
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




namespace CosTail{
    typedef enum paramIndex {
      epsilon		= 0,
      bv		= 1,
      rc		= 2,
      wc		= 3,
      wci		= 4,
      rcut		= 5
    } paramIndex_t;
    
    __host__ void initParameter (ScalorType * param, 
				 ScalorType epsilon_,
				 ScalorType bv_,
				 ScalorType wc_);
    __host__ ScalorType calRcut (const ScalorType * param);
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
nbForce (const mdNBInteraction_t ftype,
	 const ScalorType * param,
	 ScalorType diffx, ScalorType diffy, ScalorType diffz,
	 ScalorType *fx,   ScalorType *fy,   ScalorType *fz)
{
  if (ftype == mdForceLennardJones6_12){
    LennardJones6_12::force (param, diffx, diffy, diffz, fx, fy, fz);
  }
  if (ftype == mdForceLennardJones6_12_cap){
    LennardJones6_12_cap::force (param,  diffx, diffy, diffz, fx, fy, fz);
  }
  if (ftype == mdForceCosTail){
    CosTail::force (param,  diffx, diffy, diffz, fx, fy, fz);
  }
  if (ftype == mdForceNBNull){
    *fx = *fy = *fz = 0.f;
  } 
  // switch (ftype){
  // case mdForceLennardJones6_12:
  //     LennardJones6_12::force (param, diffx, diffy, diffz, fx, fy, fz);
  //     break;
  // case mdForceLennardJones6_12_cap:
  //     LennardJones6_12_cap::force (param,  diffx, diffy, diffz, fx, fy, fz);
  //     break;
  // case mdForceCosTail:
  //     CosTail::force (param,  diffx, diffy, diffz, fx, fy, fz);
  //     break;
  // default:
  //     *fx = *fy = *fz = 0.f;
  //     break;
  // }
}

__device__ void
nbForcePoten (const mdNBInteraction_t ftype,
	      const ScalorType * param,
	      ScalorType diffx, ScalorType diffy, ScalorType diffz,
	      ScalorType *fx,   ScalorType *fy,   ScalorType *fz,
	      ScalorType *dp)
{
  if (ftype == mdForceLennardJones6_12){
    *dp = LennardJones6_12::forcePoten (param, diffx, diffy, diffz, fx, fy, fz);
  }
  if (ftype == mdForceLennardJones6_12_cap){
    *dp = LennardJones6_12_cap::forcePoten (param,  diffx, diffy, diffz, fx, fy, fz);
  }
  if (ftype == mdForceCosTail){
    *dp = CosTail::forcePoten (param,  diffx, diffy, diffz, fx, fy, fz);
  }
  if (ftype == mdForceNBNull){
    *fx = *fy = *fz = 0.f;
    *dp = 0.f;
  }  
  // switch (ftype){
  // case mdForceLennardJones6_12:
  //     *dp =  LennardJones6_12::forcePoten (param, diffx, diffy, diffz, fx, fy, fz);
  //     break;
  // case mdForceLennardJones6_12_cap:
  //     *dp =  LennardJones6_12_cap::forcePoten (param,  diffx, diffy, diffz, fx, fy, fz);
  //     break;
  // case mdForceCosTail:
  //     *dp =  CosTail::forcePoten (param,  diffx, diffy, diffz, fx, fy, fz);
  //     break;
  // default:
  //     *fx = *fy = *fz = 0.f;
  //     *dp =  0.f;
  //     break;
  // }
}


inline __host__  ScalorType 
calRcut (const mdNBInteraction_t ftype,
	 const ScalorType * param)
{
  switch (ftype){
  case mdForceNBNull:
      return 0.f;
  case mdForceLennardJones6_12:
      return LennardJones6_12::calRcut (param);
  case mdForceLennardJones6_12_cap:
      return LennardJones6_12_cap::calRcut (param);
  case mdForceCosTail:
      return CosTail::calRcut (param);
  default:
      throw MDExcptUndefinedNBForceType ("calRcut");
  }
}

















__device__ ForceIndexType AtomNBForceTable::calForceIndex (ForceIndexType * table, 
							   const IndexType Ntype,
							   const TypeType atom0, 
							   const TypeType atom1)
{
  IndexType i, j;
  atom0 <= atom1 ? (i = atom0, j = atom1) : (i = atom1, j = atom0);
  return table[i * Ntype + j - ((i*(i+1)) >> 1)];
  // return table[i * Ntype + j - i];
}
__device__ ForceIndexType AtomNBForceTable::calForceIndex (volatile ForceIndexType * table, 
							   const IndexType Ntype,
							   const TypeType atom0, 
							   const TypeType atom1)
{
  IndexType i, j;
  atom0 <= atom1 ? (i = atom0, j = atom1) : (i = atom1, j = atom0);
  return table[i * Ntype + j - ((i*(i+1)) >> 1)];
  // return table[i * Ntype + j - i];
}


inline __host__ IndexType AtomNBForceTable::calDataLength (const IndexType Ntype)
{
  return ((Ntype * (Ntype + 1)) >> 1);
}

__device__ IndexType AtomNBForceTable::dCalDataLength (const IndexType Ntype)
{
  return ((Ntype * (Ntype + 1)) >> 1);
}

inline __host__ int AtomNBForceTable::deviceInitTable (const ForceIndexType * htable,
						       ForceIndexType ** dtable,
						       const IndexType Ntype)
{
  cudaMalloc ((void**)&dtable, calDataLength(Ntype) * sizeof(ForceIndexType));
  cudaMemcpy (*dtable, htable, calDataLength(Ntype) * sizeof(ForceIndexType),
	      cudaMemcpyHostToDevice);
  checkCUDAError ("AtomNBForceTable::deviceInitTable");
  return 0;
}

inline __host__ void AtomNBForceTable::setTableItem (ForceIndexType * table,
						     const IndexType Ntype,
						     const TypeType atom0, 
						     const TypeType atom1,
						     const ForceIndexType forceIdx)
{
  IndexType i, j;
  atom0 <= atom1 ? (i = atom0, j = atom1) : (i = atom1, j = atom0);
  table[i * Ntype + j - ((i*(i+1)) >> 1)] = forceIdx;
  // table[i * Ntype + j - i] = forceIdx;
}

__device__ mdNBInteraction_t NBForceSetting::getType (ForceIndexType forceIdx,
						      const mdNBInteraction_t * typeTable)
{
  return typeTable[forceIdx];
}



// __host__ int NBForceSetting::reallocHostNBForceSetting (HostNBForceSetting * hsetting, 
// 							size_t sizef, size_t sizep)
// {
//   hsetting->type = (mdNBInteraction_t *) realloc (hsetting->type, 
// 						  sizef * sizeof(mdNBInteraction_t));
//   hsetting->param = (ScalorType *) realloc (hsetting->param, 
// 					    sizep * sizeof(ScalorType));
//   hsetting->paramPosi = (IndexType *) realloc (hsetting->paramPosi, 
// 					       sizef * sizeof(IndexType));
//   if (hsetting->type == NULL || 
//       hsetting->param == NULL || 
//       hsetting->paramPosi == NULL){
//     return 1;
//   }
//   else return 0;
// }








__device__ void  NBForceSetting::getParam (ForceIndexType forceIdx,
					   ScalorType * param, 
					   IndexType  * paramPosi,
					   ScalorType ** param_ptr)
{
  *param_ptr = &param[paramPosi[(forceIdx)]];
}

inline __host__ void LennardJones6_12::initParameter (ScalorType * param, 
						      ScalorType epsilon_,
						      ScalorType sigma_,
						      ScalorType shift_,
						      ScalorType rcut_)
{
  param[epsilon] = epsilon_;
  param[sigma] = sigma_;
  param[shift] = shift_;
  param[rcut]  = rcut_;
}

inline __host__ ScalorType LennardJones6_12::calRcut (const ScalorType * param)
{
  return param[rcut];
}

__device__ void LennardJones6_12::force (const ScalorType * param,
					 ScalorType diffx,
					 ScalorType diffy,
					 ScalorType diffz,
					 ScalorType *fx, 
					 ScalorType *fy,
					 ScalorType *fz)
{
  ScalorType dr2 = diffx*diffx + diffy*diffy + diffz*diffz;
  ScalorType boolScalor;
  if (dr2 > param[rcut]*param[rcut]) boolScalor = 0.f;
  else boolScalor = - 24.f;
  ScalorType ri2 = __frcp_rn(dr2);
  ScalorType sri2 = param[sigma] * param[sigma] * ri2;
  ScalorType sri6 = sri2*sri2*sri2;
  ScalorType scalor = boolScalor * param[epsilon] * (2.f * (sri6*sri6) - sri6) * ri2;
  *fx = diffx * scalor;
  *fy = diffy * scalor;
  *fz = diffz * scalor;

  // ScalorType dr2 = diffx*diffx + diffy*diffy + diffz*diffz;
  // ScalorType boolScalor;
  // if (dr2 > param[rcut]*param[rcut]) boolScalor = 0.f;
  // else boolScalor = 24.f;
  // Scalortype ri2 = __frcp_rn(dr2);
  // ScalorType sri2 = param[sigma] * param[sigma] * ri2;
  // ScalorType sri4 = sri2*sri2;
  // ScalorType sri8 = sri4*sri4;
  // ScalorType scalor = boolScalor * param[epsilon] * (2.f * (sri2*sri4*sri8) - sri8);
  // *fx = diffx * scalor;
  // *fy = diffy * scalor;
  // *fz = diffz * scalor;
}
  
__device__ ScalorType LennardJones6_12::forcePoten (const ScalorType * param,
						    ScalorType diffx,
						    ScalorType diffy,
						    ScalorType diffz,
						    ScalorType *fx, 
						    ScalorType *fy,
						    ScalorType *fz)
{
  ScalorType dr2 = diffx*diffx + diffy*diffy + diffz*diffz;
  ScalorType boolScalor;
  if (dr2 > param[rcut]*param[rcut]) boolScalor = 0.f;
  else boolScalor = 4.f;
  ScalorType ri2 = 1.f/dr2;
  ScalorType sri2 = param[sigma] * param[sigma] * ri2;
  ScalorType sri6 = sri2*sri2*sri2;
  ScalorType scalor = - boolScalor * param[epsilon] * (12.f * sri6*sri6 - 6.f * sri6) * ri2;
  *fx = diffx * scalor;
  *fy = diffy * scalor;
  *fz = diffz * scalor;
  return boolScalor * param[epsilon] * (sri6*sri6 - sri6 + param[shift]);
}
  


inline __host__ ScalorType LennardJones6_12_cap::
lj6_12_originForce (ScalorType * param, ScalorType r)
{
  ScalorType dr2 = r*r;
  if (dr2 > param[rcut]*param[rcut]) {
    return 0.f;
  }
  ScalorType ri = 1.f/r;
  ScalorType ri2 = ri*ri;
  ScalorType sri2 = param[sigma]*param[sigma] * ri2;
  ScalorType sri6 = sri2*sri2*sri2;
  return 4.f * param[epsilon] * (12.f * sri6*sri6 - 6.f * sri6) * ri;
}

inline __host__ ScalorType LennardJones6_12_cap::
lj6_12_originPoten (ScalorType * param, ScalorType r)
{
  ScalorType dr2 = r*r;
  if (dr2 > param[rcut]*param[rcut]) {
    return 0.f;
  }
  ScalorType ri = 1.f/r;
  ScalorType ri2 = ri*ri;
  ScalorType sri2 = param[sigma]*param[sigma] * ri2;
  ScalorType sri6 = sri2*sri2*sri2;
  return 4.f * param[epsilon] * (sri6*sri6 - sri6 + param[shift]);
}

inline __host__ ScalorType LennardJones6_12_cap::
lj6_12_findCapR (ScalorType * param)
{
  ScalorType point = 0.95f * param[rcut];
  while (lj6_12_originForce(param, point) < param[cap]){
    point *= 0.95f;
  }
  ScalorType lower = point;
  ScalorType upper = param[rcut];
  while (fabs(upper - lower) > 1e-4){
    ScalorType mid = 0.5f * (upper + lower);
    if (lj6_12_originForce(param, mid) > param[cap]) {
      lower = mid;
    }
    else {
      upper = mid;
    }
  }
  return 0.5f * (upper + lower);
}


inline __host__ void LennardJones6_12_cap::initParameter (ScalorType * param,
							  ScalorType epsilon_,
							  ScalorType sigma_,
							  ScalorType shift_,
							  ScalorType rcut_,
							  ScalorType cap_)
{
  param[epsilon] = epsilon_;
  param[sigma] = sigma_;
  param[shift] = shift_;
  param[rcut] = rcut_;
  param[cap] = cap_;
  param[capR] = lj6_12_findCapR (param);
  param[pcapR] = lj6_12_originPoten(param, param[capR]);
}

inline __host__ ScalorType LennardJones6_12_cap::calRcut (const ScalorType * param)
{
  return param[rcut];
}

__device__ void LennardJones6_12_cap::force (const ScalorType * param,
					     ScalorType diffx,
					     ScalorType diffy,
					     ScalorType diffz,
					     ScalorType *fx, 
					     ScalorType *fy, 
					     ScalorType *fz)
{
  ScalorType dr2 = diffx*diffx + diffy*diffy + diffz*diffz;
  if (dr2 > param[rcut]*param[rcut]) {
    *fx = *fy = *fz = 0.f;
    return;
  }
  if (dr2 == 0.f){
    *fx = - param[cap];
    *fy = 0.f;
    *fz = 0.f;
    return;
  }
  if (dr2 < param[capR]*param[capR]){
    ScalorType s = 1.f / sqrtf(diffx*diffx + diffy*diffy + diffz*diffz);
    *fx = -s * param[cap] * diffx;
    *fy = -s * param[cap] * diffy;
    *fz = -s * param[cap] * diffz;
    return;
  }
  ScalorType ri2 = 1.f/dr2;
  ScalorType sri2 = param[sigma]*param[sigma] * ri2;
  ScalorType sri6 = sri2*sri2*sri2;
  ScalorType scalor = 4.f * param[epsilon] * (12.f * sri6*sri6 - 6.f * sri6) * ri2;
  *fx = -diffx * scalor;
  *fy = -diffy * scalor;
  *fz = -diffz * scalor;
}


__device__ ScalorType LennardJones6_12_cap::forcePoten (const ScalorType * param,
							ScalorType diffx,
							ScalorType diffy,
							ScalorType diffz,
							ScalorType *fx, 
							ScalorType *fy, 
							ScalorType *fz)
{
  ScalorType dr2 = diffx*diffx + diffy*diffy + diffz*diffz;
  if (dr2 > param[rcut]*param[rcut]) {
    *fx = *fy = *fz = 0.f;
    return 0.f;
  }
  if (dr2 == 0.f){
    *fx = -param[cap];
    *fy = 0.f;
    *fz = 0.f;
    return 0.f;
  }
  if (dr2 < param[capR]*param[capR]){
    ScalorType s = 1.f / sqrtf(diffx*diffx + diffy*diffy + diffz*diffz);
    *fx = -s * param[cap] * diffx;
    *fy = -s * param[cap] * diffy;
    *fz = -s * param[cap] * diffz;
    return -param[cap] * (sqrtf(dr2) - param[capR]) + param[pcapR] ;
  }
  ScalorType ri2 = 1.f/dr2;
  ScalorType sri2 = param[sigma]*param[sigma] * ri2;
  ScalorType sri6 = sri2*sri2*sri2;
  ScalorType scalor = 4.f * param[epsilon] * (12.f * sri6*sri6 - 6.f * sri6) * ri2;
  *fx = -diffx * scalor;
  *fy = -diffy * scalor;
  *fz = -diffz * scalor;
  return 4.f * param[epsilon] * (sri6*sri6 - sri6 + param[shift]);
}





inline __host__ void CosTail::initParameter (ScalorType * param, 
					     ScalorType epsilon_,
					     ScalorType bv_,
					     ScalorType wc_)
{
  param[bv] = bv_;
  param[epsilon] = epsilon_;
  param[wc] = wc_;
  
  if (param[wc] != 0.f)   param[wci] = 1.f / param[wc];
  else param[wci] = 0.f;
  
  param[rc] = param[bv] * powf (2.f, 1.f / 6.f);
  param[rcut] = param[rc] + param[wc];
}

inline __host__ ScalorType CosTail::calRcut (const ScalorType * param)
{
  return param[rcut];
}  
  
__device__ void CosTail::force (const ScalorType * param,
				ScalorType diffx,
				ScalorType diffy,
				ScalorType diffz,
				ScalorType *fx, 
				ScalorType *fy,
				ScalorType *fz)
{
  ScalorType dr2 = diffx*diffx + diffy*diffy + diffz*diffz;
  ScalorType fscalor = 0.f;
  ScalorType ri = 1.f / sqrtf (dr2);

  if (dr2 <= param[rc] * param[rc]){
    ScalorType sri = param[bv] * ri;
    ScalorType sri2 = sri * sri;
    ScalorType sri6 = sri2*sri2*sri2;
    fscalor = - 24.f * param[epsilon] * (2.f * sri6*sri6 - sri6) * ri * ri;
  }
  else if (dr2 < param[rcut] * param[rcut]) {
    ScalorType tmp = M_PIF * param[wci];
    ScalorType term = 0.5f *  (dr2 * ri - param[rc]) * tmp;
    fscalor = param[epsilon] * tmp *
    	__cosf(term) * __sinf(term) * ri;
  }

  *fx = diffx * fscalor;
  *fy = diffy * fscalor;
  *fz = diffz * fscalor;
}
  
__device__ ScalorType CosTail::forcePoten (const ScalorType * param,
					   ScalorType diffx,
					   ScalorType diffy,
					   ScalorType diffz,
					   ScalorType *fx, 
					   ScalorType *fy,
					   ScalorType *fz)
{
  ScalorType dr2 = diffx*diffx + diffy*diffy + diffz*diffz;
  ScalorType rvalue = 0.f;
  ScalorType fscalor = 0.f;
  ScalorType ri = 1.f / sqrtf (dr2);

  if (dr2 <= param[rc] * param[rc]){
    ScalorType sri = param[bv] * ri;
    ScalorType sri2 = sri * sri;
    ScalorType sri6 = sri2*sri2*sri2;
    fscalor = - 24.f * param[epsilon] * (2.f * sri6*sri6 - sri6) * ri * ri;
    rvalue = 4.f * param[epsilon] * (sri6*sri6 - sri6);
    if (param[wc] == 0.f) rvalue += param[epsilon];
    else rvalue += 0.f;
  }
  else if (dr2 < param[rcut] * param[rcut]) {
    ScalorType term = 0.5f * M_PIF * (dr2 * ri - param[rc]) * param[wci];
    ScalorType cost = __cosf(term);
    fscalor = param[epsilon] * M_PIF * param[wci] *
	cost * __sinf(term) * ri;
    rvalue = - param[epsilon] * cost * cost;
  }

  *fx = diffx * fscalor;
  *fy = diffy * fscalor;
  *fz = diffz * fscalor;

  return rvalue;
}

#endif
