#ifndef __BoxGeometry_h_wanghan__
#define __BoxGeometry_h_wanghan__

#include "common.h"
#include "MDSystem.h"

namespace RectangularBoxGeometry{
    struct RectangularBox 
    {
      VectorType size;
      VectorType sizei;
    };

    enum BoxDirection {
      mdRectBoxDirectionX	= 1,
      mdRectBoxDirectionY	= 2,
      mdRectBoxDirectionZ	= 4
    } ;
    typedef int BoxDirection_t;
    
    __host__ void setBoxSize (ScalorType x, ScalorType y, ScalorType z,
			      RectangularBox *box);
    __device__ void moveParticleToBox (
	RectangularBox box,
	ScalorType * x, ScalorType * y, ScalorType * z,
	IntScalorType * noix, 
	IntScalorType * noiy, 
	IntScalorType * noiz);
    __device__ void moveParticleToBox_1image(
	ScalorType boxx,
	ScalorType *x,  
	IntScalorType *noix);
    __device__ void shortestImage (RectangularBox box,
				   ScalorType * x, ScalorType * y, ScalorType * z);
    __device__ void shortestImage (const ScalorType boxL, const ScalorType boxLi,
				   ScalorType * x);
    // __device__ void normalizeSystem (RectangularBox box, DeviceMDData * ddata);
    // __device__ void normalizeSystem (RectangularBox box, 
    // 				     IndexType numAtom,
    // 				     ScalorType *coordx, 
    // 				     ScalorType *coordy, 
    // 				     ScalorType *coordz,
    // 				     IntScalorType *coordNoix, 
    // 				     IntScalorType *coordNoiy, 
    // 				     IntScalorType *coordNoiz);
    
    // needs ceil(numAtom/blockDim.x) blocks
    __global__ void normalizeSystem (RectangularBox box, 
				     IndexType numAtom,
#ifndef COORD_IN_ONE_VEC
				     ScalorType *coordx, 
				     ScalorType *coordy, 
				     ScalorType *coordz,
#else
				     CoordType * coord,
#endif
				     IntScalorType *coordNoix, 
				     IntScalorType *coordNoiy, 
				     IntScalorType *coordNoiz);
}
;

inline __host__ void RectangularBoxGeometry::setBoxSize (
    ScalorType x, ScalorType y, ScalorType z,
    RectangularBoxGeometry::RectangularBox *box)
{
  box->size.x = x;
  box->size.y = y;
  box->size.z = z;
  box->sizei.x = 1/x;
  box->sizei.y = 1/y;
  box->sizei.z = 1/z;
}

__device__ void RectangularBoxGeometry::moveParticleToBox(
    RectangularBoxGeometry::RectangularBox rectBox,
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

__device__ void RectangularBoxGeometry::moveParticleToBox_1image(
    ScalorType boxx,
    ScalorType *x,  
    IntScalorType *noix)
{
  IntScalorType tmp;
  if (*x >= boxx){
    tmp = 1;
    *noix += tmp;
    *x -= tmp * boxx;
  }
  else if (*x < 0) {
    tmp = -1;
    *noix += tmp;
    *x -= tmp * boxx;
  }
}


  

__device__ void RectangularBoxGeometry::shortestImage (
    RectangularBoxGeometry::RectangularBox box,
    ScalorType * x, ScalorType * y, ScalorType * z)
{
  *x -= floorf(*x * box.sizei.x + 0.5f) * box.size.x;
  *y -= floorf(*y * box.sizei.y + 0.5f) * box.size.y;
  *z -= floorf(*z * box.sizei.z + 0.5f) * box.size.z;
  // if (*x >  0.5f * box.size.x) *x -= box.size.x;
  // if (*x < -0.5f * box.size.x) *x += box.size.x;
  // if (*y >  0.5f * box.size.y) *y -= box.size.y;
  // if (*y < -0.5f * box.size.y) *y += box.size.y;
  // if (*z >  0.5f * box.size.z) *z -= box.size.z;
  // if (*z < -0.5f * box.size.z) *z += box.size.z;
}

__device__ void RectangularBoxGeometry::shortestImage (
    const ScalorType boxL, const ScalorType boxLi,
    ScalorType * x)
{
  *x -= floorf(*x * boxLi + 0.5f) * boxL;
}


// __device__ void RectangularBoxGeometry::normalizeSystem (
//     RectangularBoxGeometry::RectangularBox box,
//     DeviceMDData * ddata)
// {
//   IndexType bid = blockIdx.x;
//   IndexType ii = bid * blockDim.x + threadIdx.x;
  
//   if (ii < ddata->numAtom)
//     RectangularBoxGeometry::moveParticleToBox (
// 	box, 
// 	&ddata->coordx[ii], &ddata->coordy[ii], &ddata->coordz[ii], 
// 	&ddata->coordNoix[ii], &ddata->coordNoiy[ii], &ddata->coordNoiz[ii]);
// }

// __device__ void normalizeSystem (RectangularBoxGeometry::RectangularBox box, 
// 				 IndexType numAtom,
// 				 ScalorType *coordx, 
// 				 ScalorType *coordy, 
// 				 ScalorType *coordz,
// 				 IntScalorType *coordNoix, 
// 				 IntScalorType *coordNoiy, 
// 				 IntScalorType *coordNoiz)
// {
//   IndexType bid = blockIdx.x;
//   IndexType ii = bid * blockDim.x + threadIdx.x;
  
//   if (ii < numAtom)
//     RectangularBoxGeometry::moveParticleToBox (
// 	box, 
// 	&coordx[ii], &coordy[ii], &coordz[ii], 
// 	&coordNoix[ii], &coordNoiy[ii], &coordNoiz[ii]);
// }






//  *noix -= int (x * 
  // int tmp;
  // if (p.r()[0] >= 0){
  //   tmp = int(p.r()[0] * li[0] + .5);
  //   p.noi()[0] -= tmp;
  // }
  // else{
  //   tmp = int(p.r()[0] * li[0] - .5);
  //   p.noi()[0] -= tmp;
  // }
  // p.r()[0] -= tmp * boxl[0];

  // if (p.r()[1] >= 0){
  //   tmp = int(p.r()[1] * li[1] + .5);
  //   p.noi()[1] -= tmp;
  // }
  // else{
  //   tmp = int(p.r()[1] * li[1] - .5);
  //   p.noi()[1] -= tmp;
  // }
  // p.r()[1] -= tmp * boxl[1];

  // if (p.r()[2] >= 0){
  //   tmp = int(p.r()[2] * li[2] + .5);
  //   p.noi()[2] -= tmp;
  // }
  // else{
  //   tmp = int(p.r()[2] * li[2] - .5);
  //   p.noi()[2] -= tmp;
  // }
  // p.r()[2] -= tmp * boxl[2];
// }


#endif 
