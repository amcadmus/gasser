#ifndef __BoxGeometry_h_wanghan__
#define __BoxGeometry_h_wanghan__

#include "common.h"
#include <math.h> 

// #include "MDSystem.h"

namespace RectangularBoxGeometry{
  struct RectangularBox 
  {
    HostVectorType size;
    HostVectorType sizei;
  };
  
  enum BoxDirection {
    mdRectBoxDirectionX	= 1,
    mdRectBoxDirectionY	= 2,
    mdRectBoxDirectionZ	= 4
  } ;
  typedef int BoxDirection_t;
  
  void setBoxSize (ScalorType x,
		   ScalorType y,
		   ScalorType z,
		   RectangularBox * box);
  void setBoxSize (const ScalorType & x,
		   const ScalorType & y,
		   const ScalorType & z,
		   RectangularBox & box);
  void hostMoveParticleToBox (RectangularBox box,
			      ScalorType * x,
			      ScalorType * y,
			      ScalorType * z,
			      IntScalorType * noix, 
			      IntScalorType * noiy, 
			      IntScalorType * noiz);
  

#ifdef DEVICE_CODE
  static __device__ void
  moveParticleToBox (RectangularBox box,
		     ScalorType * x,
		     ScalorType * y,
		     ScalorType * z,
		     IntScalorType * noix, 
		     IntScalorType * noiy, 
		     IntScalorType * noiz);
  static __device__ void
  moveParticleToBox_1image (ScalorType boxx,
			    ScalorType *x,  
			    IntScalorType *noix);
  static __device__ void
  shortestImage (RectangularBox box,
		 ScalorType * x,
		 ScalorType * y,
		 ScalorType * z);
  static __device__ void
  shortestImage (const ScalorType boxL,
		 const ScalorType boxLi,
		 ScalorType * x);
  static __device__ void
  shortestImage (const ScalorType boxL,
		 ScalorType & x);
  // needs ceil(numAtom/blockDim.x) blocks
  __global__ void
  normalizeSystem (RectangularBox box, 
		   IndexType numAtom,
		   CoordType * coord,
		   IntScalorType *coordNoix, 
		   IntScalorType *coordNoiy, 
		   IntScalorType *coordNoiz);
#endif
}
;

inline void RectangularBoxGeometry::
hostMoveParticleToBox (RectangularBox rectBox,
		       ScalorType * x,
		       ScalorType * y,
		       ScalorType * z,
		       IntScalorType * noix, 
		       IntScalorType * noiy, 
		       IntScalorType * noiz)
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

inline void RectangularBoxGeometry::
setBoxSize (ScalorType x,
	    ScalorType y,
	    ScalorType z,
	    RectangularBox *box)
{
  box->size.x = x;
  box->size.y = y;
  box->size.z = z;
  box->sizei.x = 1/x;
  box->sizei.y = 1/y;
  box->sizei.z = 1/z;
}

inline void RectangularBoxGeometry::
setBoxSize (const ScalorType & x,
	    const ScalorType & y,
	    const ScalorType & z,
	    RectangularBox & box)
{
  box.size.x = x;
  box.size.y = y;
  box.size.z = z;
  box.sizei.x = 1/x;
  box.sizei.y = 1/y;
  box.sizei.z = 1/z;
}


#ifdef DEVICE_CODE
__device__ void RectangularBoxGeometry::
moveParticleToBox(RectangularBox rectBox,
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

__device__ void RectangularBoxGeometry::
moveParticleToBox_1image(ScalorType boxx,
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

__forceinline__ __device__ void RectangularBoxGeometry::
shortestImage (RectangularBox box,
	       ScalorType * x,
	       ScalorType * y,
	       ScalorType * z)
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

__device__ void RectangularBoxGeometry::
shortestImage (const ScalorType boxL,
	       const ScalorType boxLi,
	       ScalorType * x)
{
  *x -= floorf(*x * boxLi + 0.5f) * boxL;
  // if (*x >=  0.5f * boxL) *x -= boxL;
  // if (*x < -0.5f * boxL) *x += boxL;
}

__device__ void RectangularBoxGeometry::
shortestImage (const ScalorType boxL,
	       ScalorType & x)
{
  ScalorType tmp = 0.5f * boxL;
  if (x >= tmp) x -= boxL;
  if (x < -tmp) x += boxL;
}


#endif // DEVICE_CODE

#endif 
