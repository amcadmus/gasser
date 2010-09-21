#define DEVICE_CODE
#include "BoxGeometry.h"


#ifndef COORD_IN_ONE_VEC
__global__ void RectangularBoxGeometry::normalizeSystem (
    RectangularBoxGeometry::RectangularBox box, 
    IndexType numAtom,
    ScalorType *coordx, 
    ScalorType *coordy, 
    ScalorType *coordz,
    IntScalorType *coordNoix, 
    IntScalorType *coordNoiy, 
    IntScalorType *coordNoiz)
{
  IndexType bid = blockIdx.x;
  IndexType ii = bid * blockDim.x + threadIdx.x;
  
  if (ii < numAtom){
    RectangularBoxGeometry::moveParticleToBox_1image (
	box.size.x, &coordx[ii], &coordNoix[ii]);
    RectangularBoxGeometry::moveParticleToBox_1image (
	box.size.y, &coordy[ii], &coordNoiy[ii]);
    RectangularBoxGeometry::moveParticleToBox_1image (
	box.size.z, &coordz[ii], &coordNoiz[ii]);
  }
}
#else
__global__ void RectangularBoxGeometry::normalizeSystem (
    RectangularBoxGeometry::RectangularBox box, 
    IndexType numAtom,
    CoordType * coord,
    IntScalorType *coordNoix, 
    IntScalorType *coordNoiy, 
    IntScalorType *coordNoiz)
{
  IndexType bid = blockIdx.x;
  IndexType ii = bid * blockDim.x + threadIdx.x;
  
  if (ii < numAtom){
    RectangularBoxGeometry::moveParticleToBox_1image (
	box.size.x, &(coord[ii].x), &coordNoix[ii]);
    RectangularBoxGeometry::moveParticleToBox_1image (
	box.size.y, &(coord[ii].y), &coordNoiy[ii]);
    RectangularBoxGeometry::moveParticleToBox_1image (
	box.size.z, &(coord[ii].z), &coordNoiz[ii]);
  }
}
#endif
