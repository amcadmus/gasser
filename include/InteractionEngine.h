#ifndef __InteractionEngine_h_wanghan__
#define __InteractionEngine_h_wanghan__

#include "common.h"
#include "BoxGeometry.h"
#include "MDSystem.h"
#include "NeighborList.h"
#include "Statistic.h"
#include "Auxiliary.h"
#include "BondList.h"

using namespace RectangularBoxGeometry;

// needs ceil(numAtom/blockDim.x) blocks
__global__ void calNonBondedInteraction (const IndexType numAtom,
#ifndef COORD_IN_ONE_VEC
					 const ScalorType * coordx,
					 const ScalorType * coordy, 
					 const ScalorType * coordz,
#else
					 const CoordType * coord,
#endif
					 ScalorType * forcx,
					 ScalorType * forcy, 
					 ScalorType * forcz,
					 const TypeType * type,
					 const RectangularBox box,
					 const DeviceNeighborList nlist,
					 mdError_t * ptr_de,
					 IndexType * errorIndex,
					 ScalorType * errorScalor);
__global__ void calNonBondedInteraction (const IndexType numAtom,
#ifndef COORD_IN_ONE_VEC
					 const ScalorType * coordx,
					 const ScalorType * coordy, 
					 const ScalorType * coordz,
#else
					 const CoordType * coord,
#endif
					 ScalorType * forcx,
					 ScalorType * forcy, 
					 ScalorType * forcz,
					 const TypeType * type,
					 const RectangularBox box,
					 const DeviceNeighborList nlist,
					 ScalorType * statistic_buff0,
					 ScalorType * statistic_buff1,
					 ScalorType * statistic_buff2,
					 ScalorType * statistic_buff3,
					 mdError_t * ptr_de,
					 IndexType * errorIndex,
					 ScalorType * errorScalor);
__global__ void calBondInteraction (const IndexType numAtom,
#ifndef COORD_IN_ONE_VEC
				    const ScalorType * coordx,
				    const ScalorType * coordy, 
				    const ScalorType * coordz,
#else
				    const CoordType * coord,
#endif
				    ScalorType * forcx,
				    ScalorType * forcy, 
				    ScalorType * forcz,
				    const RectangularBox box,
				    const DeviceBondList bdlist);
__global__ void calBondInteraction (const IndexType numAtom,
#ifndef COORD_IN_ONE_VEC
				    const ScalorType * coordx,
				    const ScalorType * coordy, 
				    const ScalorType * coordz,
#else
				    const CoordType * coord,
#endif
				    ScalorType * forcx,
				    ScalorType * forcy, 
				    ScalorType * forcz,
				    const RectangularBox box,
				    const DeviceBondList bdlist,
				    ScalorType * statistic_buff0,
				    ScalorType * statistic_buff1,
				    ScalorType * statistic_buff2,
				    ScalorType * statistic_buff3,
				    mdError_t * ptr_de);

// __global__ void calNonBondedInteraction (const IndexType numAtom,
// 					 const ScalorType * coordx,
// 					 const ScalorType * coordy, 
// 					 const ScalorType * coordz,
// 					 ScalorType * forcx,
// 					 ScalorType * forcy, 
// 					 ScalorType * forcz,
// 					 const TypeType * type,
// 					 const RectangularBox box,
// 					 const DeviceNeighborList nlist,
// 					 ScalorType * nbP,
// 					 ScalorType * statistic_buff);

// needs ceil(numAtom/blockDim.x) blocks
__global__ void clearForce (const IndexType numAtom,
			    ScalorType * forcx,
			    ScalorType * forcy, 
			    ScalorType * forcz);

  // IndexType jj = 0;
  // while (__all(jj < myNumNei)){
  //   ++jj;
  //   __syncthreads();
  // }    
  // ScalorType targetx, targety, targetz;





#endif
