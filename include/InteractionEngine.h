#ifndef __InteractionEngine_h_wanghan__
#define __InteractionEngine_h_wanghan__

#include "common.h"
#include "BoxGeometry.h"
#include "MDSystem.h"
#include "NeighborList.h"
#include "Statistic.h"
#include "Auxiliary.h"
#include "BondList.h"
#include "AngleList.h"
#include "ExclusionList.h"

using namespace RectangularBoxGeometry;

// needs ceil(numAtom/blockDim.x) blocks
__global__ void
calNonBondedInteraction_all (const IndexType		numAtom,
			     const CoordType *		coord,
			     ScalorType *		forcx,
			     ScalorType *		forcy, 
			     ScalorType *		forcz,
			     const TypeType *		type,
			     const RectangularBox	box,
			     const ScalorType		rcut,
			     mdError_t *		ptr_de);
__global__ void
calNonBondedInteraction_all  (const IndexType		numAtom,
			      const CoordType *		coord,
			      ScalorType *		forcx,
			      ScalorType *		forcy, 
			      ScalorType *		forcz,
			      const TypeType *		type,
			      const RectangularBox	box,
			      const ScalorType		rcut,
			      ScalorType *		statistic_nb_buff0,
			      ScalorType *		statistic_nb_buff1,
			      ScalorType *		statistic_nb_buff2,
			      ScalorType *		statistic_nb_buff3,
			      mdError_t *		ptr_de);

// needs ceil(numAtom/blockDim.x) blocks
__global__ void
calNonBondedInteraction_neighbor (const IndexType		numAtom,
				  const CoordType *		coord,
				  ScalorType *			forcx,
				  ScalorType *			forcy, 
				  ScalorType *			forcz,
				  const TypeType *		type,
				  const RectangularBox		box,
				  const DeviceNeighborList	nlist);
__global__ void
calNonBondedInteraction_neighbor (const IndexType		numAtom,
				  const CoordType *		coord,
				  ScalorType *			forcx,
				  ScalorType *			forcy, 
				  ScalorType *			forcz,
				  const TypeType *		type,
				  const RectangularBox		box,
				  const DeviceNeighborList	nlist,
				  ScalorType *			statistic_buff0,
				  ScalorType *			statistic_buff1,
				  ScalorType *			statistic_buff2,
				  ScalorType *			statistic_buff3);
__global__ void
calNonBondedInteraction_neighbor (const IndexType		numAtom,
				  const CoordType *		coord,
				  ScalorType *			forcx,
				  ScalorType *			forcy, 
				  ScalorType *			forcz,
				  const TypeType *		type,
				  const RectangularBox		box,
				  const DeviceNeighborList	nlist,
				  const DeviceExclusionList	dexcllist,
				  const bool			sharedExclusionList);

__global__ void
calNonBondedInteraction_neighbor (const IndexType		numAtom,
				  const CoordType *		coord,
				  ScalorType *			forcx,
				  ScalorType *			forcy, 
				  ScalorType *			forcz,
				  const TypeType *		type,
				  const RectangularBox		box,
				  const DeviceNeighborList	nlist,
				  const DeviceExclusionList	dexcllist,
				  const bool			sharedExclusionList,
				  ScalorType *			statistic_buff0,
				  ScalorType *			statistic_buff1,
				  ScalorType *			statistic_buff2,
				  ScalorType *			statistic_buff3);

// needs NCell blocks
__global__ void
calNonBondedInteraction_cell (const IndexType		numAtom,
			      const CoordType *		coord,
			      ScalorType *		forcx,
			      ScalorType *		forcy, 
			      ScalorType *		forcz,
			      const TypeType *		type,
			      const RectangularBox	box,
			      const DeviceCellList	clist,
			      const ScalorType		rcut,
			      mdError_t *		ptr_de);
__global__ void
calNonBondedInteraction_cell (const IndexType		numAtom,
			      const CoordType *		coord,
			      ScalorType *		forcx,
			      ScalorType *		forcy, 
			      ScalorType *		forcz,
			      const TypeType *		type,
			      const RectangularBox	box,
			      const DeviceCellList	clist,
			      const ScalorType		rcut,
			      ScalorType *		statistic_nb_buff0,
			      ScalorType *		statistic_nb_buff1,
			      ScalorType *		statistic_nb_buff2,
			      ScalorType *		statistic_nb_buff3,
			      mdError_t *		ptr_de);


// __global__ void
// calNonBondedInteraction (const IndexType		numAtom,
// 			 const CoordType *		coord,
// 			 ScalorType *			forcx,
// 			 ScalorType *			forcy, 
// 			 ScalorType *			forcz,
// 			 const TypeType *		type,
// 			 const RectangularBox		box,
// 			 DeviceCellList			clist,
// 			 const ScalorType		rcut,
// 			 DeviceNeighborList		nlist,
// 			 mdError_t *			ptr_de);
// __global__ void
// calNonBondedInteraction (const IndexType		numAtom,
// 			 const CoordType *		coord,
// 			 ScalorType *			forcx,
// 			 ScalorType *			forcy, 
// 			 ScalorType *			forcz,
// 			 const TypeType *		type,
// 			 const RectangularBox		box,
// 			 DeviceCellList			clist,
// 			 const ScalorType		rcut,
// 			 DeviceNeighborList		nlist,
// 			 ScalorType *			statistic_nb_buff0,
// 			 ScalorType *			statistic_nb_buff1,
// 			 ScalorType *			statistic_nb_buff2,
// 			 ScalorType *			statistic_nb_buff3,
// 			 mdError_t *			ptr_de);


// needs ceil(numAtom/blockDim.x) blocks
__global__ void
calTwinRangeCorrection_all    (const IndexType		numAtom,
			       const CoordType *	coord,
			       ScalorType *		forcx,
			       ScalorType *		forcy, 
			       ScalorType *		forcz,
			       const TypeType *		type,
			       const RectangularBox	box,
			       const ScalorType		rcut1,
			       const ScalorType		rcut2,
			       ScalorType *		statistic_nb_buff0,
			       ScalorType *		statistic_nb_buff1,
			       ScalorType *		statistic_nb_buff2,
			       ScalorType *		statistic_nb_buff3,
			       mdError_t *		ptr_de);
// needs NCell blocks
__global__ void
calTwinRangeCorrection_cell (const IndexType		numAtom,
			     const CoordType *		coord,
			     ScalorType *		forcx,
			     ScalorType *		forcy, 
			     ScalorType *		forcz,
			     const TypeType *		type,
			     const RectangularBox	box,
			     const DeviceCellList	clist,
			     const ScalorType		rcut1,
			     const ScalorType		rcut2,
			     ScalorType *		statistic_nb_buff0,
			     ScalorType *		statistic_nb_buff1,
			     ScalorType *		statistic_nb_buff2,
			     ScalorType *		statistic_nb_buff3,
			     mdError_t *		ptr_de);

// needs ceil(numAtom/blockDim.x) blocks
__global__ void
buildNeighborListCalTwinRangeCorr_all (const IndexType		numAtom,
				       const CoordType *	coord,
				       ScalorType *		forcx,
				       ScalorType *		forcy, 
				       ScalorType *		forcz,
				       const TypeType *		type,
				       const RectangularBox	box,
				       const ScalorType		rcut1,
				       const ScalorType		rcut2,
				       DeviceNeighborList	nlist,
				       ScalorType *		statistic_nb_buff0,
				       ScalorType *		statistic_nb_buff1,
				       ScalorType *		statistic_nb_buff2,
				       ScalorType *		statistic_nb_buff3,
				       mdError_t *		ptr_de);
// needs NCell blocks
__global__ void
buildNeighborListCalTwinRangeCorr_cell (const IndexType		numAtom,
					const CoordType *	coord,
					ScalorType *		forcx,
					ScalorType *		forcy, 
					ScalorType *		forcz,
					const TypeType *	type,
					const RectangularBox	box,
					const DeviceCellList	clist,
					const ScalorType	rcut1,
					const ScalorType	rcut2,
					DeviceNeighborList	nlist,
					ScalorType *		statistic_nb_buff0,
					ScalorType *		statistic_nb_buff1,
					ScalorType *		statistic_nb_buff2,
					ScalorType *		statistic_nb_buff3,
					mdError_t *		ptr_de);



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
__global__ void calAngleInteraction (const IndexType numAtom,
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
				     const DeviceAngleList bdlist);
__global__ void calAngleInteraction (const IndexType numAtom,
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
				     const DeviceAngleList bdlist,
				     ScalorType * statistic_buff,
				     ScalorType * statistic_b_buff1,
				     ScalorType * statistic_b_buff2,
				     ScalorType * statistic_b_buff3,
				     mdError_t * ptr_de);

// __global__ void calNonAngleedInteraction (const IndexType numAtom,
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


__global__ void calNonBondedInteraction (
    const IndexType numAtom,
    const CoordType * coord,
    ScalorType * forcx,
    ScalorType * forcy, 
    ScalorType * forcz,
    const TypeType * type,
    const RectangularBox box,
    DeviceCellList clist,
    mdError_t * ptr_de);
__global__ void calNonBondedInteraction (
    const IndexType numAtom,
    const CoordType * coord,
    ScalorType * forcx,
    ScalorType * forcy, 
    ScalorType * forcz,
    const TypeType * type,
    const RectangularBox box,
    DeviceCellList clist,
    ScalorType * statistic_nb_buff0,
    ScalorType * statistic_nb_buff1,
    ScalorType * statistic_nb_buff2,
    ScalorType * statistic_nb_buff3,
    mdError_t * ptr_de);


__global__ void
widomDeltaPoten_NVT (const IndexType		numTestParticle,
		     const CoordType *		coordTestParticle,
		     const TypeType *		typeTestParticle,
		     const IndexType		numAtom,
		     const CoordType *		coord,
		     const TypeType *		type,
		     const RectangularBox	box,
		     DeviceCellList		clist,
		     ScalorType *		statistic_nb_buff0,
		     mdError_t *		ptr_de);
__global__ void
widomDeltaPoten_allPair_NVT (const IndexType		numTestParticle,
			     const CoordType *		coordTestParticle,
			     const TypeType *		typeTestParticle,
			     const IndexType		numAtom,
			     const CoordType *		coord,
			     const TypeType *		type,
			     const RectangularBox	box,
			     const ScalorType		rlist,
			     ScalorType *		statistic_nb_buff0,
			     mdError_t *		ptr_de);


#endif
