#ifndef __Integrator_h_wanghan__
#define __Integrator_h_wanghan__

#include "common.h"
#include "MDSystem.h"
#include "Statistic.h"

// needs 1 thread
__global__ void initRemoveTranslationalFreedom ();

// needs ceil(numAtom/blockDim.x) blocks
__global__ void prepareRemoveTranslationalFreedom (IndexType numAtom,
						   ScalorType * mass,
						   ScalorType * velox,
						   ScalorType * veloy,
						   ScalorType * veloz,
						   ScalorType * buffx,
						   ScalorType * buffy,
						   ScalorType * buffz);
// needs ceil(numAtom/blockDim.x) blocks
__global__ void removeFreedom (IndexType numAtom,
			       ScalorType * velox, 
			       ScalorType * veloy,
			       ScalorType * veloz,
			       ScalorType totalMassi,
			       ScalorType * sums);


// needs ceil(numAtom/blockDim.x) blocks
__global__ void leapFrog1Step (const IndexType numAtom,
			       const ScalorType * massi,
#ifndef COORD_IN_ONE_VEC
			       ScalorType * coordx,
			       ScalorType * coordy, 
			       ScalorType * coordz,
#else
			       CoordType * coord,
#endif
			       ScalorType * velox,
			       ScalorType * veloy, 
			       ScalorType * veloz,
			       const ScalorType * forcx,
			       const ScalorType * forcy, 
			       const ScalorType * forcz,
			       const ScalorType dt);
__global__ void leapFrog1Step (const IndexType numAtom,
			       const ScalorType * mass,
			       const ScalorType * massi,
#ifndef COORD_IN_ONE_VEC
			       ScalorType * coordx,
			       ScalorType * coordy, 
			       ScalorType * coordz,
#else
			       CoordType * coord,
#endif
			       ScalorType * velox,
			       ScalorType * veloy, 
			       ScalorType * veloz,
			       const ScalorType * forcx,
			       const ScalorType * forcy, 
			       const ScalorType * forcz,
			       const ScalorType dt,
			       ScalorType * statistic_buffxx,
			       ScalorType * statistic_buffyy,
			       ScalorType * statistic_buffzz);
__global__ void leapFrogStepX (const IndexType numAtom,
			       const ScalorType * massi,
#ifndef COORD_IN_ONE_VEC
			       ScalorType * coordx,
			       ScalorType * coordy, 
			       ScalorType * coordz,
#else
			       CoordType * coord,
#endif
			       const ScalorType * velox,
			       const ScalorType * veloy, 
			       const ScalorType * veloz,
			       const ScalorType dt);
__global__ void leapFrogStepV (const IndexType numAtom,
			       const ScalorType * massi,
			       ScalorType * velox,
			       ScalorType * veloy, 
			       ScalorType * veloz,
			       const ScalorType * forcx,
			       const ScalorType * forcy, 
			       const ScalorType * forcz,
			       const ScalorType dt);
__global__ void leapFrogStepV (const IndexType numAtom,
			       const ScalorType * mass,
			       const ScalorType * massi,
			       ScalorType * velox,
			       ScalorType * veloy, 
			       ScalorType * veloz,
			       const ScalorType * forcx,
			       const ScalorType * forcy, 
			       const ScalorType * forcz,
			       const ScalorType dt,
			       ScalorType * statistic_buffxx,
			       ScalorType * statistic_buffyy,
			       ScalorType * statistic_buffzz);

__global__ void leapFrogStepV_VCouple (const IndexType numAtom,
				       const ScalorType * massi,
				       ScalorType * velox,
				       ScalorType * veloy, 
				       ScalorType * veloz,
				       const ScalorType * forcx,
				       const ScalorType * forcy, 
				       const ScalorType * forcz,
				       const ScalorType lambda0,
				       const ScalorType lambda1,
				       const ScalorType lambda2,
				       const ScalorType dt);
__global__ void leapFrogStepV_VCouple (const IndexType numAtom,
				       const ScalorType * mass,
				       const ScalorType * massi,
				       ScalorType * velox,
				       ScalorType * veloy, 
				       ScalorType * veloz,
				       const ScalorType * forcx,
				       const ScalorType * forcy, 
				       const ScalorType * forcz,
				       const ScalorType lambda0,
				       const ScalorType lambda1,
				       const ScalorType lambda2,
				       const ScalorType dt,
				       ScalorType * statistic_buffxx,
				       ScalorType * statistic_buffyy,
				       ScalorType * statistic_buffzz);


// needs ceil(numAtom/blockDim.x) blocks
__global__ void velocityVerlet_part1 (const IndexType numAtom,
				      const ScalorType * massi,
#ifndef COORD_IN_ONE_VEC
				      ScalorType * coordx,
				      ScalorType * coordy, 
				      ScalorType * coordz,
#else
				      CoordType * coord,
#endif
				      ScalorType * velox,
				      ScalorType * veloy, 
				      ScalorType * veloz,
				      const ScalorType * forcx,
				      const ScalorType * forcy, 
				      const ScalorType * forcz,
				      const ScalorType dt);
__global__ void velocityVerlet_part2 (const IndexType numAtom,
				      const ScalorType * massi,
				      ScalorType * velox,
				      ScalorType * veloy, 
				      ScalorType * veloz,
				      const ScalorType * forcx,
				      const ScalorType * forcy, 
				      const ScalorType * forcz,
				      const ScalorType dt);
__global__ void velocityVerlet_part2 (const IndexType numAtom,
				      const ScalorType * mass,
				      const ScalorType * massi,
				      ScalorType * velox,
				      ScalorType * veloy, 
				      ScalorType * veloz,
				      const ScalorType * forcx,
				      const ScalorType * forcy, 
				      const ScalorType * forcz,
				      const ScalorType dt,
				      ScalorType * statistic_buffxx,
				      ScalorType * statistic_buffyy,
				      ScalorType * statistic_buffzz);
__global__ void velocityVerlet_part2a (const IndexType numAtom,
				      const ScalorType * mass,
				      const ScalorType * massi,
				      ScalorType * velox,
				      ScalorType * veloy, 
				      ScalorType * veloz,
				      const ScalorType * forcx,
				      const ScalorType * forcy, 
				      const ScalorType * forcz,
				      const ScalorType dt,
				      ScalorType * statistic_buff);
// needs ceil(numAtom/blockDim.x) blocks
__global__ void velocityRescale_rescale (const IndexType numAtom,
					 ScalorType * velox,
					 ScalorType * veloy, 
					 ScalorType * veloz,
					 const ScalorType alpha);
__global__ void velocityRescale_rescale (const IndexType numAtom,
					 const ScalorType * mass,
					 ScalorType * velox,
					 ScalorType * veloy, 
					 ScalorType * veloz,
					 const ScalorType alpha,
					 ScalorType * statistic_buffxx,
					 ScalorType * statistic_buffyy,
					 ScalorType * statistic_buffzz);



#endif
