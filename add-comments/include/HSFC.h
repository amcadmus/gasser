#ifndef __HSFC_h_wanghan__
#define __HSFC_h_wanghan__

#include "common.h"

typedef IndexType *	HSFCMap1dto3d;
typedef IndexType *	HSFCMap3dto1d;;

void map1dto3d (HSFCMap1dto3d map, IndexType index,
		IndexType * i, IndexType * j, IndexType * k);
void map3dto1d (HSFCMap3dto1d map, 
		IndexType nx, IndexType ny, IndexType nz,
		IndexType i,  IndexType j,  IndexType k,
		IndexType *index);

#ifndef CPP_FILE
__device__ void map1dto3d_d (HSFCMap1dto3d map, IndexType index,
			     IndexType * i, IndexType * j, IndexType * k);
__device__ void map3dto1d_d (HSFCMap3dto1d map, 
			       IndexType nx, IndexType ny, IndexType nz,
		       IndexType i,  IndexType j,  IndexType k,
		       IndexType *index);
#endif

void initHostHSFCMap1dto3d (HSFCMap1dto3d *map, 
			    IndexType nx, IndexType ny, IndexType nz );
void initHostHSFCMap3dto1d (HSFCMap3dto1d *map,
			    HSFCMap1dto3d map1dto3d,
			    IndexType nx, IndexType ny, IndexType nz );
void initDeviceHSFCMap1dto3d (HSFCMap1dto3d *dmap,
			      HSFCMap1dto3d hmap,
			      IndexType nx, IndexType ny, IndexType nz );
void initDeviceHSFCMap3dto1d (HSFCMap3dto1d * dmap,
			      HSFCMap3dto1d hmap,
			      IndexType nx, IndexType ny, IndexType nz );
			      

////////////////////////////////////////////////////////////
// inline implementation
////////////////////////////////////////////////////////////

inline void map1dto3d (HSFCMap1dto3d map, IndexType index,
		       IndexType * i, IndexType * j, IndexType * k)

{
  *i = map[index * 3];
  *j = map[index * 3 + 1];
  *k = map[index * 3 + 2];
}

inline void map3dto1d (HSFCMap3dto1d map, 
		       IndexType nx, IndexType ny, IndexType nz,
		       IndexType i,  IndexType j,  IndexType k,
		       IndexType *index)
{
  *index = map[i * ny * nz + j * nz + k];
}

#ifndef CPP_FILE
__device__ void map1dto3d_d (HSFCMap1dto3d map, IndexType index,
		       IndexType * i, IndexType * j, IndexType * k)

{
  *i = map[index * 3];
  *j = map[index * 3 + 1];
  *k = map[index * 3 + 2];
}

__device__ void map3dto1d_d (HSFCMap3dto1d map, 
		       IndexType nx, IndexType ny, IndexType nz,
		       IndexType i,  IndexType j,  IndexType k,
		       IndexType *index)
{
  *index = map[i * ny * nz + j * nz + k];
}
#endif

#endif 
