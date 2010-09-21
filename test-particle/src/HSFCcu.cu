#include "HSFC.h"
#include "common.h"


void initDeviceHSFCMap1dto3d (HSFCMap1dto3d *dmap,
			      HSFCMap1dto3d hmap,
			      IndexType nx, IndexType ny, IndexType nz )
{
  cudaMalloc ((void**)dmap, sizeof(IndexType) *nx*ny*nz * 3);
  cudaMemcpy (*dmap, hmap, sizeof(IndexType) *nx*ny*nz * 3,
	      cudaMemcpyHostToDevice);
}

void initDeviceHSFCMap3dto1d (HSFCMap3dto1d *dmap,
			      HSFCMap3dto1d hmap,
			      IndexType nx, IndexType ny, IndexType nz )
{
  cudaMalloc ((void**)dmap, sizeof(IndexType) *nx*ny*nz * 3);
  cudaMemcpy (*dmap, hmap, sizeof(IndexType) *nx*ny*nz * 3,
	      cudaMemcpyHostToDevice);
}


