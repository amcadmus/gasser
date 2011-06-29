#define DEVICE_CODE

#include "AngleList.h"

void initDeviceAngleList (DeviceAngleList & dbdlist)
{
  dbdlist.malloced = false;
  dbdlist.stride = 0;
  dbdlist.maxNumAngle = 0;
}

void destroyDeviceAngleList(DeviceAngleList &dbdlist )
{
  if (dbdlist.malloced) {
    cudaFree (dbdlist.angleNeighborIndex);
    cudaFree (dbdlist.angleIndex);
    cudaFree (dbdlist.anglePosi);
    cudaFree (dbdlist.numAngle);
    checkCUDAError ("destroyDeviceAngleList");
    dbdlist.malloced = false;
    dbdlist.stride = 0;
    dbdlist.maxNumAngle = 0;
  }
}

void mallocDeviceAngleList (const HostAngleList & hbdlist,
			    DeviceAngleList & dbdlist)
{
  if (dbdlist.malloced){
    destroyDeviceAngleList (dbdlist);
  }
  
  dbdlist.stride = hbdlist.stride;
  dbdlist.maxNumAngle = hbdlist.maxNumAngle;
  
  cudaMalloc ((void**)&(dbdlist.angleNeighborIndex), 
	      sizeof(IndexType) * hbdlist.stride * hbdlist.maxNumAngle * 2);
  checkCUDAError ("buildDeviceAngleList malloc angleNeighborIndex");
  cudaMalloc ((void**)&(dbdlist.angleIndex),
	      sizeof(TypeType) * hbdlist.stride * hbdlist.maxNumAngle);
  checkCUDAError ("buildDeviceAngleList malloc angleIndex");
  cudaMalloc ((void**)&(dbdlist.anglePosi),
	      sizeof(TypeType) * hbdlist.stride * hbdlist.maxNumAngle);
  checkCUDAError ("buildDeviceAngleList malloc anglePosi");
  cudaMalloc ((void**)&(dbdlist.numAngle),
	      sizeof(IndexType) * hbdlist.stride);
  checkCUDAError ("buildDeviceAngleList malloc numAngle");
  
  dbdlist.malloced = true;
}

void copyDeviceAngleList (const HostAngleList & hbdlist,
			  DeviceAngleList & dbdlist)
{
  if (dbdlist.malloced){
    cudaMemcpy (dbdlist.angleNeighborIndex, hbdlist.angleNeighborIndex,
		sizeof(IndexType) * hbdlist.stride * hbdlist.maxNumAngle * 2,
		cudaMemcpyHostToDevice);
    checkCUDAError ("buildDeviceAngleList cpy host angleNeighborIndex to device");
    cudaMemcpy (dbdlist.angleIndex, hbdlist.angleIndex,
		sizeof(IndexType) * hbdlist.stride * hbdlist.maxNumAngle,
		cudaMemcpyHostToDevice);
    checkCUDAError ("buildDeviceAngleList cpy host angleIndex to device");
    cudaMemcpy (dbdlist.anglePosi, hbdlist.anglePosi,
		sizeof(IndexType) * hbdlist.stride * hbdlist.maxNumAngle,
		cudaMemcpyHostToDevice);
    checkCUDAError ("buildDeviceAngleList cpy host anglePosi to device");
    cudaMemcpy (dbdlist.numAngle, hbdlist.numAngle,
		sizeof(IndexType) * hbdlist.stride,
		cudaMemcpyHostToDevice);
    checkCUDAError ("buildDeviceAngleList cpy host numAngle to device");
  }
}

void copyDeviceAngleList (const DeviceAngleList & dbdlist1,
			  DeviceAngleList & dbdlist)
{
  if (dbdlist.malloced && dbdlist1.malloced){
    cudaMemcpy (dbdlist.angleNeighborIndex, dbdlist1.angleNeighborIndex,
		sizeof(IndexType) * dbdlist1.stride * dbdlist1.maxNumAngle * 2,
		cudaMemcpyDeviceToDevice);
    checkCUDAError ("buildDeviceAngleList cpy device angleNeighborIndex to device");
    cudaMemcpy (dbdlist.angleIndex, dbdlist1.angleIndex,
		sizeof(IndexType) * dbdlist1.stride * dbdlist1.maxNumAngle,
		cudaMemcpyDeviceToDevice);
    checkCUDAError ("buildDeviceAngleList cpy device angleIndex to device");
    cudaMemcpy (dbdlist.anglePosi, dbdlist1.anglePosi,
		sizeof(IndexType) * dbdlist1.stride * dbdlist1.maxNumAngle,
		cudaMemcpyDeviceToDevice);
    checkCUDAError ("buildDeviceAngleList cpy device anglePosi to device");
    cudaMemcpy (dbdlist.numAngle, dbdlist1.numAngle,
		sizeof(IndexType) * dbdlist1.stride,
		cudaMemcpyDeviceToDevice);
    checkCUDAError ("buildDeviceAngleList cpy device numAngle to device");
  }
}


HostAngleList::HostAngleList ()
{
  stride = 0;
  maxNumAngle = 0;
  angleNeighborIndex = NULL;
  angleIndex = NULL;
  anglePosi = NULL;
  numAngle = NULL;
}

HostAngleList::~HostAngleList()
{
  clearMem();
}


void HostAngleList::
clearMem ()
{
  stride = 0;
  maxNumAngle = 0;
  freeAPointer ((void**)&angleNeighborIndex);
  freeAPointer ((void**)&angleIndex);
  freeAPointer ((void**)&anglePosi);
  freeAPointer ((void**)&numAngle);
}

void HostAngleList::
clearAngle ()
{
  for (unsigned i = 0; i < stride; ++i){
    for (unsigned j = 0; j < maxNumAngle; ++j){
      angleNeighborIndex[2 * j * stride + i] = MaxIndexValue;
      angleNeighborIndex[2 * j * stride + i + stride] = MaxIndexValue;
      angleIndex        [j * stride + i] = MaxIndexValue;
      anglePosi		[j * stride + i] = MaxIndexValue;
    }
    numAngle[i] = 0;
  }
}


void HostAngleList::reinit (const IndexType & stride_,
			    const IndexType & maxNumAngle_)
{
  clearMem();
  
  stride = stride_;
  maxNumAngle = maxNumAngle_;
  
  angleNeighborIndex = (IndexType *) malloc (sizeof(IndexType) * stride * 2 * maxNumAngle);
  if (angleNeighborIndex == NULL){
    throw MDExcptFailedMallocOnHost ("HostAngleList::init", "angleNeighborIndex",
				     sizeof(IndexType) * stride * 2 * maxNumAngle);
  }
  angleIndex = (IndexType *) malloc
      (sizeof(IndexType *) * stride * maxNumAngle);
  if (angleIndex == NULL){
    throw MDExcptFailedMallocOnHost ("HostAngleList::init", "angleIndex",
				     sizeof(IndexType *) * stride * maxNumAngle);
  }
  anglePosi = (IndexType *) malloc
      (sizeof(IndexType *) * stride * maxNumAngle);
  if (anglePosi == NULL){
    throw MDExcptFailedMallocOnHost ("HostAngleList::init", "anglePosi",
				     sizeof(IndexType *) * stride * maxNumAngle);
  }
  numAngle = (IndexType *) malloc (sizeof(IndexType) * stride);
  if (numAngle == NULL){
    throw MDExcptFailedMallocOnHost ("HostAngleList::init", "numAngle",
				     sizeof(IndexType) * stride);
  }
  
  clearAngle ();
}


void HostAngleList::
addAngle (const IndexType &ii,
	  const IndexType &jj,
	  const IndexType &kk,
	  const IndexType &angleIndex_,
	  const IndexType &anglePosi_)
{
  IndexType index = numAngle[ii] * stride + ii;
  IndexType index1 = 2*numAngle[ii] *stride + ii;
  angleNeighborIndex[index1] = jj;
  angleNeighborIndex[index1+stride] = kk;
  angleIndex[index] = angleIndex_;
  anglePosi[index] = anglePosi_;
  numAngle[ii] ++;
}

