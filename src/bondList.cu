#define DEVICE_CODE

#include "BondList.h"

void initDeviceBondList (DeviceBondList & dbdlist)
{
  dbdlist.malloced = false;
  dbdlist.stride = 0;
  dbdlist.maxNumBond = 0;
}

void destroyDeviceBondList(DeviceBondList &dbdlist )
{
  if (dbdlist.malloced) {
    cudaFree (dbdlist.bondNeighborIndex);
    cudaFree (dbdlist.bondIndex);
    cudaFree (dbdlist.numBond);
    checkCUDAError ("destroyDeviceBondList");
    dbdlist.malloced = false;
    dbdlist.stride = 0;
    dbdlist.maxNumBond = 0;
  }
}

void mallocDeviceBondList (const HostBondList & hbdlist,
			   DeviceBondList & dbdlist)
{
  if (dbdlist.malloced){
    destroyDeviceBondList (dbdlist);
  }
  
  dbdlist.stride = hbdlist.stride;
  dbdlist.maxNumBond = hbdlist.maxNumBond;
  
  cudaMalloc ((void**)&(dbdlist.bondNeighborIndex), 
	      sizeof(IndexType) * hbdlist.stride * hbdlist.maxNumBond);
  checkCUDAError ("buildDeviceBondList malloc bondNeighborIndex");
  cudaMalloc ((void**)&(dbdlist.bondIndex),
	      sizeof(TypeType) * hbdlist.stride * hbdlist.maxNumBond);
  checkCUDAError ("buildDeviceBondList malloc bondIndex");
  cudaMalloc ((void**)&(dbdlist.numBond),
	      sizeof(IndexType) * hbdlist.stride);
  checkCUDAError ("buildDeviceBondList malloc numBond");
  
  dbdlist.malloced = true;
}

void copyDeviceBondList (const HostBondList & hbdlist,
			 DeviceBondList & dbdlist)
{
  if (dbdlist.malloced){
    cudaMemcpy (dbdlist.bondNeighborIndex, hbdlist.bondNeighborIndex,
		sizeof(IndexType) * hbdlist.stride * hbdlist.maxNumBond,
		cudaMemcpyHostToDevice);
    checkCUDAError ("buildDeviceBondList cpy host bondNeighborIndex to device");
    cudaMemcpy (dbdlist.bondIndex, hbdlist.bondIndex,
		sizeof(IndexType) * hbdlist.stride * hbdlist.maxNumBond,
		cudaMemcpyHostToDevice);
    checkCUDAError ("buildDeviceBondList cpy host bondIndex to device");
    cudaMemcpy (dbdlist.numBond, hbdlist.numBond,
		sizeof(IndexType) * hbdlist.stride,
		cudaMemcpyHostToDevice);
    checkCUDAError ("buildDeviceBondList cpy host numBond to device");
  }
}

void copyDeviceBondList (const DeviceBondList & dbdlist1,
			 DeviceBondList & dbdlist)
{
  if (dbdlist.malloced && dbdlist1.malloced){
    cudaMemcpy (dbdlist.bondNeighborIndex, dbdlist1.bondNeighborIndex,
		sizeof(IndexType) * dbdlist1.stride * dbdlist1.maxNumBond,
		cudaMemcpyDeviceToDevice);
    checkCUDAError ("buildDeviceBondList cpy device bondNeighborIndex to device");
    cudaMemcpy (dbdlist.bondIndex, dbdlist1.bondIndex,
		sizeof(IndexType) * dbdlist1.stride * dbdlist1.maxNumBond,
		cudaMemcpyDeviceToDevice);
    checkCUDAError ("buildDeviceBondList cpy device bondIndex to device");
    cudaMemcpy (dbdlist.numBond, dbdlist1.numBond,
		sizeof(IndexType) * dbdlist1.stride,
		cudaMemcpyDeviceToDevice);
    checkCUDAError ("buildDeviceBondList cpy device numBond to device");
  }
}
  

HostBondList::HostBondList ()
{
  stride = 0;
  maxNumBond = 0;
  bondNeighborIndex = NULL;
  bondIndex = NULL;
  numBond = NULL;
}

HostBondList::~HostBondList()
{
  clearMem();
}


void HostBondList::
clearMem ()
{
  stride = 0;
  maxNumBond = 0;
  freeAPointer ((void**)&bondNeighborIndex);
  freeAPointer ((void**)&bondIndex);
  freeAPointer ((void**)&numBond);
}

void HostBondList::
clearBond ()
{
  for (unsigned i = 0; i < stride; ++i){
    for (unsigned j = 0; j < maxNumBond; ++j){
      bondNeighborIndex[j * stride + i] = MaxIndexValue;
      bondIndex        [j * stride + i] = MaxIndexValue;
    }
    numBond[i] = 0;
  }
}


void HostBondList::reinit (const IndexType & stride_,
			   const IndexType & maxNumBond_)
{
  clearMem();
  
  stride = stride_;
  maxNumBond = maxNumBond_;
  
  bondNeighborIndex = (IndexType *) malloc (sizeof(IndexType) * stride * maxNumBond);
  if (bondNeighborIndex == NULL){
    throw MDExcptFailedMallocOnHost ("HostBondList::init", "bondNeighborIndex",
				     sizeof(IndexType) * stride * maxNumBond);
  }
  bondIndex = (IndexType *) malloc
      (sizeof(IndexType *) * stride * maxNumBond);
  if (bondIndex == NULL){
    throw MDExcptFailedMallocOnHost ("HostBondList::init", "bondIndex",
				     sizeof(IndexType *) * stride * maxNumBond);
  }
  numBond = (IndexType *) malloc (sizeof(IndexType) * stride);
  if (numBond == NULL){
    throw MDExcptFailedMallocOnHost ("HostBondList::init", "numBond",
				     sizeof(IndexType) * stride);
  }
  
  clearBond ();
}


void HostBondList::addBond (const IndexType & ii,
			    const IndexType & jj,
			    const IndexType & bondIndex_)
{
  bondNeighborIndex[numBond[ii] * stride + ii] = jj;
  bondIndex[numBond[ii] * stride + ii] = bondIndex_;
  numBond[ii] ++;
  // bondNeighborIndex[numBond[jj] * stride + jj] = ii;
  // bondIndex[numBond[jj] * stride + jj] = bondIndex_;
  // numBond[jj] ++;
}

