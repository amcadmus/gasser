#define DEVICE_CODE

#include "ExclusionList.h"

void
initDeviceExclusionList   (DeviceExclusionList & dexcllist)
{
  dexcllist.malloced = false;
  dexcllist.stride = 0;
  dexcllist.maxNumExclusion = 0;
}

void
destroyDeviceBondList     (DeviceExclusionList & dexcllist)
{
  if (dexcllist.malloced){
    cudaFree (dexcllist.exclusionNeighborIndex);
    cudaFree (dexcllist.numExclusion);
    checkCUDAError ("destroyDeviceBondList");
    dexcllist.malloced = false;
    dexcllist.stride = 0;
    dexcllist.maxNumExclusion = 0;    
  }
}


void
mallocDeviceExclusionList (const HostExclusionList & hexcllist,
				DeviceExclusionList & dexcllist)
{
  destroyDeviceBondList (dexcllist);

  dexcllist.stride = hexcllist.stride;
  dexcllist.maxNumExclusion = hexcllist.maxNumExclusion;

  size_t memSize = dexcllist.stride * sizeof(IndexType);
  cudaMalloc ((void**)&(dexcllist.exclusionNeighborIndex),
	      memSize * dexcllist.maxNumExclusion);
  checkCUDAError ("mallocDeviceExclusionList exclusionNeighborIndex");
  cudaMalloc ((void**)&(dexcllist.numExclusion), memSize);
  checkCUDAError ("mallocDeviceExclusionList numExclusion");

  dexcllist.malloced = true;
}

void
copyDeviceExclusionList   (const HostExclusionList & hexcllist,
			   DeviceExclusionList & dexcllist) 
{
  if (dexcllist.malloced){
    size_t memSize = dexcllist.stride * sizeof(IndexType);
    cudaMemcpy (dexcllist.exclusionNeighborIndex,
		hexcllist.exclusionNeighborIndex,
		memSize * dexcllist.maxNumExclusion,
		cudaMemcpyHostToDevice);
    checkCUDAError ("copyDeviceExclusionList exclusionNeighborIndex");
    cudaMemcpy (dexcllist.numExclusion,
		hexcllist.numExclusion,
		memSize,
		cudaMemcpyHostToDevice);
    checkCUDAError ("copyDeviceExclusionList numExclusion");
  }
}


void
copyDeviceExclusionList   (const DeviceExclusionList & dexcllist1,
			   DeviceExclusionList & dexcllist) 
{
  if (dexcllist.malloced){
    size_t memSize = dexcllist.stride * sizeof(IndexType);
    cudaMemcpy (dexcllist.exclusionNeighborIndex,
		dexcllist1.exclusionNeighborIndex,
		memSize * dexcllist.maxNumExclusion,
		cudaMemcpyDeviceToDevice);
    checkCUDAError ("copyDeviceExclusionList exclusionNeighborIndex");
    cudaMemcpy (dexcllist.numExclusion,
		dexcllist1.numExclusion,
		memSize,
		cudaMemcpyDeviceToDevice);
    checkCUDAError ("copyDeviceExclusionList numExclusion");
  }
}
		

void HostExclusionList::
clearMem ()
{
  stride = 0;
  maxNumExclusion = 0;
  freeAPointer ((void**)&exclusionNeighborIndex);
  freeAPointer ((void**)&numExclusion);
}

HostExclusionList::
HostExclusionList()
    : stride(0), maxNumExclusion(0),
      exclusionNeighborIndex(NULL),
      numExclusion(NULL)
{
}

HostExclusionList::
~HostExclusionList()
{
  clearMem();
}

void HostExclusionList::
reinit (const IndexType & stride_,
	const IndexType & maxNumExclusion_)
{
  clearMem();
  
  stride = stride_;
  maxNumExclusion = maxNumExclusion_;

  size_t memSize = stride * sizeof(IndexType);
  exclusionNeighborIndex = (IndexType *) malloc (memSize * maxNumExclusion);
  if (exclusionNeighborIndex == NULL){
    throw MDExcptFailedMallocOnHost ("HostExclusionList::init",
				     "exclusionNeighborIndex",
				     memSize * maxNumExclusion);
  }
  numExclusion = (IndexType *) malloc (memSize);
  if (numExclusion == NULL){
    throw MDExcptFailedMallocOnHost ("HostExclusionList::init",
				     "numExclusion",
				     memSize);
  }

  clearExclusion();
}

void HostExclusionList::
clearExclusion ()
{
  for (IndexType i = 0; i < maxNumExclusion; ++i){
    for (IndexType j = 0; j < stride; ++j){
      exclusionNeighborIndex[i*stride + j] = MaxIndexValue;
    }
  }
  for (IndexType j = 0; j < stride; ++j){
    numExclusion[j] = 0;
  }
}


void HostExclusionList::
addExclusion (const IndexType & ii,
		   const IndexType & jj)
{
  exclusionNeighborIndex[numExclusion[ii] * stride + ii] = jj;
  numExclusion[ii] ++;
}


ExclusionList::
ExclusionList ()
{
  initDeviceExclusionList(dexcllist);
  initDeviceExclusionList(bkdexcllist);
}

ExclusionList::
~ExclusionList()
{
  destroyDeviceBondList (dexcllist);
  destroyDeviceBondList (bkdexcllist);
}

void ExclusionList::
reinit (const MDSystem & sysData,
	const Topology::System & sysTop,
	const SystemNonBondedInteraction & sysNbInter)
{
  IndexType stride = sysData.hdata.numAtom;
  IndexType maxNumExclusion = sysNbInter.maxNumberOfExclusion();
  
  hexcllist.reinit (stride, maxNumExclusion);

  IndexType shift0 = 0;
  for (unsigned i = 0; i < sysTop.molecules.size(); ++i){
    shift0 = sysTop.indexShift[i];
    IndexType molSize = sysTop.molecules[i].size();
    for (unsigned j = 0; j < sysTop.numbers[i]; ++j){
      IndexType shift1 = j * molSize;
      IndexType indexSift = shift0 + shift1;
      for (unsigned k = 0; k < molSize; ++k){
	for (unsigned l = 0;
	     l < sysNbInter.exclusionNeighborIndex[i][k].size();
	     ++l){
	  hexcllist.addExclusion (indexSift + k,
				  indexSift + sysNbInter.exclusionNeighborIndex[i][k][l]);
	}
      }
    }
  }

  destroyDeviceBondList (dexcllist);
  destroyDeviceBondList (bkdexcllist);
  mallocDeviceExclusionList (hexcllist, dexcllist);
  mallocDeviceExclusionList (hexcllist, bkdexcllist);
  copyDeviceExclusionList (hexcllist, dexcllist);
}

static __global__ void
Reshuffle_reshuffleExclusionList (const IndexType numAtom,
				  const IndexType * excllistData2,
				  const IndexType * excllistNumB2,
				  const IndexType stride,
				  const IndexType * idxTable,
				  IndexType * excllistData,
				  IndexType * excllistNumB)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  
  if (ii < numAtom){
    IndexType toid = idxTable[ii];
    excllistNumB[toid] = excllistNumB2[ii];
    for (IndexType jj = 0; jj < excllistNumB2[ii]; ++jj){
      excllistData[jj * stride + toid] = 
	  idxTable[excllistData2[jj * stride + ii]];
    }
  }
}

void ExclusionList::
reshuffle (const IndexType * indexTable,
	   const IndexType & numAtom,
	   MDTimer * timer) 
{
  if (timer != NULL) timer->tic(mdTimeReshuffleSystem);

  copyDeviceExclusionList (dexcllist, bkdexcllist);

  dim3 myBlockDim, atomGridDim;
  myBlockDim.y = 1;
  myBlockDim.z = 1;
  myBlockDim.x = DefaultNThreadPerBlock;
  IndexType nob;
  if (numAtom % myBlockDim.x == 0){
    nob = numAtom / myBlockDim.x;
  } else {
    nob = numAtom / myBlockDim.x + 1;
  }
  atomGridDim = toGridDim (nob);

  Reshuffle_reshuffleExclusionList
      <<<atomGridDim,myBlockDim>>> (
  	  numAtom,
  	  bkdexcllist.exclusionNeighborIndex,
	  bkdexcllist.numExclusion,
  	  dexcllist.stride,
  	  indexTable,
  	  dexcllist.exclusionNeighborIndex, 
	  dexcllist.numExclusion);
  checkCUDAError ("ExclusionList::reshuffle, reshuffle exclusion list");
  
  if (timer != NULL) timer->toc(mdTimeReshuffleSystem);
}

