#define DEVICE_CODE

#include "MDSystem.h"
#include "common.h"

#include "MDSystem.h"
#include <stdlib.h>

HostMDData::HostMDData()
{
  numAtom = 0;
  numMem = 0;
  NFreedom = 0;
  totalMass = totalMassi = 0.f;
#ifndef COORD_IN_ONE_VEC
  coordx = coordy = coordz = NULL;
#else
  coord = NULL;
#endif
  coordNoix = coordNoiy = coordNoiz = NULL;
  velox = veloy = veloz = NULL;
  forcx = forcy = forcz = NULL;
  type = NULL;
  mass = massi = NULL;
  charge = NULL;
  atomName = NULL;
  atomIndex = NULL;
  resdName = NULL;
  resdIndex = NULL;
}

HostMDData::~HostMDData()
{
  destroyHostMDData (this);
}

DeviceMDData::DeviceMDData()
{
  malloced = false;
  numAtom = 0;
  numMem = 0;
  NFreedom = 0;
  totalMass = totalMassi = 0.f;
}

DeviceMDData::~DeviceMDData()
{
  destroyDeviceMDData(this);
}  


__host__ void mallocHostMDData (IndexType numAtom, IndexType expectedMaxNumAtom,
				HostMDData * hdata)
{
  hdata->numAtom = numAtom;
  hdata->numMem = expectedMaxNumAtom;
  size_t sizef = hdata->numMem * sizeof(ScalorType);
#ifdef COORD_IN_ONE_VEC
  size_t sizecoord =hdata->numMem * sizeof(CoordType);
#endif
  size_t sizei = hdata->numMem * sizeof(IntScalorType);
  size_t sizec = hdata->numMem * sizeof(char) * StringSize;
  size_t sizeidx = hdata->numMem * sizeof(IndexType);

#ifndef COORD_IN_ONE_VEC
  hdata->coordx = (ScalorType *) malloc (sizef);
  if (hdata->coordx == NULL) throw MDExcptFailedMallocOnHost("hdata->coordx", sizef);
  hdata->coordy = (ScalorType *) malloc (sizef);
  if (hdata->coordy == NULL) throw MDExcptFailedMallocOnHost("hdata->coordy", sizef);
  hdata->coordz = (ScalorType *) malloc (sizef);
  if (hdata->coordz == NULL) throw MDExcptFailedMallocOnHost("hdata->coordz", sizef);
#else
  hdata->coord = (CoordType *) malloc (sizecoord);
  if (hdata->coord == NULL) throw (MDExcptFailedMallocOnHost("hdata->coord", sizecoord));
#endif
  
  hdata->coordNoix = (IntScalorType *) malloc (sizei);
  if (hdata->coordNoix == NULL) throw MDExcptFailedMallocOnHost ("hdata->coordNoix", sizei);
  hdata->coordNoiy = (IntScalorType *) malloc (sizei);
  if (hdata->coordNoiy == NULL) throw MDExcptFailedMallocOnHost ("hdata->coordNoiy", sizei);
  hdata->coordNoiz = (IntScalorType *) malloc (sizei);
  if (hdata->coordNoiz == NULL) throw MDExcptFailedMallocOnHost ("hdata->coordNoiz", sizei);
  
  hdata->velox = (ScalorType *) malloc (sizef);
  if (hdata->velox == NULL) throw MDExcptFailedMallocOnHost("hdata->velox", sizef);
  hdata->veloy = (ScalorType *) malloc (sizef);
  if (hdata->veloy == NULL) throw MDExcptFailedMallocOnHost("hdata->veloy", sizef);
  hdata->veloz = (ScalorType *) malloc (sizef);
  if (hdata->veloz == NULL) throw MDExcptFailedMallocOnHost("hdata->veloz", sizef);
  
  hdata->forcx = (ScalorType *) malloc (sizef);
  if (hdata->forcx == NULL) throw MDExcptFailedMallocOnHost("hdata->forcx", sizef);
  hdata->forcy = (ScalorType *) malloc (sizef);
  if (hdata->forcy == NULL) throw MDExcptFailedMallocOnHost("hdata->forcy", sizef);
  hdata->forcz = (ScalorType *) malloc (sizef);
  if (hdata->forcz == NULL) throw MDExcptFailedMallocOnHost("hdata->forcz", sizef);

  hdata->type = (TypeType *) malloc (hdata->numMem * sizeof(TypeType));
  if (hdata->type == NULL) {
    throw MDExcptFailedMallocOnHost("hdata->type", hdata->numMem * sizeof(TypeType));
  }
  hdata->mass = (ScalorType *) malloc (sizef);
  if (hdata->mass == NULL) throw MDExcptFailedMallocOnHost("hdata->mass", sizef);
  hdata->massi = (ScalorType *) malloc (sizef);
  if (hdata->massi == NULL) throw MDExcptFailedMallocOnHost("hdata->massi", sizef);
  hdata->charge = (ScalorType *) malloc(sizef);
  if (hdata->charge == NULL) throw MDExcptFailedMallocOnHost("hdata->charge", sizef);

  hdata->atomName = (char *) malloc (sizec);
  if (hdata->atomName == NULL) throw MDExcptFailedMallocOnHost("hdata->atomName", sizec);
  hdata->atomIndex = (IndexType *) malloc (sizeidx);
  if (hdata->atomIndex == NULL) throw MDExcptFailedMallocOnHost("hdata->atomIndex", sizeidx);
  hdata->resdName = (char *) malloc (sizec);
  if (hdata->resdName == NULL) throw MDExcptFailedMallocOnHost("hdata->resdName", sizec);
  hdata->resdIndex = (IndexType *) malloc (sizeidx); 
  if (hdata->resdIndex == NULL) throw MDExcptFailedMallocOnHost("hdata->resdIndex", sizeidx);
}

__host__ void lazyInitHostMDData (HostMDData * hdata)
{
  hdata->NFreedom = 3 * hdata->numAtom;
  for (IndexType i = 0; i < hdata->numAtom; ++i){
#ifndef COORD_IN_ONE_VEC
    hdata->coordx[i] = 0.;
    hdata->coordy[i] = 0.;
    hdata->coordz[i] = 0.;
#else
    hdata->coord[i].x = 0.f;
    hdata->coord[i].y = 0.f;
    hdata->coord[i].z = 0.f;
#endif
    
    hdata->coordNoix[i] = 0;
    hdata->coordNoiy[i] = 0;
    hdata->coordNoiz[i] = 0;

    hdata->velox[i] = 0.;
    hdata->veloy[i] = 0.;
    hdata->veloz[i] = 0.;

    hdata->forcx[i] = 0.;
    hdata->forcy[i] = 0.;
    hdata->forcz[i] = 0.;
    
    hdata->type[i] = (TypeType)(i);
    hdata->mass[i] = 1.f;
    hdata->massi[i] = 1.f;
    hdata->charge[i] = 0.f;

    strcpy (&(hdata->atomName)[i*8], "AName");
    strcpy (&(hdata->resdName)[i*8], "RName");
    hdata->resdIndex[i] = i;
    hdata->atomIndex[i] = i;
  }
  initMass (hdata);
}

__host__ void initMass (HostMDData * hdata)
{
  for (IndexType i = 0; i < hdata->numAtom; ++i){
    hdata->massi[i] = 1./(hdata->mass[i]);
  }
  hdata->totalMass = 0;
  for (IndexType i = 0; i < hdata->numAtom; ++i){
    hdata->totalMass += hdata->mass[i];
  }
  hdata->totalMassi = 1./hdata->totalMass;  
}

__host__ void destroyHostMDData (HostMDData * hdata)
{
#ifndef COORD_IN_ONE_VEC
  freeAPointer ((void**)&hdata->coordx);
  freeAPointer ((void**)&hdata->coordy);
  freeAPointer ((void**)&hdata->coordz);
#else
  freeAPointer ((void**)&hdata->coord);
#endif

  freeAPointer ((void**)&hdata->coordNoix);
  freeAPointer ((void**)&hdata->coordNoiy);
  freeAPointer ((void**)&hdata->coordNoiz);
  
  freeAPointer ((void**)&hdata->velox);
  freeAPointer ((void**)&hdata->veloy);
  freeAPointer ((void**)&hdata->veloz);
  
  freeAPointer ((void**)&hdata->forcx);
  freeAPointer ((void**)&hdata->forcy);
  freeAPointer ((void**)&hdata->forcz);

  freeAPointer ((void**)&hdata->type);
  freeAPointer ((void**)&hdata->mass);
  freeAPointer ((void**)&hdata->massi);
  freeAPointer ((void**)&hdata->charge);

  freeAPointer ((void**)&hdata->atomName);
  freeAPointer ((void**)&hdata->atomIndex);
  freeAPointer ((void**)&hdata->resdName);
  freeAPointer ((void**)&hdata->resdIndex);
}




__host__ void cpyDeviceMDDataToHost (const DeviceMDData * ddata,
				     HostMDData * hdata)
{
  size_t sizef = ddata->numAtom * sizeof(ScalorType);
  size_t sizei = ddata->numAtom * sizeof(IntScalorType);
#ifdef COORD_IN_ONE_VEC
  size_t sizecoord =hdata->numMem * sizeof(CoordType);
#endif
  hdata->numAtom = ddata->numAtom;
  hdata->totalMass = ddata->totalMass;
  hdata->totalMassi = ddata->totalMassi;
  hdata->NFreedom = ddata->NFreedom;

#ifndef COORD_IN_ONE_VEC
  cudaMemcpy (hdata->coordx, ddata->coordx, sizef, cudaMemcpyDeviceToHost);
  cudaMemcpy (hdata->coordy, ddata->coordy, sizef, cudaMemcpyDeviceToHost);
  cudaMemcpy (hdata->coordz, ddata->coordz, sizef, cudaMemcpyDeviceToHost);
#else
  cudaMemcpy (hdata->coord, ddata->coord, sizecoord, cudaMemcpyDeviceToHost);
#endif
  checkCUDAError ("cpyDeviceMDDataToHost coord");
  
  cudaMemcpy (hdata->coordNoix, ddata->coordNoix, sizei, cudaMemcpyDeviceToHost);
  cudaMemcpy (hdata->coordNoiy, ddata->coordNoiy, sizei, cudaMemcpyDeviceToHost);
  cudaMemcpy (hdata->coordNoiz, ddata->coordNoiz, sizei, cudaMemcpyDeviceToHost);
  checkCUDAError ("cpyDeviceMDDataToHost coordNoi");
  
  cudaMemcpy (hdata->velox, ddata->velox, sizef, cudaMemcpyDeviceToHost);
  cudaMemcpy (hdata->veloy, ddata->veloy, sizef, cudaMemcpyDeviceToHost);
  cudaMemcpy (hdata->veloz, ddata->veloz, sizef, cudaMemcpyDeviceToHost);
  checkCUDAError ("cpyDeviceMDDataToHost velo");

  cudaMemcpy (hdata->forcx, ddata->forcx, sizef, cudaMemcpyDeviceToHost);
  cudaMemcpy (hdata->forcy, ddata->forcy, sizef, cudaMemcpyDeviceToHost);
  cudaMemcpy (hdata->forcz, ddata->forcz, sizef, cudaMemcpyDeviceToHost);
  checkCUDAError ("cpyDeviceMDDataToHost forc");

  cudaMemcpy (hdata->type, ddata->type, ddata->numAtom * sizeof(TypeType), cudaMemcpyDeviceToHost);
  cudaMemcpy (hdata->mass, ddata->mass, sizef, cudaMemcpyDeviceToHost);
  cudaMemcpy (hdata->massi, ddata->massi, sizef, cudaMemcpyDeviceToHost);
  cudaMemcpy (hdata->charge, ddata->charge, sizef, cudaMemcpyDeviceToHost);  
  checkCUDAError ("cpyDeviceMDDataToHost other");
}

__host__ void cpyHostMDDataToDevice (const HostMDData * hdata, DeviceMDData * ddata)
{
  size_t sizef = hdata->numAtom * sizeof(ScalorType);
  size_t sizei = hdata->numAtom * sizeof(IntScalorType);
#ifdef COORD_IN_ONE_VEC
  size_t sizecoord =hdata->numMem * sizeof(CoordType);
#endif
  ddata->numAtom = hdata->numAtom;
  ddata->totalMass = hdata->totalMass;
  ddata->totalMassi = hdata->totalMassi;
  ddata->NFreedom = hdata->NFreedom;
  
#ifndef COORD_IN_ONE_VEC
  cudaMemcpy (ddata->coordx, hdata->coordx, sizef, cudaMemcpyHostToDevice);
  cudaMemcpy (ddata->coordy, hdata->coordy, sizef, cudaMemcpyHostToDevice);
  cudaMemcpy (ddata->coordz, hdata->coordz, sizef, cudaMemcpyHostToDevice);
#else
  cudaMemcpy (ddata->coord, hdata->coord, sizecoord, cudaMemcpyHostToDevice);
#endif
  checkCUDAError ("cpyHostMDDataToDevice coord");

  cudaMemcpy (ddata->coordNoix, hdata->coordNoix, sizei, cudaMemcpyHostToDevice);
  cudaMemcpy (ddata->coordNoiy, hdata->coordNoiy, sizei, cudaMemcpyHostToDevice);
  cudaMemcpy (ddata->coordNoiz, hdata->coordNoiz, sizei, cudaMemcpyHostToDevice);
  checkCUDAError ("cpyHostMDDataToDevice coordNoi");
  
  cudaMemcpy (ddata->velox, hdata->velox, sizef, cudaMemcpyHostToDevice);
  cudaMemcpy (ddata->veloy, hdata->veloy, sizef, cudaMemcpyHostToDevice);
  cudaMemcpy (ddata->veloz, hdata->veloz, sizef, cudaMemcpyHostToDevice);
  checkCUDAError ("cpyHostMDDataToDevice velo");

  cudaMemcpy (ddata->forcx, hdata->forcx, sizef, cudaMemcpyHostToDevice);
  cudaMemcpy (ddata->forcy, hdata->forcy, sizef, cudaMemcpyHostToDevice);
  cudaMemcpy (ddata->forcz, hdata->forcz, sizef, cudaMemcpyHostToDevice);
  checkCUDAError ("cpyHostMDDataToDevice forc");

  cudaMemcpy (ddata->type, hdata->type, ddata->numAtom * sizeof(TypeType), cudaMemcpyHostToDevice);
  cudaMemcpy (ddata->mass, hdata->mass, sizef, cudaMemcpyHostToDevice);
  cudaMemcpy (ddata->massi, hdata->massi, sizef, cudaMemcpyHostToDevice);
  cudaMemcpy (ddata->charge, hdata->charge, sizef, cudaMemcpyHostToDevice);   
  checkCUDAError ("cpyHostMDDataToDevice other");
}

__host__ void cpyDeviceMDDataToDevice (const DeviceMDData * hdata,
				       DeviceMDData * ddata)
{
  size_t sizef = hdata->numAtom * sizeof(ScalorType);
  size_t sizei = hdata->numAtom * sizeof(IntScalorType);
#ifdef COORD_IN_ONE_VEC
  size_t sizecoord =hdata->numMem * sizeof(CoordType);
#endif
  ddata->numAtom = hdata->numAtom;
  ddata->totalMass = hdata->totalMass;
  ddata->totalMassi = hdata->totalMassi;
  ddata->NFreedom = hdata->NFreedom;
  
#ifndef COORD_IN_ONE_VEC
  cudaMemcpy (ddata->coordx, hdata->coordx, sizef, cudaMemcpyDeviceToDevice);
  cudaMemcpy (ddata->coordy, hdata->coordy, sizef, cudaMemcpyDeviceToDevice);
  cudaMemcpy (ddata->coordz, hdata->coordz, sizef, cudaMemcpyDeviceToDevice);
#else
  cudaMemcpy (ddata->coord, hdata->coord, sizecoord, cudaMemcpyDeviceToDevice);
#endif
  checkCUDAError ("cpyDeviceMDDataToDevice coord");

  cudaMemcpy (ddata->coordNoix, hdata->coordNoix, sizei, cudaMemcpyDeviceToDevice);
  cudaMemcpy (ddata->coordNoiy, hdata->coordNoiy, sizei, cudaMemcpyDeviceToDevice);
  cudaMemcpy (ddata->coordNoiz, hdata->coordNoiz, sizei, cudaMemcpyDeviceToDevice);
  checkCUDAError ("cpyDeviceMDDataToDevice coordNoi");
  
  cudaMemcpy (ddata->velox, hdata->velox, sizef, cudaMemcpyDeviceToDevice);
  cudaMemcpy (ddata->veloy, hdata->veloy, sizef, cudaMemcpyDeviceToDevice);
  cudaMemcpy (ddata->veloz, hdata->veloz, sizef, cudaMemcpyDeviceToDevice);
  checkCUDAError ("cpyDeviceMDDataToDevice velo");

  cudaMemcpy (ddata->forcx, hdata->forcx, sizef, cudaMemcpyDeviceToDevice);
  cudaMemcpy (ddata->forcy, hdata->forcy, sizef, cudaMemcpyDeviceToDevice);
  cudaMemcpy (ddata->forcz, hdata->forcz, sizef, cudaMemcpyDeviceToDevice);
  checkCUDAError ("cpyDeviceMDDataToDevice forc");

  cudaMemcpy (ddata->type, hdata->type, ddata->numAtom * sizeof(TypeType),
	      cudaMemcpyDeviceToDevice);
  cudaMemcpy (ddata->mass, hdata->mass, sizef, cudaMemcpyDeviceToDevice);
  cudaMemcpy (ddata->massi, hdata->massi, sizef, cudaMemcpyDeviceToDevice);
  cudaMemcpy (ddata->charge, hdata->charge, sizef, cudaMemcpyDeviceToDevice);   
  checkCUDAError ("cpyDeviceMDDataToDevice other");
}


#ifdef COORD_IN_ONE_VEC
__global__ void deviceCpyTypeToCoordW (CoordType * coord,
				       const TypeType * type,
				       const IndexType N)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  if (ii < N) coord[ii].w = type[ii];
}
#endif

__host__ void initDeviceMDData (const HostMDData * hdata,
				DeviceMDData * ddata)
{
  size_t sizef = hdata->numMem * sizeof(ScalorType);
  size_t sizei = hdata->numMem * sizeof(IntScalorType);
#ifdef COORD_IN_ONE_VEC
  size_t sizecoord =hdata->numMem * sizeof(CoordType);
#endif
  ddata->numAtom = hdata->numAtom;
  ddata->numMem = hdata->numMem;
  ddata->totalMass = hdata->totalMass;
  ddata->totalMassi = hdata->totalMassi;
  ddata->NFreedom = hdata->NFreedom;

#ifndef COORD_IN_ONE_VEC
  cudaMalloc ((void**) &ddata->coordx, sizef);
  cudaMalloc ((void**) &ddata->coordy, sizef);
  cudaMalloc ((void**) &ddata->coordz, sizef);
#else
  cudaMalloc ((void**) &ddata->coord, sizecoord);
#endif
  checkCUDAError ("initDeviceMDData coord");

  cudaMalloc ((void**) &ddata->velox, sizef);
  cudaMalloc ((void**) &ddata->veloy, sizef);
  cudaMalloc ((void**) &ddata->veloz, sizef);
  checkCUDAError ("initDeviceMDData velo");

  cudaMalloc ((void**) &ddata->forcx, sizef);
  cudaMalloc ((void**) &ddata->forcy, sizef);
  cudaMalloc ((void**) &ddata->forcz, sizef);
  checkCUDAError ("initDeviceMDData forc");

  cudaMalloc ((void**) &ddata->coordNoix, sizei);
  cudaMalloc ((void**) &ddata->coordNoiy, sizei);
  cudaMalloc ((void**) &ddata->coordNoiz, sizei);
  checkCUDAError ("initDeviceMDData coordNoi");

  cudaMalloc ((void**) &ddata->type, hdata->numMem * sizeof(TypeType));
  cudaMalloc ((void**) &ddata->mass, sizef);
  cudaMalloc ((void**) &ddata->massi, sizef);
  cudaMalloc ((void**) &ddata->charge, sizef);
  checkCUDAError ("initDeviceMDData other");

  ddata->malloced = true;
  
  /* bindTextureOnDeviceMDData (ddata); */
  cpyHostMDDataToDevice (hdata, ddata);
  
#ifdef COORD_IN_ONE_VEC
  IndexType nob ;
  dim3 myBlockDim;
  myBlockDim.x = 64;
  if (ddata->numAtom % myBlockDim.x == 0){
    nob = ddata->numAtom / myBlockDim.x;
  } else {
    nob = ddata->numAtom / myBlockDim.x + 1;
  }
  dim3 atomGridDim = toGridDim (nob);
  deviceCpyTypeToCoordW<<<atomGridDim, myBlockDim>>> (
      ddata->coord, ddata->type, ddata->numAtom);
#endif
}


__host__ void destroyDeviceMDData (DeviceMDData * ddata)
{
  if (ddata->malloced){
#ifndef COORD_IN_ONE_VEC
    cudaFree (ddata->coordx);
    cudaFree (ddata->coordy);
    cudaFree (ddata->coordz);
#else
    cudaFree (ddata->coord);
#endif
    
    cudaFree (ddata->coordNoix);
    cudaFree (ddata->coordNoiy);
    cudaFree (ddata->coordNoiz);
  
    cudaFree (ddata->velox);
    cudaFree (ddata->veloy);
    cudaFree (ddata->veloz);
  
    cudaFree (ddata->forcx);
    cudaFree (ddata->forcy);
    cudaFree (ddata->forcz);

    cudaFree (ddata->type);
    cudaFree (ddata->mass);
    cudaFree (ddata->massi);
    cudaFree (ddata->charge);
    ddata->malloced = false;
  }
}


////////////////////////////////////////////////////////////
// implementation
// ////////////////////////////////////////////////////////

__device__ void
cpyDeviceMDDataElement (const DeviceMDData * ddata1,
			const IndexType indx1,
			DeviceMDData * ddata2,
			const IndexType indx2)
{
#ifndef COORD_IN_ONE_VEC
  ddata2->coordx[indx2] = ddata1->coordx[indx1];
  ddata2->coordy[indx2] = ddata1->coordy[indx1];
  ddata2->coordz[indx2] = ddata1->coordz[indx1];
#else
  (ddata2->coord[indx2]) = (ddata1->coord[indx1]);
#endif
  
  ddata2->coordNoix[indx2] = ddata1->coordNoix[indx1];
  ddata2->coordNoiy[indx2] = ddata1->coordNoiy[indx1];
  ddata2->coordNoiz[indx2] = ddata1->coordNoiz[indx1];

  ddata2->velox[indx2] = ddata1->velox[indx1];
  ddata2->veloy[indx2] = ddata1->veloy[indx1];
  ddata2->veloz[indx2] = ddata1->veloz[indx1];

  ddata2->forcx[indx2] = ddata1->forcx[indx1];
  ddata2->forcy[indx2] = ddata1->forcy[indx1];
  ddata2->forcz[indx2] = ddata1->forcz[indx1];

  ddata2->type[indx2] = ddata1->type[indx1];
  ddata2->mass[indx2] = ddata1->mass[indx1];
  ddata2->massi[indx2] = ddata1->massi[indx1];
  ddata2->charge[indx2] = ddata1->charge[indx1];

  // IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  // IndexType tid = threadIdx.x;

  // if (bid + tid == 0){
  //   ddata2->numAtom = ddata1->numAtom;
  //   ddata2->numMem  = ddata1->numMem;
  //   ddata2->totalMass = ddata1->totalMass;
  //   ddata2->totalMassi = ddata1->totalMassi;
  // }
}
