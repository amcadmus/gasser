#define DEVICE_CODE
#include "Parallel_MDData.h"
#include "Parallel_MDData_device.h"

#include "compile_error_mixcode.h"

Parallel::DeviceMDData::
DeviceMDData ()
    : numAtom_ (0), memSize_(0), malloced(false)
{
}

Parallel::DeviceMDData::
~DeviceMDData ()
{
  clearAll();
}

void Parallel::DeviceMDData::
mallocAll (const IndexType & memSize__)
{
  if (malloced){
    clearAll ();
  }
  if (memSize__ == 0) return;

  memSize_ = memSize__;
  
  size_t sizef = memSize_ * sizeof(ScalorType);
  size_t sizei = memSize_ * sizeof(IntScalorType);
  size_t sizecoord =memSize_ * sizeof(CoordType);
  size_t sizeIdx = memSize_ * sizeof(IndexType);
  
  cudaMalloc ((void**) &coord, sizecoord);
  checkCUDAError ("initDeviceMDData coord");

  cudaMalloc ((void**) &velox, sizef);
  cudaMalloc ((void**) &veloy, sizef);
  cudaMalloc ((void**) &veloz, sizef);
  checkCUDAError ("initDeviceMDData velo");

  cudaMalloc ((void**) &forcx, sizef);
  cudaMalloc ((void**) &forcy, sizef);
  cudaMalloc ((void**) &forcz, sizef);
  checkCUDAError ("initDeviceMDData forc");

  cudaMalloc ((void**) &coordNoix, sizei);
  cudaMalloc ((void**) &coordNoiy, sizei);
  cudaMalloc ((void**) &coordNoiz, sizei);
  checkCUDAError ("initDeviceMDData coordNoi");

  cudaMalloc ((void**) &globalIndex, sizeIdx);
  cudaMalloc ((void**) &type, memSize_ * sizeof(TypeType));
  cudaMalloc ((void**) &mass, sizef);
  cudaMalloc ((void**) &charge, sizef);
  checkCUDAError ("initDeviceMDData top Property");

  malloced = true;
}

void Parallel::DeviceMDData::
clearAll ()
{
  if (malloced){
    cudaFree (coord);
    
    cudaFree (coordNoix);
    cudaFree (coordNoiy);
    cudaFree (coordNoiz);
  
    cudaFree (velox);
    cudaFree (veloy);
    cudaFree (veloz);
  
    cudaFree (forcx);
    cudaFree (forcy);
    cudaFree (forcz);

    cudaFree (globalIndex);
    cudaFree (type);
    cudaFree (mass);
    cudaFree (charge);
    malloced = false;
  }
}


void Parallel::DeviceMDData::
copyFromHost (const HostMDData & hdata)
{
  if (memSize_ < hdata.numAtom()){
    clearAll();
    mallocAll(hdata.numAtom());
  }
  numAtom_ = hdata.numAtom();
  setGlobalBox (hdata.getGlobalBox());
		
  size_t sizef = numAtom_ * sizeof(ScalorType);
  size_t sizei = numAtom_ * sizeof(IntScalorType);
  size_t sizecoord = numAtom_ * sizeof(CoordType);
  size_t sizeIdx = memSize_ * sizeof(IndexType);
  
  cudaMemcpy (coord, hdata.coord, sizecoord, cudaMemcpyHostToDevice);
  checkCUDAError ("cpyHostMDDataToDevice coord");

  cudaMemcpy (coordNoix, hdata.coordNoix, sizei, cudaMemcpyHostToDevice);
  cudaMemcpy (coordNoiy, hdata.coordNoiy, sizei, cudaMemcpyHostToDevice);
  cudaMemcpy (coordNoiz, hdata.coordNoiz, sizei, cudaMemcpyHostToDevice);
  checkCUDAError ("cpyHostMDDataToDevice coordNoi");
  
  cudaMemcpy (velox, hdata.velox, sizef, cudaMemcpyHostToDevice);
  cudaMemcpy (veloy, hdata.veloy, sizef, cudaMemcpyHostToDevice);
  cudaMemcpy (veloz, hdata.veloz, sizef, cudaMemcpyHostToDevice);
  checkCUDAError ("cpyHostMDDataToDevice velo");

  cudaMemcpy (forcx, hdata.forcx, sizef, cudaMemcpyHostToDevice);
  cudaMemcpy (forcy, hdata.forcy, sizef, cudaMemcpyHostToDevice);
  cudaMemcpy (forcz, hdata.forcz, sizef, cudaMemcpyHostToDevice);
  checkCUDAError ("cpyHostMDDataToDevice forc");

  cudaMemcpy (globalIndex, hdata.globalIndex, sizeIdx, cudaMemcpyHostToDevice);
  cudaMemcpy (type, hdata.type, numAtom_ * sizeof(TypeType), cudaMemcpyHostToDevice);
  cudaMemcpy (mass, hdata.mass, sizef, cudaMemcpyHostToDevice);
  cudaMemcpy (charge, hdata.charge, sizef, cudaMemcpyHostToDevice);   
  checkCUDAError ("cpyHostMDDataToDevice other");
}


void Parallel::DeviceMDData::
copyToHost (HostMDData & hdata) const
{
  if (hdata.memSize() < numAtom_){
    hdata.reallocAll (numAtom_);
  }
  hdata.numAtom() = numAtom_;
  hdata.setGlobalBox (globalBox);
  
  size_t sizef = numAtom_ * sizeof(ScalorType);
  size_t sizei = numAtom_ * sizeof(IntScalorType);
  size_t sizecoord = numAtom_ * sizeof(CoordType);
  size_t sizeIdx = memSize_ * sizeof(IndexType);
  
  cudaMemcpy (hdata.coord, coord, sizecoord, cudaMemcpyDeviceToHost);
  checkCUDAError ("cpyDeviceMDDataToHost coord");
  
  cudaMemcpy (hdata.coordNoix, coordNoix, sizei, cudaMemcpyDeviceToHost);
  cudaMemcpy (hdata.coordNoiy, coordNoiy, sizei, cudaMemcpyDeviceToHost);
  cudaMemcpy (hdata.coordNoiz, coordNoiz, sizei, cudaMemcpyDeviceToHost);
  checkCUDAError ("cpyDeviceMDDataToHost coordNoi");
  
  cudaMemcpy (hdata.velox, velox, sizef, cudaMemcpyDeviceToHost);
  cudaMemcpy (hdata.veloy, veloy, sizef, cudaMemcpyDeviceToHost);
  cudaMemcpy (hdata.veloz, veloz, sizef, cudaMemcpyDeviceToHost);
  checkCUDAError ("cpyDeviceMDDataToHost velo");

  cudaMemcpy (hdata.forcx, forcx, sizef, cudaMemcpyDeviceToHost);
  cudaMemcpy (hdata.forcy, forcy, sizef, cudaMemcpyDeviceToHost);
  cudaMemcpy (hdata.forcz, forcz, sizef, cudaMemcpyDeviceToHost);
  checkCUDAError ("cpyDeviceMDDataToHost forc");

  cudaMemcpy (hdata.globalIndex, globalIndex, sizeIdx, cudaMemcpyDeviceToHost);
  cudaMemcpy (hdata.type, type, numAtom_ * sizeof(TypeType), cudaMemcpyDeviceToHost);
  cudaMemcpy (hdata.mass, mass, sizef, cudaMemcpyDeviceToHost);
  cudaMemcpy (hdata.charge, charge, sizef, cudaMemcpyDeviceToHost);  
  checkCUDAError ("cpyDeviceMDDataToHost other");
}


void Parallel::DeviceMDData::
copyFromDevice (const DeviceMDData & ddata)
{
  if (memSize_ < ddata.memSize_){
    clearAll();
    mallocAll (ddata.memSize());
  }
  numAtom_ = ddata.numAtom();
  setGlobalBox (ddata.getGlobalBox());
		
  size_t sizef = numAtom_ * sizeof(ScalorType);
  size_t sizei = numAtom_ * sizeof(IntScalorType);
  size_t sizecoord = numAtom_ * sizeof(CoordType);
  size_t sizeIdx = memSize_ * sizeof(IndexType);
  
  cudaMemcpy (coord, ddata.coord, sizecoord, cudaMemcpyDeviceToDevice);
  checkCUDAError ("cpyDeviceMDDataToDevice coord");

  cudaMemcpy (coordNoix, ddata.coordNoix, sizei, cudaMemcpyDeviceToDevice);
  cudaMemcpy (coordNoiy, ddata.coordNoiy, sizei, cudaMemcpyDeviceToDevice);
  cudaMemcpy (coordNoiz, ddata.coordNoiz, sizei, cudaMemcpyDeviceToDevice);
  checkCUDAError ("cpyDeviceMDDataToDevice coordNoi");
  
  cudaMemcpy (velox, ddata.velox, sizef, cudaMemcpyDeviceToDevice);
  cudaMemcpy (veloy, ddata.veloy, sizef, cudaMemcpyDeviceToDevice);
  cudaMemcpy (veloz, ddata.veloz, sizef, cudaMemcpyDeviceToDevice);
  checkCUDAError ("cpyDeviceMDDataToDevice velo");

  cudaMemcpy (forcx, ddata.forcx, sizef, cudaMemcpyDeviceToDevice);
  cudaMemcpy (forcy, ddata.forcy, sizef, cudaMemcpyDeviceToDevice);
  cudaMemcpy (forcz, ddata.forcz, sizef, cudaMemcpyDeviceToDevice);
  checkCUDAError ("cpyDeviceMDDataToDevice forc");

  cudaMemcpy (globalIndex, ddata.globalIndex, sizeIdx, cudaMemcpyDeviceToDevice);
  cudaMemcpy (type, ddata.type, numAtom_ * sizeof(TypeType), cudaMemcpyDeviceToDevice);
  cudaMemcpy (mass, ddata.mass, sizef, cudaMemcpyDeviceToDevice);
  cudaMemcpy (charge, ddata.charge, sizef, cudaMemcpyDeviceToDevice);   
  checkCUDAError ("cpyDeviceMDDataToDevice other");
}


Parallel::DeviceMDData::
DeviceMDData (const DeviceMDData & ddata)
    : numAtom_ (0), memSize_(0), malloced(false)
{
  copyFromDevice (ddata);
}


void Parallel::DeviceMDData::
initZero ()
{
  Parallel::CudaGlobal::initZeroDeviceData
      <<<memSize_ / DefaultNThreadPerBlock + 1, DefaultNThreadPerBlock>>>
      (memSize_,
       coord,
       coordNoix,
       coordNoiy,
       coordNoiz,
       velox,
       veloy,
       veloz,
       globalIndex,
       type,
       mass,
       charge);
  checkCUDAError ("DeviceMDData::initZero initZeroDeviceData");
}


__global__ void Parallel::CudaGlobal::
initZeroDeviceData(const IndexType num,
		   CoordType  * coord,
		   IntScalorType * coordNoix,
		   IntScalorType * coordNoiy,
		   IntScalorType * coordNoiz,
		   ScalorType * velox,
		   ScalorType * veloy,
		   ScalorType * veloz,
		   IndexType  * globalIndex,
		   TypeType   * type,
		   ScalorType * mass,
		   ScalorType * charge)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;

  if (ii < num){
    CoordType tmp;
    tmp.x = 0;
    tmp.y = 0;
    tmp.z = 0;
    tmp.w = MaxIndexValue;
    coord[ii] = tmp;
    coordNoix[ii] = coordNoiy[ii] = coordNoiz[ii] = 0;
    veloz[ii] = veloy[ii] = veloz[ii] = 0.f;
    globalIndex[ii] = MaxIndexValue;
    type[ii] = 0;
    mass[ii] = 0;
    charge[ii] = 0;
  }
}



  
