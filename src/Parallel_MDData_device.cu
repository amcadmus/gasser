#define DEVICE_CODE

#include "Parallel_MDData.h"

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

  size_t sizef = numAtom_ * sizeof(ScalorType);
  size_t sizei = numAtom_ * sizeof(IntScalorType);
  size_t sizecoord = numAtom_ * sizeof(CoordType);
  
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

  size_t sizef = numAtom_ * sizeof(ScalorType);
  size_t sizei = numAtom_ * sizeof(IntScalorType);
  size_t sizecoord = numAtom_ * sizeof(CoordType);
  
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

  cudaMemcpy (hdata.type, type, numAtom_ * sizeof(TypeType), cudaMemcpyDeviceToHost);
  cudaMemcpy (hdata.mass, mass, sizef, cudaMemcpyDeviceToHost);
  cudaMemcpy (hdata.charge, charge, sizef, cudaMemcpyDeviceToHost);  
  checkCUDAError ("cpyDeviceMDDataToHost other");
}

