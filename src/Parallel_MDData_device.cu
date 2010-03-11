#define DEVICE_CODE
#include "Parallel_MDData.h"
#include "Parallel_MDData_device.h"

#include "compile_error_mixcode.h"

Parallel::DeviceMDData::
DeviceMDData ()
    : numData_ (0), memSize_(0), malloced(false)
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
  size_t sizet = memSize_ * sizeof(TypeType);
  
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
  cudaMalloc ((void**) &type, sizet);
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
copyFromHost (const HostMDData & hdata,
	      const MDDataItemMask_t mask)
{
  if (memSize_ < hdata.numData()){
    clearAll();
    mallocAll (hdata.numData());
  }
  numData_ = hdata.numData();
  setGlobalBox (hdata.getGlobalBox());
		
  size_t sizef = numData_ * sizeof(ScalorType);
  size_t sizei = numData_ * sizeof(IntScalorType);
  size_t sizecoord = numData_ * sizeof(CoordType);
  size_t sizeIdx = numData_ * sizeof(IndexType);
  size_t sizet = numData_ * sizeof(TypeType);

  if (mask & MDDataItemMask_Coordinate){
    cudaMemcpy (coord, hdata.coord, sizecoord, cudaMemcpyHostToDevice);
    checkCUDAError ("cpyHostMDDataToDevice coord");
  }

  if (mask & MDDataItemMask_CoordinateNoi){
    cudaMemcpy (coordNoix, hdata.coordNoix, sizei, cudaMemcpyHostToDevice);
    cudaMemcpy (coordNoiy, hdata.coordNoiy, sizei, cudaMemcpyHostToDevice);
    cudaMemcpy (coordNoiz, hdata.coordNoiz, sizei, cudaMemcpyHostToDevice);
    checkCUDAError ("cpyHostMDDataToDevice coordNoi");
  }

  if (mask & MDDataItemMask_Velocity){
    cudaMemcpy (velox, hdata.velox, sizef, cudaMemcpyHostToDevice);
    cudaMemcpy (veloy, hdata.veloy, sizef, cudaMemcpyHostToDevice);
    cudaMemcpy (veloz, hdata.veloz, sizef, cudaMemcpyHostToDevice);
    checkCUDAError ("cpyHostMDDataToDevice velo");
  }

  if (mask & MDDataItemMask_Force){
    cudaMemcpy (forcx, hdata.forcx, sizef, cudaMemcpyHostToDevice);
    cudaMemcpy (forcy, hdata.forcy, sizef, cudaMemcpyHostToDevice);
    cudaMemcpy (forcz, hdata.forcz, sizef, cudaMemcpyHostToDevice);
    checkCUDAError ("cpyHostMDDataToDevice forc");
  }

  if (mask & MDDataItemMask_GlobalIndex){
    cudaMemcpy (globalIndex, hdata.globalIndex, sizeIdx, cudaMemcpyHostToDevice);
    checkCUDAError ("cpyHostMDDataToDevice globalIndex");
  }
  if (mask & MDDataItemMask_Type){
    cudaMemcpy (type, hdata.type, sizet, cudaMemcpyHostToDevice);
    checkCUDAError ("cpyHostMDDataToDevice type");
  }
  if (mask & MDDataItemMask_Mass){
    cudaMemcpy (mass, hdata.mass, sizef, cudaMemcpyHostToDevice);
    checkCUDAError ("cpyHostMDDataToDevice mass");
  }
  if (mask & MDDataItemMask_Charge){
    cudaMemcpy (charge, hdata.charge, sizef, cudaMemcpyHostToDevice);   
    checkCUDAError ("cpyHostMDDataToDevice charge");
  }
}


void Parallel::DeviceMDData::
copyToHost (HostMDData & hdata,
	    const MDDataItemMask_t mask) const
{
  if (hdata.memSize() < numData_){
    hdata.reallocAll (numData_);
  }
  hdata.numData_ = numData_;
  hdata.setGlobalBox (globalBox);
  
  size_t sizef = numData_ * sizeof(ScalorType);
  size_t sizei = numData_ * sizeof(IntScalorType);
  size_t sizecoord = numData_ * sizeof(CoordType);
  size_t sizeIdx = numData_ * sizeof(IndexType);
  size_t sizet = numData_ * sizeof(TypeType);

  if (mask & MDDataItemMask_Coordinate){
    cudaMemcpy (hdata.coord, coord, sizecoord, cudaMemcpyDeviceToHost);
    checkCUDAError ("cpyDeviceMDDataToHost coord");
  }

  if (mask & MDDataItemMask_CoordinateNoi){
    cudaMemcpy (hdata.coordNoix, coordNoix, sizei, cudaMemcpyDeviceToHost);
    cudaMemcpy (hdata.coordNoiy, coordNoiy, sizei, cudaMemcpyDeviceToHost);
    cudaMemcpy (hdata.coordNoiz, coordNoiz, sizei, cudaMemcpyDeviceToHost);
    checkCUDAError ("cpyDeviceMDDataToHost coordNoi");
  }

  if (mask & MDDataItemMask_Velocity){
    cudaMemcpy (hdata.velox, velox, sizef, cudaMemcpyDeviceToHost);
    cudaMemcpy (hdata.veloy, veloy, sizef, cudaMemcpyDeviceToHost);
    cudaMemcpy (hdata.veloz, veloz, sizef, cudaMemcpyDeviceToHost);
    checkCUDAError ("cpyDeviceMDDataToHost velo");
  }

  if (mask & MDDataItemMask_Force){
    cudaMemcpy (hdata.forcx, forcx, sizef, cudaMemcpyDeviceToHost);
    cudaMemcpy (hdata.forcy, forcy, sizef, cudaMemcpyDeviceToHost);
    cudaMemcpy (hdata.forcz, forcz, sizef, cudaMemcpyDeviceToHost);
    checkCUDAError ("cpyDeviceMDDataToHost forc");
  }
  
  if (mask & MDDataItemMask_GlobalIndex){
    cudaMemcpy (hdata.globalIndex, globalIndex, sizeIdx, cudaMemcpyDeviceToHost);
    checkCUDAError ("cpyDeviceMDDataToHost globalIndex");
  }
  if (mask & MDDataItemMask_Type){
    cudaMemcpy (hdata.type, type, sizet, cudaMemcpyDeviceToHost);
    checkCUDAError ("cpyDeviceMDDataToHost type");
  }
  if (mask & MDDataItemMask_Mass){
    cudaMemcpy (hdata.mass, mass, sizef, cudaMemcpyDeviceToHost);
    checkCUDAError ("cpyDeviceMDDataToHost mass");
  }
  if (mask & MDDataItemMask_Charge){
    cudaMemcpy (hdata.charge, charge, sizef, cudaMemcpyDeviceToHost);  
    checkCUDAError ("cpyDeviceMDDataToHost charge");
  }
}


void Parallel::DeviceMDData::
copyFromDevice (const DeviceMDData & ddata,
		const MDDataItemMask_t mask)
{
  if (numData_ < ddata.memSize()){
    clearAll();
    mallocAll (ddata.memSize());
  }
  numData_ = ddata.numData();
  setGlobalBox (ddata.getGlobalBox());
		
  size_t sizef = numData_ * sizeof(ScalorType);
  size_t sizei = numData_ * sizeof(IntScalorType);
  size_t sizecoord = numData_ * sizeof(CoordType);
  size_t sizeIdx = numData_ * sizeof(IndexType);
  size_t sizet = numData_ * sizeof(TypeType);

  if (mask & MDDataItemMask_Coordinate){
    cudaMemcpy (coord, ddata.coord, sizecoord, cudaMemcpyDeviceToDevice);
    checkCUDAError ("cpyDeviceMDDataToDevice coord");
  }

  if (mask & MDDataItemMask_CoordinateNoi){
    cudaMemcpy (coordNoix, ddata.coordNoix, sizei, cudaMemcpyDeviceToDevice);
    cudaMemcpy (coordNoiy, ddata.coordNoiy, sizei, cudaMemcpyDeviceToDevice);
    cudaMemcpy (coordNoiz, ddata.coordNoiz, sizei, cudaMemcpyDeviceToDevice);
    checkCUDAError ("cpyDeviceMDDataToDevice coordNoi");
  }

  if (mask & MDDataItemMask_Velocity){
    cudaMemcpy (velox, ddata.velox, sizef, cudaMemcpyDeviceToDevice);
    cudaMemcpy (veloy, ddata.veloy, sizef, cudaMemcpyDeviceToDevice);
    cudaMemcpy (veloz, ddata.veloz, sizef, cudaMemcpyDeviceToDevice);
    checkCUDAError ("cpyDeviceMDDataToDevice velo");
  }

  if (mask & MDDataItemMask_Force){
    cudaMemcpy (forcx, ddata.forcx, sizef, cudaMemcpyDeviceToDevice);
    cudaMemcpy (forcy, ddata.forcy, sizef, cudaMemcpyDeviceToDevice);
    cudaMemcpy (forcz, ddata.forcz, sizef, cudaMemcpyDeviceToDevice);
    checkCUDAError ("cpyDeviceMDDataToDevice forc");
  }

  if (mask & MDDataItemMask_GlobalIndex){
    cudaMemcpy (globalIndex, ddata.globalIndex, sizeIdx, cudaMemcpyDeviceToDevice);
    checkCUDAError ("cpyDeviceMDDataToDevice globalIndex");
  }
  if (mask & MDDataItemMask_Type){
    cudaMemcpy (type, ddata.type, sizet, cudaMemcpyDeviceToDevice);
    checkCUDAError ("cpyDeviceMDDataToDevice globalIndex");
  }
  if (mask & MDDataItemMask_Mass){
    cudaMemcpy (mass, ddata.mass, sizef, cudaMemcpyDeviceToDevice);
    checkCUDAError ("cpyDeviceMDDataToDevice globalIndex");
  }
  if (mask & MDDataItemMask_Charge){
    cudaMemcpy (charge, ddata.charge, sizef, cudaMemcpyDeviceToDevice);     
    checkCUDAError ("cpyDeviceMDDataToDevice globalIndex");
  }
}


Parallel::DeviceMDData::
DeviceMDData (const DeviceMDData & ddata)
    : numData_ (0), memSize_(0), malloced(false)
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
       forcx,
       forcy,
       forcz,
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
		   ScalorType * forcx,
		   ScalorType * forcy,
		   ScalorType * forcz,
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
    forcz[ii] = forcy[ii] = forcz[ii] = 0.f;
    globalIndex[ii] = MaxIndexValue;
    type[ii] = 0;
    mass[ii] = 0;
    charge[ii] = 0;
  }
}



  
