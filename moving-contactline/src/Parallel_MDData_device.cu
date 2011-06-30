#define DEVICE_CODE
#include "Parallel_MDData.h"
#include "Parallel_MDData_device.h"
#include "Parallel_Interface.h"
#include "Parallel_Auxiliary.h"

#include "compile_error_mixcode.h"

Parallel::DeviceMDData::
DeviceMDData ()
    : _numData(0),
      _memSize(0),
      maxNumBond(0),
      maxNumAngle(0),
      maxNumDihedral(0)
{
}

Parallel::DeviceMDData::
~DeviceMDData ()
{
  clear();
}

void Parallel::DeviceMDData::
clear ()
{
  if (memSize() != 0){    
    cudaFree (coord);
    cudaFree (coordNoi);
    checkCUDAError ("DeviceMDData::clear");
  
    cudaFree (velox);
    cudaFree (veloy);
    cudaFree (veloz);
    checkCUDAError ("DeviceMDData::clear");
  
    cudaFree (forcx);
    cudaFree (forcy);
    cudaFree (forcz);
    checkCUDAError ("DeviceMDData::clear");

    cudaFree (globalIndex);
    cudaFree (type);
    cudaFree (mass);
    cudaFree (charge);

    _memSize = 0;
    _numData = 0;
    checkCUDAError ("DeviceMDData::clear");
    
    if (maxNumBond != 0){
      cudaFree (numBond);
      cudaFree (bondNeighbor_globalIndex);
      cudaFree (bondIndex);
      maxNumBond = 0;
    }
    checkCUDAError ("DeviceMDData::clearBondTop");
    if (maxNumAngle != 0){
      cudaFree (numAngle);
      cudaFree (angleNeighbor_globalIndex);
      cudaFree (angleIndex);
      cudaFree (anglePosi);
      maxNumAngle = 0;
    }
    checkCUDAError ("DeviceMDData::clearBondTop");
    if (maxNumDihedral != 0){
      cudaFree (numDihedral);
      cudaFree (dihedralNeighbor_globalIndex);
      cudaFree (dihedralIndex);
      cudaFree (dihedralPosi);
      maxNumDihedral = 0;
    }
    checkCUDAError ("DeviceMDData::clearBondTop");
  }
}

void Parallel::DeviceMDData::
fillZero ()
{
  CoordType coord0;
  coord0.x = coord0.y = coord0.z = coord0.w = 0.f;
  CoordNoiType coordNoi0;
  coordNoi0.x = coordNoi0.y = coordNoi0.z = 0;
  IndexType numThreadsInCell = Parallel::Interface::numThreadsInCell();
  IndexType numBlock = memSize() / numThreadsInCell + 1;
  dim3 gridDim (toGridDim(numBlock));

  using namespace Parallel::Auxiliary;
  
  setValue <<<gridDim, numThreadsInCell>>> (coord, memSize(), coord0);
  setValue <<<gridDim, numThreadsInCell>>> (coordNoi, memSize(), coordNoi0);
  setValue <<<gridDim, numThreadsInCell>>> (velox, memSize(), 0.f);
  setValue <<<gridDim, numThreadsInCell>>> (veloy, memSize(), 0.f);
  setValue <<<gridDim, numThreadsInCell>>> (veloz, memSize(), 0.f);
  setValue <<<gridDim, numThreadsInCell>>> (forcx, memSize(), 0.f);
  setValue <<<gridDim, numThreadsInCell>>> (forcy, memSize(), 0.f);
  setValue <<<gridDim, numThreadsInCell>>> (forcz, memSize(), 0.f);
  setValue <<<gridDim, numThreadsInCell>>> (globalIndex, memSize(), MaxIndexValue);
  setValue <<<gridDim, numThreadsInCell>>> (type, memSize(), TypeType(0));
  setValue <<<gridDim, numThreadsInCell>>> (mass, memSize(), 0.f);
  setValue <<<gridDim, numThreadsInCell>>> (charge, memSize(), 0.f);
  
  if (maxNumBond != 0){
    setValue<<<gridDim, numThreadsInCell>>> (numBond, bondTopStride(), IndexType(0));
    numBlock = bondTopStride() * maxNumBond / numThreadsInCell + 1;
    dim3 gridDim1 (toGridDim(numBlock));
    setValue<<<gridDim1, numThreadsInCell>>> (
	bondIndex, bondTopStride() * maxNumBond, IndexType(0));
    setValue<<<gridDim1, numThreadsInCell>>> (
	bondNeighbor_globalIndex, bondTopStride() * maxNumBond, MaxIndexValue);
  }
  if (maxNumAngle != 0){
    setValue<<<gridDim, numThreadsInCell>>> (numAngle, bondTopStride(), IndexType(0));
    numBlock = bondTopStride() * maxNumAngle / numThreadsInCell + 1;
    dim3 gridDim1 (toGridDim(numBlock));
    numBlock = bondTopStride() * maxNumAngle * 2 / numThreadsInCell + 1;
    dim3 gridDim2 (toGridDim(numBlock));
    setValue<<<gridDim1, numThreadsInCell>>> (
	angleIndex, bondTopStride() * maxNumAngle, IndexType(0));
    setValue<<<gridDim1, numThreadsInCell>>> (
	anglePosi,  bondTopStride() * maxNumAngle, IndexType(0));
    setValue<<<gridDim2, numThreadsInCell>>> (
	angleNeighbor_globalIndex, bondTopStride() * maxNumAngle * 2, MaxIndexValue);
  }
  if (maxNumDihedral != 0){
    setValue<<<gridDim, numThreadsInCell>>> (numDihedral, bondTopStride(), IndexType(0));
    numBlock = bondTopStride() * maxNumDihedral / numThreadsInCell + 1;
    dim3 gridDim1 (toGridDim(numBlock));
    numBlock = bondTopStride() * maxNumDihedral * 3 / numThreadsInCell + 1;
    dim3 gridDim2 (toGridDim(numBlock));
    setValue<<<gridDim1, numThreadsInCell>>> (
	dihedralIndex, bondTopStride() * maxNumDihedral, IndexType(0));
    setValue<<<gridDim1, numThreadsInCell>>> (
	dihedralPosi,  bondTopStride() * maxNumDihedral, IndexType(0));
    setValue<<<gridDim2, numThreadsInCell>>> (
	dihedralNeighbor_globalIndex,
	bondTopStride() * maxNumDihedral * 3,
	MaxIndexValue);
  }
}


void Parallel::DeviceMDData::
easyMalloc (const IndexType memSize_,
	    const IndexType maxNumBond_,
	    const IndexType maxNumAngle_,
	    const IndexType maxNumDihedral_)
{
  clear ();

  _memSize = memSize_;
  maxNumBond = maxNumBond_;
  maxNumAngle = maxNumAngle_;
  maxNumDihedral = maxNumDihedral_;

  if (_memSize == 0) return;
  
  size_t sizef = memSize() * sizeof(ScalorType);
  size_t sizecoord =memSize() * sizeof(CoordType);
  size_t sizecoordNoi =memSize() * sizeof(CoordNoiType);
  size_t sizeIdx = memSize() * sizeof(IndexType);
  size_t sizet = memSize() * sizeof(TypeType);
  
  cudaMalloc ((void**) &coord, sizecoord);
  checkCUDAError ("initDeviceMDData coord");

  cudaMalloc ((void**) &coordNoi, sizecoordNoi);
  checkCUDAError ("initDeviceMDData coordNoi");

  cudaMalloc ((void**) &velox, sizef);
  cudaMalloc ((void**) &veloy, sizef);
  cudaMalloc ((void**) &veloz, sizef);
  checkCUDAError ("initDeviceMDData velo");

  cudaMalloc ((void**) &forcx, sizef);
  cudaMalloc ((void**) &forcy, sizef);
  cudaMalloc ((void**) &forcz, sizef);
  checkCUDAError ("initDeviceMDData forc");

  cudaMalloc ((void**) &globalIndex, sizeIdx);
  cudaMalloc ((void**) &type, sizet);
  cudaMalloc ((void**) &mass, sizef);
  cudaMalloc ((void**) &charge, sizef);
  checkCUDAError ("initDeviceMDData top Property");

  size_t size0 = sizeof(IndexType) * memSize();
  if (maxNumBond != 0){
    size_t size1 = size0 * maxNumBond;
    cudaMalloc ((void**)&numBond, size0);
    cudaMalloc ((void**)&bondIndex, size1);
    cudaMalloc ((void**)&bondNeighbor_globalIndex, size1);
    checkCUDAError ("DeviceMDData::easyMallocBondTop, bond");
  }
  if (maxNumAngle != 0){
    size_t size1 = size0 * maxNumAngle;
    size_t size2 = size0 * maxNumAngle * 2;
    cudaMalloc ((void**)&numAngle, size0);
    cudaMalloc ((void**)&angleIndex, size1);
    cudaMalloc ((void**)&anglePosi, size1);
    cudaMalloc ((void**)&angleNeighbor_globalIndex, size2);
    checkCUDAError ("DeviceMDData::easyMallocAngleTop, angle");
  }
  if (maxNumDihedral != 0){
    size_t size1 = size0 * maxNumDihedral;
    size_t size2 = size0 * maxNumDihedral * 3;
    cudaMalloc ((void**)&numDihedral, size0);
    cudaMalloc ((void**)&dihedralIndex, size1);
    cudaMalloc ((void**)&dihedralPosi, size1);
    cudaMalloc ((void**)&dihedralNeighbor_globalIndex, size2);
    checkCUDAError ("DeviceMDData::easyMallocDihedralTop, dihedral");
  }

  fillZero();
}


void Parallel::DeviceMDData::
copyFromHost (const HostMDData & hdata,
	      const MDDataItemMask_t mask)
{
  // if (!mask) return;
  IndexType expectedNumBond(0), expectedNumAngle(0), expectedNumDihedral(0);
  bool copyBond = (mask & MDDataItemMask_Bond);
  bool copyAngle = (mask & MDDataItemMask_Angle);
  bool copyDihedral = (mask & MDDataItemMask_Dihedral);  
  if (copyBond){
    expectedNumBond = hdata.getMaxNumBond();
  }
  if (copyAngle){
    expectedNumAngle = hdata.getMaxNumAngle();
  }
  if (copyDihedral){
    expectedNumDihedral = hdata.getMaxNumDihedral();
  }

  if (memSize() != hdata.memSize() ||
      (copyBond && (maxNumBond != hdata.maxNumBond)) ||
      (copyAngle && (maxNumAngle != hdata.maxNumAngle)) ||
      (copyDihedral && (maxNumDihedral != hdata.maxNumDihedral)) ){
    easyMalloc (hdata.memSize(), expectedNumBond, expectedNumAngle, expectedNumDihedral);
  }
  
  numData() = hdata.numData();
  setGlobalBox (hdata.getGlobalBox());
		
  size_t sizef = hdata.numData() * sizeof(ScalorType);
  size_t sizecoord = hdata.numData() * sizeof(CoordType);
  size_t sizecoordNoi = hdata.numData() * sizeof(CoordNoiType);
  size_t sizeIdx = hdata.numData() * sizeof(IndexType);
  size_t sizet = hdata.numData() * sizeof(TypeType);

  if (mask & MDDataItemMask_Coordinate){
    cudaMemcpy (coord, hdata.coord, sizecoord, cudaMemcpyHostToDevice);
    checkCUDAError ("cpyHostMDDataToDevice coord");
  }

  if (mask & MDDataItemMask_CoordinateNoi){
    cudaMemcpy (coordNoi, hdata.coordNoi, sizecoordNoi, cudaMemcpyHostToDevice);
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

  size_t size0 = sizeof(IndexType) * hdata.memSize();

  if (copyBond && maxNumBond != 0){
    cudaMemcpy (numBond, hdata.numBond, size0, cudaMemcpyHostToDevice);
    size_t size1 = size0 * maxNumBond;
    cudaMemcpy (bondIndex, hdata.bondIndex, size1, cudaMemcpyHostToDevice);
    cudaMemcpy (bondNeighbor_globalIndex,
		hdata.bondNeighbor_globalIndex,
		size1,
		cudaMemcpyHostToDevice);
  }
  if (copyAngle && maxNumAngle != 0){
    cudaMemcpy (numAngle, hdata.numAngle, size0, cudaMemcpyHostToDevice);
    size_t size1 = size0 * maxNumAngle;
    size_t size2 = size0 * maxNumAngle * 2;
    cudaMemcpy (angleIndex, hdata.angleIndex, size1, cudaMemcpyHostToDevice);
    cudaMemcpy (anglePosi, hdata.anglePosi, size1, cudaMemcpyHostToDevice);
    cudaMemcpy (angleNeighbor_globalIndex,
		hdata.angleNeighbor_globalIndex,
		size2,
		cudaMemcpyHostToDevice);
  }
  if (copyDihedral && maxNumDihedral != 0){
    cudaMemcpy (numDihedral, hdata.numDihedral, size0, cudaMemcpyHostToDevice);
    size_t size1 = size0 * maxNumDihedral;
    size_t size2 = size0 * maxNumDihedral * 3;
    cudaMemcpy (dihedralIndex, hdata.dihedralIndex, size1, cudaMemcpyHostToDevice);
    cudaMemcpy (dihedralPosi, hdata.dihedralPosi, size1, cudaMemcpyHostToDevice);
    cudaMemcpy (dihedralNeighbor_globalIndex,
		hdata.dihedralNeighbor_globalIndex,
		size2,
		cudaMemcpyHostToDevice);
  }
}


void Parallel::DeviceMDData::
copyToHost (HostMDData & hdata,
	    const MDDataItemMask_t mask) const
{
  // if (!mask) return;
  IndexType expectedNumBond(0), expectedNumAngle(0), expectedNumDihedral(0);
  bool copyBond = (mask & MDDataItemMask_Bond);
  bool copyAngle = (mask & MDDataItemMask_Angle);
  bool copyDihedral = (mask & MDDataItemMask_Dihedral);  
  if (copyBond){
    expectedNumBond = maxNumBond;
  }
  if (copyAngle){
    expectedNumAngle = maxNumAngle;
  }
  if (copyDihedral){
    expectedNumDihedral = maxNumDihedral;
  }

  if (memSize() != hdata.memSize() ||
      (copyBond && (maxNumBond != hdata.maxNumBond)) ||
      (copyAngle && (maxNumAngle != hdata.maxNumAngle)) ||
      (copyDihedral && (maxNumDihedral != hdata.maxNumDihedral)) ){
    hdata.easyMalloc (memSize(), expectedNumBond, expectedNumAngle, expectedNumDihedral);
  }

  hdata.setGlobalBox (getGlobalBox());
  hdata.numData() = numData();
  
  size_t sizef = numData() * sizeof(ScalorType);
  size_t sizecoord = numData() * sizeof(CoordType);
  size_t sizecoordNoi = numData() * sizeof(CoordNoiType);
  size_t sizeIdx = numData() * sizeof(IndexType);
  size_t sizet = numData() * sizeof(TypeType);

  if (mask & MDDataItemMask_Coordinate){
    cudaMemcpy (hdata.coord, coord, sizecoord, cudaMemcpyDeviceToHost);
    checkCUDAError ("cpyDeviceMDDataToHost coord");
  }

  if (mask & MDDataItemMask_CoordinateNoi){
    cudaMemcpy (hdata.coordNoi, coordNoi, sizecoordNoi, cudaMemcpyDeviceToHost);
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

  size_t size0 = sizeof(IndexType) * memSize();
  if (expectedNumBond != 0){
    cudaMemcpy (hdata.numBond, numBond, size0, cudaMemcpyDeviceToHost);
    size_t size1 = size0 * maxNumBond;
    cudaMemcpy (hdata.bondIndex, bondIndex, size1, cudaMemcpyDeviceToHost);
    cudaMemcpy (hdata.bondNeighbor_globalIndex,
		bondNeighbor_globalIndex,
		size1,
		cudaMemcpyDeviceToHost);
  }
  if (expectedNumAngle != 0){
    cudaMemcpy (hdata.numAngle, numAngle, size0, cudaMemcpyDeviceToHost);
    size_t size1 = size0 * maxNumAngle;
    size_t size2 = size0 * maxNumAngle * 2;
    cudaMemcpy (hdata.angleIndex, angleIndex, size1, cudaMemcpyDeviceToHost);
    cudaMemcpy (hdata.anglePosi, anglePosi, size1, cudaMemcpyDeviceToHost);
    cudaMemcpy (hdata.angleNeighbor_globalIndex,
		angleNeighbor_globalIndex,
		size2,
		cudaMemcpyDeviceToHost);
  }
  if (expectedNumDihedral != 0){
    cudaMemcpy (hdata.numDihedral, numDihedral, size0, cudaMemcpyDeviceToHost);
    size_t size1 = size0 * maxNumDihedral;
    size_t size2 = size0 * maxNumDihedral * 3;
    cudaMemcpy (hdata.dihedralIndex, dihedralIndex, size1, cudaMemcpyDeviceToHost);
    cudaMemcpy (hdata.dihedralPosi, dihedralPosi, size1, cudaMemcpyDeviceToHost);
    cudaMemcpy (hdata.dihedralNeighbor_globalIndex,
		dihedralNeighbor_globalIndex,
		size2,
		cudaMemcpyDeviceToHost);
  }
}


void Parallel::DeviceMDData::
copyFromDevice (const DeviceMDData & ddata,
		const MDDataItemMask_t mask)
{
  // if (!mask) return;
  IndexType expectedNumBond(0), expectedNumAngle(0), expectedNumDihedral(0);
  bool copyBond = (mask & MDDataItemMask_Bond);
  bool copyAngle = (mask & MDDataItemMask_Angle);
  bool copyDihedral = (mask & MDDataItemMask_Dihedral);  
  if (copyBond){
    expectedNumBond = ddata.getMaxNumBond();
  }
  if (copyAngle){
    expectedNumAngle = ddata.getMaxNumAngle();
  }
  if (copyDihedral){
    expectedNumDihedral = ddata.getMaxNumDihedral();
  }

  if (memSize() != ddata.memSize() ||
      (copyBond && (maxNumBond != ddata.maxNumBond)) ||
      (copyAngle && (maxNumAngle != ddata.maxNumAngle)) ||
      (copyDihedral && (maxNumDihedral != ddata.maxNumDihedral)) ){
    easyMalloc (ddata.memSize(), expectedNumBond, expectedNumAngle, expectedNumDihedral);
  }

  numData() = ddata.numData();
  setGlobalBox (ddata.getGlobalBox());
		
  size_t sizef = ddata.numData() * sizeof(ScalorType);
  size_t sizecoord = ddata.numData() * sizeof(CoordType);
  size_t sizecoordNoi = ddata.numData() * sizeof(CoordNoiType);
  size_t sizeIdx = ddata.numData() * sizeof(IndexType);
  size_t sizet = ddata.numData() * sizeof(TypeType);

  if (mask & MDDataItemMask_Coordinate){
    cudaMemcpy (coord, ddata.coord, sizecoord, cudaMemcpyDeviceToDevice);
    checkCUDAError ("cpyDeviceMDDataToDevice coord");
  }

  if (mask & MDDataItemMask_CoordinateNoi){
    cudaMemcpy (coordNoi, ddata.coordNoi, sizecoordNoi, cudaMemcpyDeviceToDevice);
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

  size_t size0 = sizeof(IndexType) * ddata.memSize();

  if (copyBond && maxNumBond != 0){
    cudaMemcpy (numBond, ddata.numBond, size0, cudaMemcpyDeviceToDevice);
    size_t size1 = size0 * maxNumBond;
    cudaMemcpy (bondIndex, ddata.bondIndex, size1, cudaMemcpyDeviceToDevice);
    cudaMemcpy (bondNeighbor_globalIndex,
		ddata.bondNeighbor_globalIndex,
		size1,
		cudaMemcpyDeviceToDevice);
  }
  if (copyAngle && maxNumAngle != 0){
    cudaMemcpy (numAngle, ddata.numAngle, size0, cudaMemcpyDeviceToDevice);
    size_t size1 = size0 * maxNumAngle;
    size_t size2 = size0 * maxNumAngle * 2;
    cudaMemcpy (angleIndex, ddata.angleIndex, size1, cudaMemcpyDeviceToDevice);
    cudaMemcpy (anglePosi, ddata.anglePosi, size1, cudaMemcpyDeviceToDevice);
    cudaMemcpy (angleNeighbor_globalIndex,
		ddata.angleNeighbor_globalIndex,
		size2,
		cudaMemcpyDeviceToDevice);
  }
  if (copyDihedral && maxNumDihedral != 0){
    cudaMemcpy (numDihedral, ddata.numDihedral, size0, cudaMemcpyDeviceToDevice);
    size_t size1 = size0 * maxNumDihedral;
    size_t size2 = size0 * maxNumDihedral * 3;
    cudaMemcpy (dihedralIndex, ddata.dihedralIndex, size1, cudaMemcpyDeviceToDevice);
    cudaMemcpy (dihedralPosi, ddata.dihedralPosi, size1, cudaMemcpyDeviceToDevice);
    cudaMemcpy (dihedralNeighbor_globalIndex,
		ddata.dihedralNeighbor_globalIndex,
		size2,
		cudaMemcpyDeviceToDevice);
  }
}


Parallel::DeviceMDData::
DeviceMDData (const DeviceMDData & ddata)
    : _numData(0),
      _memSize(0), 
      maxNumBond(0),
      maxNumAngle(0),
      maxNumDihedral(0)
{
  copyFromDevice (ddata, MDDataItemMask_All);
}


// void Parallel::DeviceMDData::
// initZero ()
// {
//   Parallel::CudaGlobal::initZeroDeviceData
//       <<<memSize_ / DefaultNThreadPerBlock + 1, DefaultNThreadPerBlock>>>
//       (memSize_,
//        coord,
//        coordNoi,
//        velox,
//        veloy,
//        veloz,
//        forcx,
//        forcy,
//        forcz,
//        globalIndex,
//        type,
//        mass,
//        charge);
//   checkCUDAError ("DeviceMDData::initZero initZeroDeviceData");
// }


// __global__ void Parallel::CudaGlobal::
// initZeroDeviceData(const IndexType num,
// 		   CoordType  * coord,
// 		   CoordNoiType * coordNoi,
// 		   ScalorType * velox,
// 		   ScalorType * veloy,
// 		   ScalorType * veloz,
// 		   ScalorType * forcx,
// 		   ScalorType * forcy,
// 		   ScalorType * forcz,
// 		   IndexType  * globalIndex,
// 		   TypeType   * type,
// 		   ScalorType * mass,
// 		   ScalorType * charge)
// {
//   IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
//   IndexType tid = threadIdx.x;
//   IndexType ii = tid + bid * blockDim.x;

//   if (ii < num){
//     CoordType tmp;
//     tmp.x = 0;
//     tmp.y = 0;
//     tmp.z = 0;
//     tmp.w = -1;
//     coord[ii] = tmp;
//     coordNoi[ii].x = coordNoi[ii].y = coordNoi[ii].z = 0;
//     veloz[ii] = veloy[ii] = veloz[ii] = 0.f;
//     forcz[ii] = forcy[ii] = forcz[ii] = 0.f;
//     globalIndex[ii] = MaxIndexValue;
//     type[ii] = 0;
//     mass[ii] = 0;
//     charge[ii] = 0;
//   }
// }


void Parallel::GlobalHostMDData::
initWriteData_xtcFile (const char * filename, float prec)
{
  xdfile = NULL;
  xdfile = xdrfile_open (filename, "w");
  if (xdfile == NULL){
    MDExcptCannotOpenFile ("MDSystem::initWriteXtc", filename);
  }
  for (unsigned i = 0; i < 3; ++i){
    for (unsigned j = 0; j < 3; ++j){
      xdbox[i][j] = 0.f;
    }	      
  }
  xdx = (rvec *) malloc (sizeof(rvec) * numData());
  if (xdx == NULL){
    MDExcptFailedMallocOnHost ("MDSystem::initWriteXtc", "xdx", sizeof(rvec) * numData());
  }
  xdprec = prec;
}


void Parallel::GlobalHostMDData::
writeData_xtcFile (int step, float time)
{
  for (IndexType i = 0; i < numData(); ++i){
    xdx[i][0] = coord[i].x;
    xdx[i][1] = coord[i].y;
    xdx[i][2] = coord[i].z;
  }
  xdbox[0][0] = globalBox.size.x;
  xdbox[1][1] = globalBox.size.y;
  xdbox[2][2] = globalBox.size.z;
  write_xtc (xdfile, numData(), step, time, xdbox, xdx, xdprec);
}

void Parallel::GlobalHostMDData::
endWriteData_xtcFile ()
{
  free (xdx);
  xdrfile_close(xdfile);
}

void Parallel::DeviceMDData::
mallocFromDevice (const DeviceMDData & ddata)
{
  setGlobalBox (ddata.getGlobalBox());
  easyMalloc (ddata.memSize(), ddata.getMaxNumBond(), ddata.getMaxNumAngle(),
	      ddata.getMaxNumDihedral());
  _numData = 0;
}

void Parallel::DeviceMDData::
mallocFromHost (const HostMDData & hdata)
{
  setGlobalBox (hdata.getGlobalBox());
  easyMalloc (hdata.memSize(), hdata.getMaxNumBond(), hdata.getMaxNumAngle(),
	      hdata.getMaxNumDihedral());
  _numData = 0;
}


void Parallel::DeviceMDData::
mallocToHost (HostMDData & hdata) const 
{
  hdata.setGlobalBox (getGlobalBox());
  hdata.easyMalloc (memSize(), getMaxNumBond(), getMaxNumAngle(), getMaxNumDihedral());
  hdata.numData() = 0;
}









