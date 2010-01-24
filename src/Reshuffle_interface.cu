#include "Reshuffle_interface.h"
#include "Reshuffle.h"
#include "NonBondedInteraction.h"
#include "HSFC.h"

__global__ void init_bkNlistData (const IndexType numAtom, const IndexType length,
				  const IndexType stride,
				  IndexType * bknlistData,
				  IndexType * backMapTable);

HSFCMap1dto3d dmap1dto3d;
HSFCMap3dto1d dmap3dto1d;

void Reshuffle::init (const MDSystem & sys,
		      const NeighborList & nlist, 
		      const IndexType & NTread)
{
  myBlockDim.y = 1;
  myBlockDim.z = 1;
  myBlockDim.x = NTread;
  cellGridDim = nlist.cellGridDim;

  initDeviceMDData(&sys.hdata, &imageSystem);
  
  IndexType nob;
  if (sys.hdata.numAtom % myBlockDim.x == 0){
    nob = sys.hdata.numAtom / myBlockDim.x;
  } else {
    nob = sys.hdata.numAtom / myBlockDim.x + 1;
  }
  atomGridDim = toGridDim (nob);

  nob = cellGridDim.x * cellGridDim.y;;
  cudaMalloc ((void**)&posiBuff, sizeof(IndexType)*nob);
  cudaMalloc ((void**)&bknlistData, 
	      sizeof(IndexType) * nlist.dnlist.stride * nlist.dnlist.listLength);
  cudaMalloc ((void**)&bkNBForceIndex,
	      sizeof(ForceIndexType)* nlist.dnlist.stride * nlist.dnlist.listLength);
  cudaMalloc ((void**)&bkNneighbor,
	      sizeof(IndexType) * nlist.dnlist.stride);
  cudaMalloc ((void**)&idxTable, sizeof(IndexType) * sys.hdata.numAtom);
  cudaMalloc ((void**)&backMapTable, sizeof(IndexType) * sys.hdata.numAtom);
  cudaMalloc ((void**)&backMapTableBuff, sizeof(IndexType) * sys.hdata.numAtom);
  // bond list
  cudaMalloc ((void**)&bkBondListData, 
	      sizeof(IndexType) * sys.bdlist.dbdlist.stride * sys.bdlist.dbdlist.listLength);
  cudaMalloc ((void**)&bkBondListBondIndex, 
	      sizeof(ForceIndexType) * sys.bdlist.dbdlist.stride * sys.bdlist.dbdlist.listLength);
  cudaMalloc ((void**)&bkBondListNumB,
	      sizeof(IndexType) * sys.bdlist.dbdlist.stride);
  // angle list
  cudaMalloc ((void**)&bkAngleListNei,
	      sizeof(IndexType) * sys.anglelist.danglelist.stride *
	      sys.anglelist.danglelist.listLength * 2);
  cudaMalloc ((void**)&bkAngleListPosi,
	      sizeof(IndexType) * sys.anglelist.danglelist.stride *
	      sys.anglelist.danglelist.listLength);
  cudaMalloc ((void**)&bkAngleListAngleIndex,
	      sizeof(ForceIndexType) * sys.anglelist.danglelist.stride *
	      sys.anglelist.danglelist.listLength);
  cudaMalloc ((void**)&bkAngleListNangle,
	      sizeof(IndexType) * sys.anglelist.danglelist.stride);
#ifndef COORD_IN_ONE_VEC
  cudaMalloc ((void**)&bkNlistJudgeBuffx, sizeof(ScalorType)*sys.hdata.numAtom);
  cudaMalloc ((void**)&bkNlistJudgeBuffy, sizeof(ScalorType)*sys.hdata.numAtom);
  cudaMalloc ((void**)&bkNlistJudgeBuffz, sizeof(ScalorType)*sys.hdata.numAtom);
#else
  cudaMalloc ((void**)&bkNlistJudgeBuff, sizeof(CoordType) * sys.hdata.numAtom);
#endif
  checkCUDAError ("Reshuffle::init allocation");
  
  init_bkNlistData <<<atomGridDim, myBlockDim>>>
      (sys.hdata.numAtom, nlist.dnlist.listLength, nlist.dnlist.stride,
       bknlistData, backMapTable);

  // init HSFC
  // HSFCMap3dto1d hmap3dto1d;
  // HSFCMap1dto3d hmap1dto3d;
  // initHostHSFCMap1dto3d (&hmap1dto3d,
  // 			 nlist.dclist.NCell.x, 
  // 			 nlist.dclist.NCell.y,
  // 			 nlist.dclist.NCell.z);
  // initHostHSFCMap3dto1d (&hmap3dto1d,
  // 			 hmap1dto3d,
  // 			 nlist.dclist.NCell.x, 
  // 			 nlist.dclist.NCell.y,
  // 			 nlist.dclist.NCell.z);
  // initDeviceHSFCMap1dto3d (&dmap1dto3d,
  // 			   hmap1dto3d,
  // 			   nlist.dclist.NCell.x, 
  // 			   nlist.dclist.NCell.y,
  // 			   nlist.dclist.NCell.z);
  // initDeviceHSFCMap3dto1d (&dmap3dto1d,
  // 			   hmap3dto1d,
  // 			   nlist.dclist.NCell.x, 
  // 			   nlist.dclist.NCell.y,
  // 			   nlist.dclist.NCell.z);
  // checkCUDAError ("Reshuffle::init HSFC");
}

__global__ void init_bkNlistData (const IndexType numAtom, const IndexType length,
				  const IndexType stride,
				  IndexType * bknlistData,
				  IndexType * backMapTable)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;

  if (ii <   numAtom){
    for (IndexType jj = 0; jj < length; ++jj){
      bknlistData[jj * stride + ii] = MaxIndexValue;
    }
    backMapTable[ii] = ii;
  }
}


Reshuffle::~Reshuffle ()
{
  cudaFree(posiBuff);
  cudaFree(bknlistData);
  cudaFree(bkNBForceIndex);
  cudaFree(bkNneighbor);
  cudaFree(backMapTableBuff);
  cudaFree(backMapTable);
  cudaFree(bkBondListData);
  cudaFree(bkBondListBondIndex);
  cudaFree(bkBondListNumB);
#ifndef COORD_IN_ONE_VEC
  cudaFree(bkNlistJudgeBuffx);
  cudaFree(bkNlistJudgeBuffy);
  cudaFree(bkNlistJudgeBuffz);
#else
  cudaFree(bkNlistJudgeBuff);
#endif
  checkCUDAError ("Reshuffle::~Reshuffle");
}

void Reshuffle::shuffleSystem (MDSystem & sys,
			       NeighborList & nlist,
			       MDTimer * timer)
{
  if (nlist.mode == AllPairBuilt){
    return;
  }

  cellGridDim = nlist.cellGridDim;
  
  if (timer != NULL) timer->tic(mdTimeReshuffleSystem);
  IndexType nob = cellGridDim.x * cellGridDim.y;
  // possible streams
  // Reshuffle_backupSystem
  //     <<<atomGridDim, myBlockDim>>> (
  // 	  sys.ddata, imageSystem);

  Reshuffle_backupDeviceMDData_part1 
      <<<atomGridDim, myBlockDim>>> (
	  sys.ddata.numAtom,
#ifndef COORD_IN_ONE_VEC
	  sys.ddata.coordx, sys.ddata.coordy, sys.ddata.coordz,
#else
	  sys.ddata.coord,
#endif
	  sys.ddata.coordNoix, sys.ddata.coordNoiy, sys.ddata.coordNoiz,
#ifndef COORD_IN_ONE_VEC
	  imageSystem.coordx, imageSystem.coordy, imageSystem.coordz,
#else
	  imageSystem.coord,
#endif
	  imageSystem.coordNoix, imageSystem.coordNoiy, imageSystem.coordNoiz);
  Reshuffle_backupDeviceMDData_part2
      <<<atomGridDim, myBlockDim>>> (
	  sys.ddata.numAtom,
	  sys.ddata.velox, sys.ddata.veloy, sys.ddata.veloz,
	  sys.ddata.forcx, sys.ddata.forcy, sys.ddata.forcz,
	  imageSystem.velox, imageSystem.veloy, imageSystem.veloz,
	  imageSystem.forcx, imageSystem.forcy, imageSystem.forcz);
  Reshuffle_backupDeviceMDData_part3 
      <<<atomGridDim, myBlockDim>>> (
	  sys.ddata.numAtom,
	  sys.ddata.type, sys.ddata.mass, sys.ddata.massi, sys.ddata.charge,
	  imageSystem.type, imageSystem.mass, imageSystem.massi, imageSystem.charge);
  checkCUDAError ("Reshuffle::shuffleSystem, back up system");
  Reshuffle_backupNeighborLists
      <<<atomGridDim, myBlockDim,
      2 * myBlockDim.x * sizeof(IndexType)>>> (
	  sys.ddata.numAtom,
	  nlist.dnlist.data, nlist.dnlist.forceIndex,
	  nlist.dnlist.stride, nlist.dnlist.Nneighbor,
	  bknlistData, bkNBForceIndex, bkNneighbor);
#ifndef COORD_IN_ONE_VEC
  Reshuffle_backupScalorTypeBuff
      <<<atomGridDim, myBlockDim>>> (
	  sys.ddata.numAtom, nlist.backupCoordx, bkNlistJudgeBuffx);
  Reshuffle_backupScalorTypeBuff
      <<<atomGridDim, myBlockDim>>> (
	  sys.ddata.numAtom, nlist.backupCoordy, bkNlistJudgeBuffy);
  Reshuffle_backupScalorTypeBuff
      <<<atomGridDim, myBlockDim>>> (
	  sys.ddata.numAtom, nlist.backupCoordz, bkNlistJudgeBuffz);
#else
  Reshuffle_backupCoord
      <<<atomGridDim, myBlockDim>>> (
	  sys.ddata.numAtom, nlist.backupCoord, bkNlistJudgeBuff);
#endif
  checkCUDAError ("Reshuffle::shuffleSystem, back up neighbor list");
  // bond list
  Reshuffle_backupBondList
      <<<atomGridDim, myBlockDim>>> (
  	  sys.ddata.numAtom,
  	  sys.bdlist.dbdlist.data, 
	  sys.bdlist.dbdlist.bondIndex, 
	  sys.bdlist.dbdlist.Nbond,
  	  sys.bdlist.dbdlist.stride, 
	  sys.bdlist.dbdlist.listLength,
  	  bkBondListData,
	  bkBondListBondIndex,
	  bkBondListNumB);
  checkCUDAError ("Reshuffle::shuffleSystem, back up bond list");
  // angle list
  Reshuffle_backupAngleList
      <<<atomGridDim, myBlockDim>>>(
	  sys.ddata.numAtom,
	  sys.anglelist.danglelist.angleNei,
	  sys.anglelist.danglelist.myPosi,
	  sys.anglelist.danglelist.angleIndex,
	  sys.anglelist.danglelist.Nangle,
  	  sys.anglelist.danglelist.stride, 
	  sys.anglelist.danglelist.listLength,
	  bkAngleListNei,
	  bkAngleListPosi,
	  bkAngleListAngleIndex,
	  bkAngleListNangle);
  checkCUDAError ("Reshuffle::shuffleSystem, back up angle list");
  Reshuffle_backupBackMapTable 
      <<<atomGridDim, myBlockDim>>> (
	  sys.ddata.numAtom,
	  backMapTable,
	  backMapTableBuff);
  checkCUDAError ("Reshuffle::shuffleSystem, back up backMapTable");

  
  Reshuffle_calPosiList
      <<<1, 1>>> (
   	  nlist.dclist.numbers, nob, posiBuff);
  checkCUDAError ("Reshuffle::shuffleSystem, cal posi");
  Reshuffle_calIndexTable 
      <<<cellGridDim, myBlockDim>>> (
   	  nlist.dclist.data, posiBuff, idxTable);
  checkCUDAError ("Reshuffle::shuffleSystem, cal idxTable");
/*  Reshuffle_calPosiList_HSFC
    <<<1, 32>>> (
    nlist.dclist, posiBuff, dmap3dto1d);
    checkCUDAError ("Reshuffle::shuffleSystem, cal posi");
    Reshuffle_calIndexTable_HSFC
    <<<cellGridDim, myBlockDim>>> (
    nlist.dclist, posiBuff, idxTable, dmap3dto1d);
    checkCUDAError ("Reshuffle::shuffleSystem, cal idxTable");
*/

  Reshuffle_calBackMapTable
      <<<atomGridDim, myBlockDim>>> (
	  sys.ddata.numAtom,
	  backMapTableBuff, idxTable,
	  backMapTable);
  checkCUDAError ("Reshuffle::shuffleSystem, cal BackMapTable");
  Reshuffle_reshuffleNeighborList
      <<<atomGridDim, myBlockDim,
      2 * myBlockDim.x * sizeof(IndexType)>>> (
	  sys.ddata.numAtom,
	  bknlistData, bkNBForceIndex, nlist.dnlist.stride, bkNneighbor,
	  idxTable,
	  nlist.dnlist.data, nlist.dnlist.forceIndex,
	  nlist.dnlist.Nneighbor);
#ifndef COORD_IN_ONE_VEC
  Reshuffle_reshuffleScalorTypeBuff
      <<<atomGridDim, myBlockDim>>> (
	  sys.ddata.numAtom, bkNlistJudgeBuffx, idxTable, nlist.backupCoordx);
  Reshuffle_reshuffleScalorTypeBuff
      <<<atomGridDim, myBlockDim>>> (
	  sys.ddata.numAtom, bkNlistJudgeBuffy, idxTable, nlist.backupCoordy);
  Reshuffle_reshuffleScalorTypeBuff
      <<<atomGridDim, myBlockDim>>> (
	  sys.ddata.numAtom, bkNlistJudgeBuffz, idxTable, nlist.backupCoordz);
#else
  Reshuffle_reshuffleCoord
      <<<atomGridDim, myBlockDim>>> (
	  sys.ddata.numAtom, bkNlistJudgeBuff, idxTable, nlist.backupCoord);
#endif
  checkCUDAError ("Reshuffle::shuffleSystem, reshuffle neighbor list");
  Reshuffle_reshuffleBondList
      <<<atomGridDim,myBlockDim>>> (
  	  sys.ddata.numAtom,
  	  bkBondListData, 
	  bkBondListBondIndex, 
	  bkBondListNumB,
  	  sys.bdlist.dbdlist.stride,
	  sys.bdlist.dbdlist.listLength, 
  	  idxTable,
  	  sys.bdlist.dbdlist.data, 
	  sys.bdlist.dbdlist.bondIndex, 
	  sys.bdlist.dbdlist.Nbond);
  checkCUDAError ("Reshuffle::shuffleSystem, reshuffle bond list");
  Reshuffle_reshuffleAngleList
      <<<atomGridDim,myBlockDim>>> (
	  sys.ddata.numAtom,
	  bkAngleListNei,
	  bkAngleListPosi,
	  bkAngleListAngleIndex,
	  bkAngleListNangle,
  	  sys.anglelist.danglelist.stride,
	  sys.anglelist.danglelist.listLength, 
  	  idxTable,
	  sys.anglelist.danglelist.angleNei,
	  sys.anglelist.danglelist.myPosi,
	  sys.anglelist.danglelist.angleIndex,
	  sys.anglelist.danglelist.Nangle);
  checkCUDAError ("Reshuffle::shuffleSystem, reshuffle angle list");
  Reshuffle_reshuffleCellList
      <<<cellGridDim, myBlockDim>>> (
	  nlist.dclist.data, idxTable, posiBuff);
  checkCUDAError ("Reshuffle::shuffleSystem, reshuffle cell list");
  // Reshuffle_reshuffleDeviceMDData
  //     <<<atomGridDim, myBlockDim>>> (
  // 	  imageSystem, idxTable, sys.ddata);
  Reshuffle_reshuffleDeviceMDData_part1 
      <<<atomGridDim, myBlockDim>>> (
	  sys.ddata.numAtom,
#ifndef COORD_IN_ONE_VEC
	  imageSystem.coordx, 
	  imageSystem.coordy, 
	  imageSystem.coordz,
#else
	  imageSystem.coord,
#endif
	  imageSystem.coordNoix, 
	  imageSystem.coordNoiy, 
	  imageSystem.coordNoiz,
	  idxTable,
#ifndef COORD_IN_ONE_VEC
	  sys.ddata.coordx,
	  sys.ddata.coordy,
	  sys.ddata.coordz,
#else
	  sys.ddata.coord,
#endif
	  sys.ddata.coordNoix,
	  sys.ddata.coordNoiy, 
	  sys.ddata.coordNoiz);
  Reshuffle_reshuffleDeviceMDData_part2
      <<<atomGridDim, myBlockDim>>> (
	  sys.ddata.numAtom,
	  imageSystem.velox, imageSystem.veloy, imageSystem.veloz,
	  imageSystem.forcx, imageSystem.forcy, imageSystem.forcz,
	  idxTable,
	  sys.ddata.velox, sys.ddata.veloy, sys.ddata.veloz,
	  sys.ddata.forcx, sys.ddata.forcy, sys.ddata.forcz);
  Reshuffle_reshuffleDeviceMDData_part3
      <<<atomGridDim, myBlockDim>>> (
	  sys.ddata.numAtom,
	  imageSystem.type, 
	  imageSystem.mass, 
	  imageSystem.massi, 
	  imageSystem.charge,
	  idxTable,
	  sys.ddata.type, 
	  sys.ddata.mass,
	  sys.ddata.massi, 
	  sys.ddata.charge);
  checkCUDAError ("Reshuffle::shuffleSystem, reshuffle");
  if (timer != NULL) timer->toc(mdTimeReshuffleSystem);
}


void Reshuffle::recoverMDData (const DeviceMDData & currentData,
			       DeviceMDData & recoveredData,
			       MDTimer * timer)
{
  if (timer != NULL) timer->tic(mdTimeDataIO);
  Reshuffle_reshuffleDeviceMDData_part1 
      <<<atomGridDim, myBlockDim>>> (
	  currentData.numAtom,
#ifndef COORD_IN_ONE_VEC
	  currentData.coordx, currentData.coordy, currentData.coordz,
#else
	  currentData.coord,
#endif
	  currentData.coordNoix, currentData.coordNoiy, currentData.coordNoiz,
	  backMapTable,
#ifndef COORD_IN_ONE_VEC
	  recoveredData.coordx, recoveredData.coordy, recoveredData.coordz,
#else
	  recoveredData.coord,
#endif
	  recoveredData.coordNoix, recoveredData.coordNoiy, recoveredData.coordNoiz);
  Reshuffle_reshuffleDeviceMDData_part2
      <<<atomGridDim, myBlockDim>>> (
	  currentData.numAtom,
	  currentData.velox, currentData.veloy, currentData.veloz,
	  currentData.forcx, currentData.forcy, currentData.forcz,
	  backMapTable,
	  recoveredData.velox, recoveredData.veloy, recoveredData.veloz,
	  recoveredData.forcx, recoveredData.forcy, recoveredData.forcz);
  Reshuffle_reshuffleDeviceMDData_part3
      <<<atomGridDim, myBlockDim>>> (
	  currentData.numAtom,
	  currentData.type, currentData.mass, currentData.massi, currentData.charge,
	  backMapTable,
	  recoveredData.type, recoveredData.mass, 
	  recoveredData.massi, recoveredData.charge);
  checkCUDAError ("Reshuffle::recoveredData, recover");
  if (timer != NULL) timer->toc(mdTimeDataIO);
}

void Reshuffle::recoverMDDataToHost (MDSystem & sys   ,
				     MDTimer * timer)
{
  recoverMDData (sys.ddata, sys.recoveredDdata, timer);
  if (timer != NULL) timer->tic(mdTimeDataTransfer);
  cpyDeviceMDDataToHost (&(sys.recoveredDdata), &(sys.hdata));
  if (timer != NULL) timer->toc(mdTimeDataTransfer);
}




// __global__ void Reshuffle_backupSystem (const DeviceMDData data1,
// 					DeviceMDData data2)
// {
//   IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
//   IndexType tid = threadIdx.x;
//   IndexType ii = tid + bid * blockDim.x;

//   if (ii < data1.numAtom )
//     cpyDeviceMDDataElement (&data1, ii, &data2, ii);
// }


#ifndef COORD_IN_ONE_VEC
__global__ void Reshuffle_backupDeviceMDData_part1 (IndexType numAtom,
						    const ScalorType * coordx1,
						    const ScalorType * coordy1,
						    const ScalorType * coordz1,
						    const IntScalorType * coordNoix1,
						    const IntScalorType * coordNoiy1,
						    const IntScalorType * coordNoiz1,
						    ScalorType * coordx2,
						    ScalorType * coordy2,
						    ScalorType * coordz2,
						    IntScalorType * coordNoix2,
						    IntScalorType * coordNoiy2,
						    IntScalorType * coordNoiz2)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;

  if (ii < numAtom){
    coordx2[ii] = coordx1[ii];
    coordy2[ii] = coordy1[ii];
    coordz2[ii] = coordz1[ii];
    coordNoix2[ii] = coordNoix1[ii];
    coordNoiy2[ii] = coordNoiy1[ii];
    coordNoiz2[ii] = coordNoiz1[ii];
  }
}
#else
__global__ void Reshuffle_backupDeviceMDData_part1 (IndexType numAtom,
						    const CoordType * coord1,
						    const IntScalorType * coordNoix1,
						    const IntScalorType * coordNoiy1,
						    const IntScalorType * coordNoiz1,
						    CoordType * coord2,
						    IntScalorType * coordNoix2,
						    IntScalorType * coordNoiy2,
						    IntScalorType * coordNoiz2)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;

  if (ii < numAtom){
    coord2[ii] = coord1[ii];
    coordNoix2[ii] = coordNoix1[ii];
    coordNoiy2[ii] = coordNoiy1[ii];
    coordNoiz2[ii] = coordNoiz1[ii];
  }
}
#endif

__global__ void Reshuffle_backupDeviceMDData_part2 (IndexType numAtom,
						    const ScalorType * velox1,
						    const ScalorType * veloy1,
						    const ScalorType * veloz1,
						    const ScalorType * forcx1,
						    const ScalorType * forcy1,
						    const ScalorType * forcz1,
						    ScalorType * velox2,
						    ScalorType * veloy2,
						    ScalorType * veloz2,
						    ScalorType * forcx2,
						    ScalorType * forcy2,
						    ScalorType * forcz2)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;

  if (ii < numAtom){
    velox2[ii] = velox1[ii];
    veloy2[ii] = veloy1[ii];
    veloz2[ii] = veloz1[ii];
    forcx2[ii] = forcx1[ii];
    forcy2[ii] = forcy1[ii];
    forcz2[ii] = forcz1[ii];
  }
}
						       

__global__ void Reshuffle_backupDeviceMDData_part3 (IndexType numAtom,
						    const TypeType * type1,
						    const ScalorType * mass1,
						    const ScalorType * massi1,
						    const ScalorType * charge1,
						    TypeType * type2,
						    ScalorType * mass2,
						    ScalorType * massi2,
						    ScalorType * charge2)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;

  if (ii < numAtom){
    type2[ii] = type1[ii];
    mass2[ii] = mass1[ii];
    massi2[ii] = massi1[ii];
    charge2[ii] = charge1[ii];
  }
}

#ifdef COORD_IN_ONE_VEC
__global__ void Reshuffle_backupCoord (IndexType numAtom,
				       const CoordType * buff,
				       CoordType * bkbuff)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  if (ii < numAtom){
    bkbuff[ii] = buff[ii];
  }
}
#endif

__global__ void Reshuffle_backupScalorTypeBuff (IndexType numAtom,
						const ScalorType * buff,
						ScalorType * bkbuff)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  if (ii < numAtom){
    bkbuff[ii] = buff[ii];
  }
}
  
__global__ void Reshuffle_backupIntScalorTypeBuff (IndexType numAtom,
						   const IntScalorType * buff,
						   IntScalorType * bkbuff)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  if (ii < numAtom){
    bkbuff[ii] = buff[ii];
  }
}
  
__global__ void Reshuffle_backupTypeTypeBuff (IndexType numAtom,
					      const TypeType * buff,
					      TypeType * bkbuff)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  if (ii < numAtom){
    bkbuff[ii] = buff[ii];
  }
}

__global__ void Reshuffle_backupIndexTypeBuff (IndexType numAtom,
					       const IndexType * buff,
					       IndexType * bkbuff)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  if (ii < numAtom){
    bkbuff[ii] = buff[ii];
  }
}  


__global__ void Reshuffle_backupNeighborLists (const IndexType numAtom,
					      const IndexType * nlistData1,
					      const ForceIndexType * nbForceIndex1,
					      const IndexType stride,
					      const IndexType * Nneighbor1,
					      IndexType * nlistData2,
					      ForceIndexType * nbForceIndex2,
					      IndexType * Nneighbor2)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;

  // __shared__ volatile IndexType myNumbers [MaxThreadsPerBlock * 2];
  extern __shared__ volatile IndexType myNumbers[];
  
  IndexType N;
  if ((bid + 1) * blockDim.x < numAtom) N = blockDim.x;
  else if (bid * blockDim.x >= numAtom) N = 0;
  else N = numAtom - bid * blockDim.x;

  myNumbers[tid] = 0;
  myNumbers[tid + blockDim.x] = 0;
  if (ii < numAtom){
    Nneighbor2[ii] = myNumbers[tid] = Nneighbor1[ii];
  }
  __syncthreads();
  IndexType maxNum = maxVectorBlockBuffer (myNumbers, N);
  __syncthreads();

  for (IndexType jj = 0; jj < maxNum; ++jj){
    if (jj < myNumbers[tid] && ii < numAtom ){
      nlistData2   [jj * stride + ii] = nlistData1   [jj * stride + ii];
      nbForceIndex2[jj * stride + ii] = nbForceIndex1[jj * stride + ii];
    }
  }
}


__global__ void
Reshuffle_backupBondList (const IndexType numAtom,
			  const IndexType * bdlistData,
			  const ForceIndexType * bdlistBondIndex,
			  const IndexType * bdlistNumB,
			  const IndexType stride,
			  const IndexType listLength,
			  IndexType * bdlistData2,
			  ForceIndexType * bdlistBondIndex2,
			  IndexType * bdlistNumB2)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;

  if (ii < numAtom){
    bdlistNumB2[ii] = bdlistNumB[ii];
    for (IndexType jj = 0; jj < listLength; ++jj){
      bdlistData2     [jj * stride + ii] = bdlistData     [jj * stride + ii];
      bdlistBondIndex2[jj * stride + ii] = bdlistBondIndex[jj * stride + ii];
    }
  }
}

__global__ void
Reshuffle_backupAngleList (const IndexType numAtom,
			   const IndexType * angleListNei,
			   const IndexType * angleListPosi,
			   const ForceIndexType * angleListAngleIndex,
			   const IndexType * anglelistNangle,
			   const IndexType stride,
			   const IndexType listLength,
			   IndexType * bkAngleListNei,
			   IndexType * bkAngleListPosi,
			   ForceIndexType * bkAngleListAngleIndex,
			   IndexType * bkAngleListNangle)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;

  if (ii < numAtom){
    bkAngleListNangle[ii] = anglelistNangle[ii];
    for (IndexType jj = 0; jj < listLength; ++ jj){
      bkAngleListNei[((jj<<1)+0)*stride + ii] =angleListNei[((jj<<1)+0)*stride + ii];
      bkAngleListNei[((jj<<1)+1)*stride + ii] =angleListNei[((jj<<1)+1)*stride + ii];
      bkAngleListAngleIndex[jj * stride + ii] =angleListAngleIndex[jj * stride + ii];
      bkAngleListPosi[jj * stride + ii] = angleListPosi[jj * stride + ii];
    }
  }
}



__global__ void Reshuffle_backupBackMapTable (const IndexType numAtom,
					      const IndexType * backMapTable,
					      IndexType * backMapTableBuff)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;

  if (ii < numAtom){
    backMapTableBuff[ii] = backMapTable[ii];
  }
}

  


__global__ void Reshuffle_calPosiList_HSFC  (const DeviceCellList clist,
					     IndexType * posiList,
					     HSFCMap3dto1d dmap3dto1d)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  if (bid > 0) return ;
  IndexType tid = threadIdx.x;
  IndexType ncell = clist.NCell.x * clist.NCell.y * clist.NCell.z;
  
  IndexType p = 0;
  while (p < ncell){
    IndexType i, j, k;
    IndexType posi;
    IndexType cellid = p + tid;
    if (cellid >= ncell) break;
    D1toD3 (clist, cellid, i, j, k);
    map3dto1d_d (dmap3dto1d,
	       clist.NCell.x, clist.NCell.y, clist.NCell.z,
	       i, j, k,
	       & posi);
    // printf ("%03d %03d %03d  %03d\n", i, j, k, posi);
    posiList[posi] = clist.numbers[cellid];
    p += blockDim.x;
  }

  __syncthreads();
  if (tid + bid == 0){
    IndexType previousNatom;
    IndexType thisNatom = posiList[0];
    posiList[0] = 0;
    for (IndexType ii = 1; ii < ncell; ++ii){
      previousNatom = thisNatom;
      thisNatom = posiList[ii];
      posiList[ii] = posiList[ii-1] + previousNatom;
    }
  }
}

__global__ void Reshuffle_calIndexTable_HSFC  (const DeviceCellList clist,
					       const IndexType * posiList,
					       IndexType * idxTable,
					       HSFCMap3dto1d dmap3dto1d)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  
  IndexType i, j, k;
  IndexType posi;
  D1toD3 (clist, bid, i, j, k);
  map3dto1d_d (dmap3dto1d,
	     clist.NCell.x, clist.NCell.y, clist.NCell.z,
	     i, j, k,
	     & posi);
  IndexType fromid = clist.data[bid * blockDim.x + tid];
  IndexType toid = posiList[posi] + tid;

  if (fromid != MaxIndexValue){
    idxTable[fromid] = toid;
  }
}

__global__ void Reshuffle_calPosiList  (const IndexType * cellNumbers,
					const IndexType nob,
					IndexType * posiList)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;

  if (tid + bid == 0){
    posiList[0] = 0;
    for (IndexType ii = 1; ii < nob; ++ii){
      posiList[ii] = posiList[ii-1] + cellNumbers[ii-1];
    }
  }
}

__global__ void Reshuffle_calIndexTable (const IndexType * clistData,
					 const IndexType * posiList,
					 IndexType * idxTable)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;

  IndexType fromid = clistData[bid * blockDim.x + tid];
  IndexType toid = posiList[bid] + tid;
  if (fromid != MaxIndexValue){
    idxTable[fromid] = toid;
  }
}

__global__ void Reshuffle_calBackMapTable (const IndexType numAtom,
					   const IndexType * backMapTableBuff,
					   const IndexType * idxTable,
					   IndexType *backMapTable)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  if (ii < numAtom){
    backMapTable[idxTable[ii]] = backMapTableBuff[ii];
  }
}
  

__global__ void Reshuffle_reshuffleCellList (IndexType * clistData,
					     const IndexType * posiList)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;

  IndexType fromid = clistData[bid * blockDim.x + tid];
  IndexType toid = posiList[bid] + tid;
  if (fromid != MaxIndexValue){
    clistData[bid * blockDim.x + tid] = toid;
  }
}

__global__ void Reshuffle_reshuffleCellList (IndexType * clistData,
					     const IndexType * idxTable,
					     const IndexType * posiList)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;

  IndexType fromid = clistData[bid * blockDim.x + tid];
  if (fromid != MaxIndexValue){
    IndexType toid = idxTable[fromid];
    clistData[bid * blockDim.x + tid] = toid;
  }
}


__global__ void Reshuffle_reshuffleNeighborList (const IndexType numAtom,
						 const IndexType * nlistData1,
						 const ForceIndexType* nbForceIndex1,
						 const IndexType stride,
						 const IndexType * Nneighbor1,
						 const IndexType * idxTable,
						 IndexType * nlistData2,
						 ForceIndexType * nbForceIndex2,
						 IndexType * Nneighbor2)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;
  // fromid is  ii

  extern __shared__ volatile IndexType myNumbers [];
  // __shared__ volatile IndexType myNumbers [MaxThreadsPerBlock * 2];

  IndexType toid;
  if (ii < numAtom)
    toid = idxTable[ii];
  IndexType myNum = 0;
  
  myNumbers[tid] = 0;
  myNumbers[tid + blockDim.x] = 0;
  if (ii < numAtom){
    myNum = myNumbers[tid] = Nneighbor1[ii];
  }
  __syncthreads();
  IndexType maxNum = maxVectorBlockBuffer (myNumbers, blockDim.x);
  __syncthreads();

  if (ii < numAtom){
    Nneighbor2[toid] = Nneighbor1[ii];
  }
  for (unsigned jj = 0; jj < maxNum; ++jj){
    if (jj < myNum && ii < numAtom){
      nlistData2[jj * stride + toid] = idxTable[nlistData1[jj * stride + ii]];
      nbForceIndex2[jj * stride + toid] = nbForceIndex1[jj * stride + ii];
    }
  }
}

__global__ void
Reshuffle_reshuffleBondList (const IndexType numAtom,
			     const IndexType * bdlistData2,
			     const ForceIndexType * bdlistBondIndex2,
			     const IndexType * bdlistNumB2,
			     const IndexType stride,
			     const IndexType listLength,
			     const IndexType * idxTable,
			     IndexType * bdlistData,
			     ForceIndexType * bdlistBondIndex,
			     IndexType * bdlistNumB)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;

  if (ii < numAtom){
    IndexType toid = idxTable[ii];
    bdlistNumB[toid] = bdlistNumB2[ii];
    for (IndexType jj = 0; jj < bdlistNumB2[ii]; ++jj){
      bdlistData    [jj * stride + toid] = 
	  idxTable[bdlistData2[jj * stride + ii]];
      bdlistBondIndex[jj * stride + toid] = 
	  bdlistBondIndex2[jj * stride + ii];
    }
  }
}

__global__ void
Reshuffle_reshuffleAngleList (const IndexType numAtom,
			      const IndexType * bkAngleListNei,
			      const IndexType * bkAngleListPosi,
			      const ForceIndexType * bkAngleListAngleIndex,
			      const IndexType * bkAngleListNangle,
			      const IndexType stride,
			      const IndexType listLength,
			      const IndexType * idxTable,
			      IndexType * angleListNei,
			      IndexType * angleListPosi,
			      ForceIndexType * angleListAngleIndex,
			      IndexType * angleListNangle)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;

  if (ii < numAtom){
    IndexType toid = idxTable[ii];
    angleListNangle[toid] = bkAngleListNangle[ii];
    for (IndexType jj = 0; jj < bkAngleListNangle[ii]; ++ jj){
      angleListNei[((jj<<1)+0)*stride + toid] =
	  idxTable[bkAngleListNei[((jj<<1)+0)*stride + ii]];
      angleListNei[((jj<<1)+1)*stride + toid] =
	  idxTable[bkAngleListNei[((jj<<1)+1)*stride + ii]];
      angleListPosi[jj * stride + toid] = bkAngleListPosi[jj * stride + ii];
      angleListAngleIndex[jj * stride + toid] =
	  bkAngleListAngleIndex[jj * stride + ii];
    }
  }
}



// __global__ void Reshuffle_reshuffleDeviceMDData (const DeviceMDData ddata,
// 						 const IndexType * idxTable,
// 						 DeviceMDData ddata2)
// {
//   IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
//   IndexType tid = threadIdx.x;
//   IndexType ii = tid + bid * blockDim.x;

//   if (ii < ddata2.numAtom){
//     cpyDeviceMDDataElement (&ddata1, ii, &ddata2, idxTable[ii]);
//   }
// }


__global__ void Reshuffle_reshuffleScalorTypeBuff (IndexType numAtom,
						   const ScalorType * bkbuff,
						   const IndexType * idxTable,
						   ScalorType * buff)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  if (ii < numAtom){
    buff[idxTable[ii]] = bkbuff[ii];
  }
}

#ifdef COORD_IN_ONE_VEC
__global__ void Reshuffle_reshuffleCoord (IndexType numAtom,
					  const CoordType * bkbuff,
					  const IndexType * idxTable,
					  CoordType * buff)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  if (ii < numAtom){
    buff[idxTable[ii]] = bkbuff[ii];
  }
}
#endif

__global__ void Reshuffle_reshuffleTypeTypeBuff (IndexType numAtom,
						 const TypeType * bkbuff,
						 const IndexType * idxTable,
						 TypeType * buff)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  if (ii < numAtom){
    buff[idxTable[ii]] = bkbuff[ii];
  }
}

__global__ void Reshuffle_reshuffleIntScalorTypeBuff (IndexType numAtom,
						      const IntScalorType * bkbuff,
						      const IndexType * idxTable,
						      IntScalorType * buff)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  if (ii < numAtom){
    buff[idxTable[ii]] = bkbuff[ii];
  }
}

__global__ void Reshuffle_reshuffleIndexTypeBuff (IndexType numAtom,
						  const IndexType * bkbuff,
						  const IndexType * idxTable,
						  IndexType * buff)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  if (ii < numAtom){
    buff[idxTable[ii]] = bkbuff[ii];
  }
}
  
#ifndef COORD_IN_ONE_VEC
__global__ void Reshuffle_reshuffleDeviceMDData_part1 (IndexType numAtom,
						       const ScalorType * coordx2,
						       const ScalorType * coordy2,
						       const ScalorType * coordz2,
						       const IntScalorType * coordNoix2,
						       const IntScalorType * coordNoiy2,
						       const IntScalorType * coordNoiz2,
						       const IndexType * idxTable,
						       ScalorType * coordx1,
						       ScalorType * coordy1,
						       ScalorType * coordz1,
						       IntScalorType * coordNoix1,
						       IntScalorType * coordNoiy1,
						       IntScalorType * coordNoiz1)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;
  
  if (ii < numAtom){
    IndexType toid = idxTable[ii];
    coordx1[toid] = coordx2[ii];
    coordy1[toid] = coordy2[ii];
    coordz1[toid] = coordz2[ii];
    coordNoix1[toid] = coordNoix2[ii];
    coordNoiy1[toid] = coordNoiy2[ii];
    coordNoiz1[toid] = coordNoiz2[ii];
  }
}
#else
__global__ void Reshuffle_reshuffleDeviceMDData_part1 (IndexType numAtom,
						       const CoordType * coord2,
						       const IntScalorType * coordNoix2,
						       const IntScalorType * coordNoiy2,
						       const IntScalorType * coordNoiz2,
						       const IndexType * idxTable,
						       CoordType * coord1,
						       IntScalorType * coordNoix1,
						       IntScalorType * coordNoiy1,
						       IntScalorType * coordNoiz1)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;
  
  if (ii < numAtom){
    IndexType toid = idxTable[ii];
    coord1[toid] = coord2[ii];
    coordNoix1[toid] = coordNoix2[ii];
    coordNoiy1[toid] = coordNoiy2[ii];
    coordNoiz1[toid] = coordNoiz2[ii];
  }
}
#endif


__global__ void Reshuffle_reshuffleDeviceMDData_part2 (IndexType numAtom,
						       const ScalorType * velox2,
						       const ScalorType * veloy2,
						       const ScalorType * veloz2,
						       const ScalorType * forcx2,
						       const ScalorType * forcy2,
						       const ScalorType * forcz2,
						       const IndexType * idxTable,
						       ScalorType * velox1,
						       ScalorType * veloy1,
						       ScalorType * veloz1,
						       ScalorType * forcx1,
						       ScalorType * forcy1,
						       ScalorType * forcz1)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;

  if (ii < numAtom){
    IndexType toid = idxTable[ii];
    velox1[toid] = velox2[ii];
    veloy1[toid] = veloy2[ii];
    veloz1[toid] = veloz2[ii];
    forcx1[toid] = forcx2[ii];
    forcy1[toid] = forcy2[ii];
    forcz1[toid] = forcz2[ii];
  }
}
						       

__global__ void Reshuffle_reshuffleDeviceMDData_part3 (IndexType numAtom,
						       const TypeType * type2,
						       const ScalorType * mass2,
						       const ScalorType * massi2,
						       const ScalorType * charge2,
						       const IndexType * idxTable,
						       TypeType * type1,
						       ScalorType * mass1,
						       ScalorType * massi1,
						       ScalorType * charge1)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;

  if (ii < numAtom){
    IndexType toid = idxTable[ii];
    type1[toid] = type2[ii];
    mass1[toid] = mass2[ii];
    massi1[toid] = massi2[ii];
    charge1[toid] = charge2[ii];
  }
}
						       
