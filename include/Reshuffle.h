#ifndef __Reshuffle_h_wanghan__
#define __Reshuffle_h_wanghan__

#include "common.h"
#include "HSFC.h"

// needs ceil(numAtom/blockDim.x) blocks
// __global__ void Reshuffle_backupSystem (const DeviceMDData data1,
// 					DeviceMDData data2);
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
						    IntScalorType * coordNoiz2);
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
						    ScalorType * forcz2);
__global__ void Reshuffle_backupDeviceMDData_part3 (IndexType numAtom,
						    const TypeType * type1,
						    const ScalorType * mass1,
						    const ScalorType * massi1,
						    const ScalorType * charge1,
						    TypeType * type2,
						    ScalorType * mass2,
						    ScalorType * massi2,
						    ScalorType * charge2);
// needs ceil(numAtom/blockDim.x) blocks
__global__ void Reshuffle_backupScalorTypeBuff (IndexType numAtom,
						const ScalorType * buff,
						ScalorType * bkbuff);
__global__ void Reshuffle_backupIntScalorTypeBuff (IndexType numAtom,
						   const IntScalorType * buff,
						   IntScalorType * bkbuff);
__global__ void Reshuffle_backupTypeTypeBuff (IndexType numAtom,
					      const TypeType * buff,
					      TypeType * bkbuff);
__global__ void Reshuffle_backupIndexTypeBuff (IndexType numAtom,
					       const IndexType * buff,
					       IndexType * bkbuff);
__global__ void Reshuffle_reshuffleScalorTypeBuff (IndexType numAtom,
						   const ScalorType * bkbuff,
						   const IndexType * idxTable,
						   ScalorType * buff);
__global__ void Reshuffle_reshuffleIntScalorTypeBuff (IndexType numAtom,
						      const IntScalorType * bkbuff,
						      const IndexType * idxTable,
						      IntScalorType * buff);
__global__ void Reshuffle_reshuffleTypeTypeBuff (IndexType numAtom,
						 const TypeType * bkbuff,
						 const IndexType * idxTable,
						 TypeType * buff);
__global__ void Reshuffle_reshuffleIndexTypeBuff (IndexType numAtom,
						  const IndexType * bkbuff,
						  const IndexType * idxTable,
						  IndexType * buff);

// needs ceil(numAtom/blockDim.x) blocks
__global__ void Reshuffle_backupNeighborLists (const IndexType numAtom,
					      const IndexType * nlistData1,
					      const ForceIndexType * nbForceIndex1,
					      const IndexType stride,
					      const IndexType * Nneighbor1,
					      IndexType * nlistData2,
					      ForceIndexType * nbForceIndex2,
					      IndexType * Nneighbor2);
// needs ceil(numAtom/blockDim.x) blocks
__global__ void Reshuffle_backupBondList (const IndexType numAtom,
					  const IndexType * bdlistData,
					  const ForceIndexType * bdlistBondIndex,
					  const IndexType * bdlistNumB,
					  const IndexType stride,
					  const IndexType listLength,
					  IndexType * bdlistData2,
					  ForceIndexType * bdlistBondIndex2,
					  IndexType * bdlistNumB2);
// needs ceil(numAtom/blockDim.x) blocks
__global__ void Reshuffle_backupBackMapTable (const IndexType numAtom,
					      const IndexType * backMapTable,
					      IndexType * backMapTableBuff);

// needs 1 thread
__global__ void Reshuffle_calPosiList  (const IndexType * cellNumbers,
					const IndexType nob,
					IndexType * posiList);
// needs NCell
__global__ void Reshuffle_calIndexTable (const IndexType * clistData,
					 const IndexType * posiList,
					 IndexType * idxTable);
// something like 16/32 threads
__global__ void Reshuffle_calPosiList_HSFC  (const DeviceCellList clist,
					     IndexType * posiList,
					     HSFCMap3dto1d dmap3dto1d);
// needs NCell
__global__ void Reshuffle_calIndexTable_HSFC  (const DeviceCellList clist,
					       const IndexType * posiList,
					       IndexType * idxTable,
					       HSFCMap3dto1d dmap3dto1d);
// needs NCell
__global__ void Reshuffle_reshuffleCellList (IndexType * clistData,
					     const IndexType * posiList);
__global__ void Reshuffle_reshuffleCellList (IndexType * clistData,
					     const IndexType * idxTable,
					     const IndexType * posiList);
// needs ceil(numAtom/blockDim.x) blocks
__global__ void Reshuffle_calBackMapTable (const IndexType numAtom,
					   const IndexType * backMapTableBuff,
					   const IndexType * idxTable,
					   IndexType *backMapTable);
// needs ceil(numAtom/blockDim.x) blocks
__global__ void Reshuffle_reshuffleNeighborList (const IndexType numAtom,
						 const IndexType * nlistData2,
						 const ForceIndexType* nbForceIndex2,
						 const IndexType stride,
						 const IndexType * Nneighbor2,
						 const IndexType * idxTable,
						 IndexType * nlistData1,
						 IndexType * Nneighbor1,
						 ForceIndexType * nbForceIndex1);
// needs ceil(numAtom/blockDim.x) blocks
__global__ void Reshuffle_reshuffleBondList (const IndexType numAtom,
					     const IndexType * bdlistData2,
					     const ForceIndexType * bdlistBondIndex2,
					     const IndexType * bdlistNumB2,
					     const IndexType stride,
					     const IndexType listLength,
					     const IndexType * idxTable,
					     IndexType * bdlistData,
					     ForceIndexType  * bdlistBondIndex,
					     IndexType * bdlistNumB);
// needs ceil(numAtom/blockDim.x) blocks
// __global__ void Reshuffle_reshuffleDeviceMDData (const DeviceMDData ddata1,
// 						 const IndexType * idxTable,
// 						 DeviceMDData ddata2);
__global__ void Reshuffle_reshuffleDeviceMDData_part1 (IndexType numAtom,
						       const ScalorType * coordx1,
						       const ScalorType * coordy1,
						       const ScalorType * coordz1,
						       const IntScalorType * coordNoix1,
						       const IntScalorType * coordNoiy1,
						       const IntScalorType * coordNoiz1,
						       const IndexType * idxTable,
						       ScalorType * coordx2,
						       ScalorType * coordy2,
						       ScalorType * coordz2,
						       IntScalorType * coordNoix2,
						       IntScalorType * coordNoiy2,
						       IntScalorType * coordNoiz2);
__global__ void Reshuffle_reshuffleDeviceMDData_part2 (IndexType numAtom,
						       const ScalorType * velox1,
						       const ScalorType * veloy1,
						       const ScalorType * veloz1,
						       const ScalorType * forcx1,
						       const ScalorType * forcy1,
						       const ScalorType * forcz1,
						       const IndexType * idxTable,
						       ScalorType * velox2,
						       ScalorType * veloy2,
						       ScalorType * veloz2,
						       ScalorType * forcx2,
						       ScalorType * forcy2,
						       ScalorType * forcz2);
__global__ void Reshuffle_reshuffleDeviceMDData_part3 (IndexType numAtom,
						       const TypeType * type1,
						       const ScalorType * mass1,
						       const ScalorType * massi1,
						       const ScalorType * charge1,
						       const IndexType * idxTable,
						       TypeType * type2,
						       ScalorType * mass2,
						       ScalorType * massi2,
						       ScalorType * charge2);





#endif
