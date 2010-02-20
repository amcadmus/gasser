#ifndef __NeighborList_h_wanghan__
#define __NeighborList_h_wanghan__

#include "common.h"
#include "MDSystem.h"
#include "BoxGeometry.h"
#include "Auxiliary.h"

using namespace RectangularBoxGeometry;

struct DeviceCellList 
{
  ScalorType rlist;		/**< radius for list building */
  IntVectorType NCell;		/**< the number of cells on each direction */
  VectorType NCelli;		/**< inverse the NCell */
  IndexType * data;		/**< matrix, the i-th line keep indexes of
				 * atom in the i-th cell */
  IndexType * numbers;		/**< vector stores the number of atoms
				 * in each cell */
  IndexType stride;		/**< stride of cell list, which should
				 * be larger than thread per block */
};

struct DeviceNeighborList
{
  ScalorType rlist;		/**< radius of neighbor list */
  IndexType listLength;		/**< max length of the neighbor list */
  IndexType * data;		/**< matrix, i-th row stores the
				 * indexes of neigbors of i-th atom*/
  IndexType * Nneighbor;	/**< vector stores the number of neigbors of
				 * each atom*/
  ForceIndexType * forceIndex;	/**< matrix stores the index of
				 * non-bonded interaction, whose type
				 * and parameters are keeped by the
				 * interaction engine. */
  IndexType stride;             /**< stride of neighbor list, which is
				 * the expected larger than or equal
				 * to the number of Atoms */
};

/** 
 * An 3D index to 1D index mapping for cells.
 * 
 * @param clist Reference to the cell list.
 * @param ix x component of cell 3D index.
 * @param iy y component of cell 3D index.
 * @param iz z component of cell 3D index.
 * 
 * @return 1D index of the corresponding cell
 */
__device__ IndexType D3toD1 (const DeviceCellList &clist, 
			     const IndexType &ix,
			     const IndexType &iy,
			     const IndexType &iz);
/** 
 * An 1D index to 3D index mapping for cells.
 * 
 * @param clist Reference of the cell list
 * @param i 1D index of the cell
 * @param x x component of cell 3D index.
 * @param y y component of cell 3D index.
 * @param z z component of cell 3D index.
 */
__device__ void D1toD3 (const DeviceCellList & clist,
			const IndexType &i, 
			IndexType &x,
			IndexType &y,
			IndexType &z);

__device__ IndexType getDeviceCellListData (const DeviceCellList & clist, 
					    const IndexType & cid,
					    const IndexType & aid);

////////////////////////////////////////////////////////////
// cell list operations
////////////////////////////////////////////////////////////
// needs NCell blocks
/** 
 * Prepare the building of device cell list. Needs at least NCell
 * blocks.
 * 
 * @param clist the cell list
 */
__global__ void prepare_naivelyBuildDeviceCellList (DeviceCellList  clist);
// needs ceil(numAtom/blockDim.x) blocks
__global__ void naivelyBuildDeviceCellList (IndexType numAtom,
#ifndef COORD_IN_ONE_VEC
					    ScalorType * coordx,
					    ScalorType * coordy,
					    ScalorType * coordz,
#else
					    CoordType * coord,
#endif
					    RectangularBox box,
					    DeviceCellList clist,
					    mdError_t * ptr_de = NULL,
					    IndexType * erridx = NULL,
					    ScalorType * errsrc = NULL);
__global__ void naivelyBuildDeviceCellList2 (IndexType numAtom,
#ifndef COORD_IN_ONE_VEC
					     ScalorType * coordx,
					     ScalorType * coordy,
					     ScalorType * coordz,
#else
					     CoordType * coord,
#endif
					     IntScalorType * coordNoix,
					     IntScalorType * coordNoiy,
					     IntScalorType * coordNoiz,
					     RectangularBox box,
					     DeviceCellList clist,
					     mdError_t * ptr_de = NULL,
					     IndexType * erridx = NULL,
					     ScalorType * errsrc = NULL);

// needs NCell blocks
__global__ void buildDeviceCellList_initBuff (IndexType * sendBuff,
					      IndexType * targetBuff);
// needs NCell blocks
__global__ void buildDeviceCellList_clearBuff (IndexType * sendBuff);
// needs NCell blocks
__global__ void buildDeviceCellList_step1 (IndexType numAtom,
#ifndef COORD_IN_ONE_VEC
					   ScalorType * coordx,
					   ScalorType * coordy,
					   ScalorType * coordz,
#else
					   CoordType * coord,
#endif
					   IntScalorType * coordNoix,
					   IntScalorType * coordNoiy,
					   IntScalorType * coordNoiz,
					   RectangularBox box,
					   DeviceCellList clist,
					   IndexType * sendBuff,
					   IndexType * targetBuff,
					   mdError_t * ptr_de = NULL,
					   IndexType * erridx = NULL,
					   ScalorType * errsrc = NULL);

// needs NCell blocks
__global__ void buildDeviceCellList_step2 (RectangularBox box,
					   DeviceCellList clist,
					   IndexType * sendBuff,
					   IndexType * targetBuff,
					   IndexType bitDeepth,
					   mdError_t * ptr_de = NULL);

////////////////////////////////////////////////////////////
// neighbor list operations
////////////////////////////////////////////////////////////
//needs ceil(numAtom/blockDim.x) blocks
__global__ void buildDeviceNeighborList_AllPair  (IndexType numAtom,
#ifndef COORD_IN_ONE_VEC
						  ScalorType * coordx,
						  ScalorType * coordy, 
						  ScalorType * coordz,
#else
						  CoordType * coord,
#endif
						  TypeType * type,
						  RectangularBox box,
						  DeviceNeighborList nlist,
						  ForceIndexType * nbForceTable,
						  IndexType NatomType,
						  bool sharednbForceTable,
						  mdError_t * ptr_de = NULL);
// needs NCell blocks
__global__ void buildDeviceNeighborList_DeviceCellList (IndexType numAtom,
#ifndef COORD_IN_ONE_VEC
							ScalorType * coordx, 
							ScalorType * coordy, 
							ScalorType * coordz,
#else
							CoordType * coord,
#endif
							TypeType * type,
							RectangularBox box,
							DeviceCellList clist,
							DeviceNeighborList nlist,
							ForceIndexType * nbForceTable,
							IndexType NatomType,
							bool sharednbForceTable,
							mdError_t * ptr_de = NULL);
// needs ceil(numAtom/blockDim.x) blocks
__global__ void judgeRebuild_backupCoord (const IndexType numAtom,
					  const ScalorType * from,
					  ScalorType * to);
// needs ceil(numAtom/blockDim.x) blocks
__global__ void judgeRebuild_judgeCoord (const IndexType numAtom,
					 const RectangularBox box,
#ifndef COORD_IN_ONE_VEC
					 const ScalorType* coordx,
					 const ScalorType* coordy,
					 const ScalorType* coordz,
					 const ScalorType* backupCoordx,
					 const ScalorType* backupCoordy,
					 const ScalorType* backupCoordz,
#else
					 const CoordType * coord,
					 const CoordType * backupCoord,
#endif
					 const ScalorType diffTol2,
					 ScalorType * buff,
					 IndexType * tag);
__global__ void judgeRebuild_judgeCoord_block (const IndexType numAtom,
					       const RectangularBox box,
#ifndef COORD_IN_ONE_VEC
					       const ScalorType* coordx,
					       const ScalorType* coordy,
					       const ScalorType* coordz,
					       const ScalorType* backupCoordx,
					       const ScalorType* backupCoordy,
					       const ScalorType* backupCoordz,
#else
					       const CoordType * coord,
					       const CoordType * backupCoord,
#endif
					       const ScalorType diffTol2,
					       IndexType * judgeRebuild_buff);


////////////////////////////////////////////////////////////
// intialize cell
////////////////////////////////////////////////////////////

// thread 0 in block 0 makes allocate for the cell list. in the end
// other blocks may not see the allocated cell list.
// __device__ void initDeviceCellList (const DeviceMDData & ddata,
// 				    const RectangularBoxGeometry::RectangularBox &box,
// 				    const ScalorType &rlist,
// 				    DeviceCellList * clist,
// 				    IndexType * suggestGridDim)
// {  
//   IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
//   IndexType tid = threadIdx.x;
//   IndexType numCell;
  
//   // the cell list dimension is exactly the block dimension
//   if (tid + bid == 0){
//     double rlisti = 1./rlist;
//     clist->NCelli.x = 1./ (clist->NCell.x = floorf(box.size.x * rlisti));
//     clist->NCelli.y = 1./ (clist->NCell.y = floorf(box.size.y * rlisti));
//     clist->NCelli.z = 1./ (clist->NCell.z = floorf(box.size.z * rlisti));
//     numCell = clist->NCell.x * clist->NCell.y * clist->NCell.z;
//     clist->rlist = rlist;
//     // suppose the number of atoms in any cell is smaller or equal
//     // to the number of threads in a block
//     clist->stride = blockDim.x;
// #ifdef COMPILE_EMU
//     clist->data = (IndexType *)malloc(sizeof(ScalorType) * numCell * clist->stride);
//     clist->numbers = (unsigned *)malloc( sizeof(unsigned) * numCell);
// #else
//     cudaMalloc ((void**)&(clist->data), sizeof(ScalorType) * numCell * clist->stride);
//     //          size                 Number of cell    stide
//     cudaMalloc ((void**)&(numbers), sizeof(unsigned) * numCell);
// #endif
//     for (IndexType i = 0; i < numCell; ++i)
//       clist->numbers[i] = 0;
//     * suggestGridDim = numCell;
//   }
//   // printf ("1 here bid is %d thid is %d\n", bid, tid);
// }

  



////////////////////////////////////////////////////////////
// implementation
// cell list operations
////////////////////////////////////////////////////////////

__device__ IndexType getDeviceCellListData (const DeviceCellList & clist, 
					    const IndexType & cid, const IndexType & pid)
{
  return clist.data[cid * clist.stride + pid];
}

__device__ IndexType D3toD1 (const DeviceCellList & clist,
			     const IndexType &ix, const IndexType &iy, const IndexType &iz)
{
  return iz + clist.NCell.z * iy + clist.NCell.z * clist.NCell.y * ix;
  // return IndexType(clist.NCell.y) * (IndexType(clist.NCell.x) * ix + iy) + iz;
}
__device__ void D1toD3 (const DeviceCellList & clist, const IndexType &i, 
			IndexType &x, IndexType &y, IndexType &z)
{
  IndexType tmp = i;
  z = tmp % (clist.NCell.z);
  tmp = (tmp - z) / clist.NCell.z;
  y = tmp % (clist.NCell.y);
  x = (tmp - y) / clist.NCell.y;
}



////////////////////////////////////////////////////////////
// neighbor list operations
////////////////////////////////////////////////////////////

// __device__ void initDeviceNeighborList (const DeviceMDData & ddata,
// 					const RectangularBoxGeometry::RectangularBox &box,
// 					const ScalorType & rlist,
// 					DeviceNeighborList * nlist,
// 					IndexType DeviceNeighborListExpansion)
// {
//   IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
//   if (threadIdx.x + bid == 0){
//     nlist->rlist = rlist;
//     nlist->stride = ddata.numAtom;
//     ScalorType density = ddata.numAtom / (box.size.x * box.size.y * box.size.z);
//     IndexType expectedNumberInList 
// 	= 4./3. * M_PI * rlist * rlist * rlist * density;
//     nlist->listLength = expectedNumberInList * DeviceNeighborListExpansion;
// #ifdef COMPILE_EMU    
//     nlist->data = (IndexType *)malloc (sizeof(IndexType) * nlist->stride * nlist->listLength);
//     nlist->Nneighbor = (IndexType *)malloc(sizeof(IndexType) * ddata.numAtom);
// #else
//     cudaMalloc ((void**)&(nlist->data), sizeof(IndexType) * nlist->stride * nlist->listLength);
//     cudaMalloc ((void**)&(nlist->Nneighbor), sizeof(IndexType) * ddata.numAtom);
// #endif
//   }
// }



__device__ void kthSort (volatile IndexType * sharedbuff, IndexType k)
{
  IndexType tid = threadIdx.x;
  __shared__ volatile IndexType predbuff[MaxThreadsPerBlock * 2];
  predbuff[tid] = getKthBit(sharedbuff[tid], k);
  predbuff[tid+blockDim.x] = 0;

  __syncthreads();
  IndexType total1 = sumVectorBlockBuffer (predbuff, blockDim.x);
  IndexType target, mydata = sharedbuff[tid];
  if (getKthBit(sharedbuff[tid], k)) {
    target = blockDim.x - predbuff[tid];
  }
  else {
    // IndexType total0 = blockDim.x - total1;
    // IndexType after0 = blockDim.x - tid - predbuff[tid];
    // target = total0 - after0;
    target = tid + predbuff[tid] - total1;
  }
  __syncthreads();
  sharedbuff[target] = mydata;
  __syncthreads();
}

__device__ IndexType headSort (volatile IndexType * sharedbuff,
			       volatile IndexType * targetCell)
{
  IndexType k = NUintBit - 1;
  IndexType tid = threadIdx.x;
  __shared__ volatile IndexType predbuff[MaxThreadsPerBlock * 2];
  predbuff[tid] = getKthBit(sharedbuff[tid], k);
  predbuff[tid+blockDim.x] = 0;

  __syncthreads();
  IndexType total1 = sumVectorBlockBuffer (predbuff, blockDim.x);
  IndexType target, mydata = sharedbuff[tid], mycell = targetCell[tid];
  if (getKthBit(sharedbuff[tid], k)) {
    target = blockDim.x - predbuff[tid];
  }
  else {
    target = tid + predbuff[tid] - total1;
  }
  __syncthreads();
  sharedbuff[target] = mydata;
  targetCell[target] = mycell;
  __syncthreads();
  return total1;
}


__device__ void sortList (volatile IndexType * sharedbuff, IndexType bitDeepth)
{
  __syncthreads();
  for (IndexType i = 0; i < bitDeepth; ++i){
    kthSort(sharedbuff, i);
  }
}




#endif
