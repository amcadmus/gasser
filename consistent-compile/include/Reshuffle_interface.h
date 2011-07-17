#ifndef __Reshuffle_interface_h_wanghan__
#define __Reshuffle_interface_h_wanghan__

#include "common.h"
#include "MDSystem_interface.h"
#include "NeighborList_interface.h"
// #include "MDTimer.h"

template<typename TYPE>
static __global__ void
Reshuffle_reshuffleArray (const TYPE * bkbuff,
			  const IndexType numAtom,
			  const IndexType * idxTable,
			  TYPE * buff)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  if (ii < numAtom){
    buff[idxTable[ii]] = bkbuff[ii];
  }
}

/// Calculate the index table according to a cell list.
/**
 * The index table is used to reshuffle objects that are
 * Reshufflable. The table tells objects how the indexes of atoms
 * mapped. After reshuffle, atoms in one cell are stored continuously.
 * The cells are stored in a natural 3D index manner.
 */

class Reshuffle
{
  bool malloced;
  void clear ();
public:
  IndexType * indexTable;	/**< The index table, on device. This
				 * is a public data member because of
				 * the restriction of CUDA. */
public:
  /** 
   * Constructor.
   * 
   * @param sys The MDSystem to reshuffle.
   */
  Reshuffle (const MDSystem & sys);
  ~Reshuffle ();
  /** 
   * Reinitializer.
   * 
   * @param sys The MDSystem to reshuffle.
   */
  void reinit (const MDSystem & sys);
  /** 
   * Calculate the index table according to the cell list.
   * 
   * @param clist The cell list.
   * @param timer Timer measuring the performance.
   * 
   * @return ture if the index table is successfully built. false
   * otherwise.
   */
  bool calIndexTable (const CellList & clist,
		      MDTimer * timer = NULL);
  // const IndexType * getIndexTable () const {return indexTable;}
};


// class Reshuffle
// {
//   IndexType * posiBuff;
//   IndexType * idxTable;
//   DeviceMDData imageSystem;
//   IndexType * bknlistData;
//   ForceIndexType * bkNBForceIndex;
//   IndexType * bkNneighbor;
//   IndexType * backMapTable;
//   IndexType * backMapTableBuff;
//   // bk bond list
//   IndexType * bkBondListData;
//   ForceIndexType * bkBondListBondIndex;
//   IndexType * bkBondListNumB;
//   // bk angle list
//   IndexType * bkAngleListNei;
//   IndexType * bkAngleListPosi;
//   ForceIndexType * bkAngleListAngleIndex;
//   IndexType * bkAngleListNangle;
// #ifndef COORD_IN_ONE_VEC
//   ScalorType * bkNlistJudgeBuffx;
//   ScalorType * bkNlistJudgeBuffy;
//   ScalorType * bkNlistJudgeBuffz;
// #else
//   CoordType * bkNlistJudgeBuff;
// #endif
//   dim3 cellGridDim;
//   dim3 atomGridDim;
//   dim3 myBlockDim;
//   void recoverMDData (const DeviceMDData & currentData,
// 		      DeviceMDData & recoveredData,
// 		      MDTimer * timer);
//   bool hasBond ;
//   bool hasAngle;
// public:
//   Reshuffle (const MDSystem & sys,
// 	     const NeighborList & nlist, 
// 	     const IndexType & NTread) 
//       {init (sys, nlist, NTread);}
//   void init (const MDSystem & sys,
// 	     const NeighborList & nlist,
// 	     const IndexType & NTread);
//   ~Reshuffle ();
//   void shuffleSystem (MDSystem & sys,
// 		      NeighborList & nlist,
// 		      MDTimer *timer=NULL);
//   void recoverMDDataToHost (MDSystem & sys,
// 			    MDTimer *timer=NULL);  
// };





#endif
