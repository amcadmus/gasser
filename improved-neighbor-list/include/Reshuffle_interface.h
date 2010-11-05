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


class Reshuffle
{
  bool malloced;
  void clear ();
public:
  IndexType * indexTable;
public:
  Reshuffle (const MDSystem & sys);
  ~Reshuffle ();
  void reinit (const MDSystem & sys);
  bool calIndexTable (const CellList & nlist,
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
