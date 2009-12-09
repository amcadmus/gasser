#ifndef __Reshuffle_interface_h_wanghan__
#define __Reshuffle_interface_h_wanghan__

#include "common.h"
#include "MDSystem_interface.h"
#include "NeighborList_interface.h"
#include "BondList_interface.h"
#include "MDTimer.h"

class Reshuffle
{
  IndexType * posiBuff;
  IndexType * idxTable;
  DeviceMDData imageSystem;
  IndexType * bknlistData;
  ForceIndexType * bkNBForceIndex;
  IndexType * bkNneighbor;
  IndexType * backMapTable;
  IndexType * backMapTableBuff;
  IndexType * bkBondListData;
  ForceIndexType * bkBondListBondIndex;
  IndexType * bkBondListNumB;
#ifndef COORD_IN_ONE_VEC
  ScalorType * bkNlistJudgeBuffx;
  ScalorType * bkNlistJudgeBuffy;
  ScalorType * bkNlistJudgeBuffz;
#else
  CoordType * bkNlistJudgeBuff;
#endif
  dim3 cellGridDim;
  dim3 atomGridDim;
  dim3 myBlockDim;
  void recoverMDData (const DeviceMDData & currentData,
		      DeviceMDData & recoveredData,
		      MDTimer * timer);
public:
  Reshuffle (const MDSystem & sys,
	     const NeighborList & nlist, 
	     const IndexType & NTread) 
      {init (sys, nlist, NTread);}
  void init (const MDSystem & sys,
	     const NeighborList & nlist,
	     const IndexType & NTread);
  ~Reshuffle ();
  void shuffleSystem (MDSystem & sys,
		      NeighborList & nlist,
		      MDTimer *timer=NULL);
  void recoverMDDataToHost (MDSystem & sys,
			    MDTimer *timer=NULL);  
};





#endif
