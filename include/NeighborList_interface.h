#ifndef __NeighborList_interface_h_wanghan__
#define __NeighborList_interface_h_wanghan__

#include "common.h"
#include "BoxGeometry.h"
#include "MDSystem.h"
#include "NeighborList.h"
#include "MDError_interface.h"
#include "MDTimer_interface.h"
#include "MDSystem_interface.h"

using namespace RectangularBoxGeometry;

typedef enum {
  NullMode, 
  AllPairBuilt, 
  CellListBuilt
} NeighborListBuiltMode;


class NeighborList
{
  bool malloced ;
  IndexType * mySendBuff;
  IndexType * myTargetBuff;
  // ScalorType * judgeRebuild_buff;
  IndexType  * judgeRebuild_tag;
  SumVector<IndexType> sum_judge;
  MDError err;
private:
  IndexType buildDeviceNeighborList_AllPair_sbuffSize;
  IndexType buildDeviceNeighborList_DeviceCellList_sbuffSize;
  IndexType buildDeviceCellList_step1_sbuffSize;
  IndexType judgeRebuild_judgeCoord_block_sbuffSize;
public:
  NeighborListBuiltMode mode;
  bool CellOnX, CellOnY, CellOnZ;
  DeviceCellList dclist;
  DeviceNeighborList dnlist;
  IndexType NatomType;
  ForceIndexType * nbForceTable;
  IndexType nbForceTableLength;
  bool sharednbForceTable;
  IndexType bitDeepth;
  dim3 cellGridDim;
  dim3 atomGridDim;
  dim3 myBlockDim;
public:
  ScalorType * backupCoordx;
  ScalorType * backupCoordy;
  ScalorType * backupCoordz;
public:
  /** 
   * Constructure
   * 
   * @param sys System on which the neighbor list will be built.
   * @param rlist The radius of neighbor list
   * @param NTread Number of thread per block. All calls to kernel will
   * use this the same value.
   * @param DeviceNeighborListExpansion The program will automaticaly calculate
   * the length of neighbor list. But some times this length is not long enough
   * to maintain all neighbors (especially for the non-homogeneous systems).
   * This user provided parameter is the ratio of actually length over the
   * calculated length.
   * @param bdir Denotes which direction of box will be cell devided. User
   * should validate this parameter.
   */
  NeighborList (const MDSystem & sys,
		const ScalorType & rlist, const IndexType & NTread,
		const IndexType & DeviceNeighborListExpansion = 5,
		const RectangularBoxGeometry::BoxDirection_t & bdir = 7)
      {init (sys, rlist, NTread, DeviceNeighborListExpansion, bdir);}
  /** 
   * Destructor
   * 
   */
  ~NeighborList ();
  /** 
   * Initialize the neighbor list
   * 
   * @param sys System on which the neighbor list will be built.
   * @param rlist The radius of neighbor list
   * @param NTread Number of thread per block. All calls to kernel will
   * use this the same value.
   * @param DeviceNeighborListExpansion The program will automaticaly calculate
   * the length of neighbor list. But some times this length is not long enough
   * to maintain all neighbors (especially for the non-homogeneous systems).
   * This user provided parameter is the ratio of actually length over the
   * calculated length.
   * @param bdir Denotes which direction of box will be cell devided. User
   * should validate this parameter.
   */
  void init (const MDSystem & sys,
	     const ScalorType & rlist, const IndexType & NTread,
	     const IndexType & DeviceNeighborListExpansion = 5,
	     const RectangularBoxGeometry::BoxDirection_t & bdir = 7);
  /** 
   * Build neighbor list
   * 
   * @param sys System on which the neighbor list will be built.
   * @param timer A timer monitors the time consumption.
   */
  void build (const MDSystem & sys,
	      MDTimer *timer = NULL);
  /** 
   * Rebuild neighbor list
   * 
   * @param sys System on which the neighbor list will be built.
   * @param timer A timer monitors the time consumption.
   */
  void reBuild (const MDSystem & sys,
		MDTimer *timer = NULL);
  /** 
   * Judge if the neighbor should be rebuild. This function will not actually
   * rebuild the list.
   * 
   * @param sys System on which the neighbor list will be built.
   * @param diffTol If the displacement of any atom hits this tolrance, the
   * the program will suggest you to rebuild the neighbor list.
   * @param timer A timer monitors the time consumption.
   * 
   * @return ture if it suggest you to rebuild the neighbor list. false otherwise.
   */
  bool judgeRebuild (const MDSystem & sys,
		     const ScalorType &diffTol,
		     MDTimer *timer=NULL);
};

	     

#endif
