class NeighborList;

#ifndef __NeighborList_interface_h_wanghan__
#define __NeighborList_interface_h_wanghan__

#include "common.h"
#include "BoxGeometry.h"
#include "Reshufflable.h"
#include "MDSystem.h"
#include "MDSystem_interface.h"
#include "NeighborList.h"
#include "MDError_interface.h"
#include "MDTimer_interface.h"
#include "SystemNonBondedInteraction.h"

#include "CellList_interface.h"

using namespace RectangularBoxGeometry;

class NeighborList : public Reshufflable
{
  MDError err;
private:
  // IndexType buildDeviceNeighborList_AllPair_sbuffSize;
  // IndexType buildDeviceNeighborList_DeviceCellList_sbuffSize;
private:
  bool mallocedDeviceNeighborList;
  bool mallocedNonBondedForceTable;
  bool initedGlobalTexture;
  void clearDeviceNeighborList();
  void clearNonBondedForce ();
  void unbindGlobalTexture ();
  void clear();
  void mallocDeviceNeighborList (const MDSystem & sys,
				 const ScalorType & DeviceNeighborListExpansion);
  void initNonBondedInteraction (const SystemNonBondedInteraction & sysNbInter);
  void bindGlobalTexture (const MDSystem & sys);
private:
  void normalizeMDSystem (const MDSystem & sys);
  void buildDeviceNeighborListCellList (const MDSystem & sys,
					const CellList & clist);
  void buildDeviceNeighborListAllPair  (const MDSystem & sys);

public:
  ScalorType		myrlist;
  DeviceNeighborList	dnlist;
  dim3			atomGridDim;
  dim3			myBlockDim;
  IndexType		NatomType;
  IndexType *		nbForceTable;
  IndexType		nbForceTableLength;
private: // reshuffle backup
  IndexType *		bkdnlistData;
  IndexType *		bkdnlistNneighbor;
  IndexType *		bkdnlistForceIndex;
  CoordType *		bkbackupCoord;
public:
  NeighborList (const SystemNonBondedInteraction & sysNbInter,
		const MDSystem & sys,
		const ScalorType & rlist,
		const IndexType & NTread,
		const ScalorType & DeviceNeighborListExpansion = 1);
  ~NeighborList();
  void reinit (const SystemNonBondedInteraction & sysNbInter,
	       const MDSystem & sys,
	       const ScalorType & rlist,
	       const IndexType & NTread,
	       const ScalorType & DeviceNeighborListExpansion = 1);
  void rebuild (const MDSystem & sys,
		const CellList & clist,
		MDTimer * timer = NULL);
  virtual void reshuffle (const IndexType * indexTable,
			  const IndexType & numAtom,
			  MDTimer *timer = NULL);
}
    ;




// typedef enum {
//   NullMode, 
//   AllPairBuilt, 
//   CellListBuilt
// } NeighborListBuiltMode;

// class WrongBuildMethod {};

// class NeighborList : public Reshufflable
// {
//   IndexType * mySendBuff;
//   IndexType * myTargetBuff;
//   // ScalorType * judgeRebuild_buff;
//   IndexType  * judgeRebuild_tag;
//   SumVector<IndexType> sum_judge;
//   MDError err;
// private:
//   IndexType buildDeviceNeighborList_AllPair_sbuffSize;
//   IndexType buildDeviceNeighborList_DeviceCellList_sbuffSize;
//   IndexType buildDeviceCellList_step1_sbuffSize;
//   IndexType judgeRebuild_judgeCoord_block_sbuffSize;
// private:
//   bool mallocedDeviceCellList;
//   bool mallocedDeviceNeighborList;
//   bool initedGlobalTexture;
//   bool mallocedJudgeStuff;
//   bool mallocedNonBondedForceTable;
//   void clearDeviceCellList();
//   void clearDeviceNeighborList();
//   void unbindGlobalTexture();
//   void clearJudgeStuff();
//   void clear();
//   void DecideNeighboringMethod (const MDSystem & sys,
// 				const ScalorType & rlist,
// 				const BoxDirection_t & bdir,
// 				const IndexType & divide,
// 				NeighborListBuiltMode & mode,
// 				IntVectorType & NCell);
//   void mallocDeviceCellList (const IntVectorType & NCell,
// 			     const HostVectorType & boxSize);
//   void mallocDeviceNeighborList (const MDSystem & sys,
// 				 const ScalorType & DeviceNeighborListExpansion);
//   void bindGlobalTexture (const MDSystem & sys);
//   void mallocJudgeStuff(const MDSystem & sys);
//   void initNonBondedInteraction (const SystemNonBondedInteraction & sysNbInter);
// private:
//   void normalizeMDSystem (const MDSystem & sys);
//   void backupJudgeCoord (const MDSystem & sys);
//   void naivelyBuildDeviceCellList (const MDSystem & sys);
//   void buildDeviceCellList (const MDSystem & sys);
//   void buildDeviceNeighborListCellList (const MDSystem & sys);
//   void buildDeviceNeighborListAllPair  (const MDSystem & sys);
// public:
//   NeighborListBuiltMode mode;
//   ScalorType myrlist;
//   BoxDirection_t mybdir;
//   IndexType mydivide;
//   DeviceCellList dclist;
//   DeviceNeighborList dnlist;
//   IndexType NatomType;
//   IndexType * nbForceTable;
//   IndexType nbForceTableLength;
//   bool sharednbForceTable;
//   IndexType bitDeepth;
//   dim3 cellGridDim;
//   dim3 atomGridDim;
//   dim3 myBlockDim;
// public:
//   CoordType * backupCoord;
// private: // reshuffle backup
//   IndexType * bkdnlistData;
//   IndexType * bkdnlistNneighbor;
//   IndexType * bkdnlistForceIndex;
//   CoordType * bkbackupCoord;
// public:
//   /** 
//    * Constructure
//    * 
//    * @param sys System on which the neighbor list will be built.
//    * @param rlist The radius of neighbor list
//    * @param NTread Number of thread per block. All calls to kernel will
//    * use this the same value.
//    * @param DeviceNeighborListExpansion The program will automaticaly calculate
//    * the length of neighbor list. But some times this length is not long enough
//    * to maintain all neighbors (especially for the non-homogeneous systems).
//    * This user provided parameter is the ratio of actually length over the
//    * calculated length.
//    * @param bdir Denotes which direction of box will be cell divided. User
//    * should validate this parameter.
//    */
//   NeighborList (const SystemNonBondedInteraction & sysNbInter,
// 		const MDSystem & sys,
// 		const ScalorType & rlist, const IndexType & NTread,
// 		const ScalorType & DeviceNeighborListExpansion = 5,
// 		const RectangularBoxGeometry::BoxDirection_t & bdir = 7,
// 		const IndexType divide = 1)
//       : mallocedDeviceCellList (false), mallocedDeviceNeighborList (false),
// 	mallocedJudgeStuff (false), initedGlobalTexture (false),
// 	mallocedNonBondedForceTable (false),
// 	mode (NullMode)
//       {init (sysNbInter, sys, rlist, NTread, DeviceNeighborListExpansion, bdir, divide);}
//   /** 
//    * Destructor
//    * 
//    */
//   ~NeighborList ();
//   /** 
//    * Initialize the neighbor list
//    * 
//    * @param sys System on which the neighbor list will be built.
//    * @param rlist The radius of neighbor list
//    * @param NTread Number of thread per block. All calls to kernel will
//    * use this the same value.
//    * @param DeviceNeighborListExpansion The program will automaticaly calculate
//    * the length of neighbor list. But some times this length is not long enough
//    * to maintain all neighbors (especially for the non-homogeneous systems).
//    * This user provided parameter is the ratio of actually length over the
//    * calculated length.
//    * @param bdir Denotes which direction of box will be cell divided. User
//    * should validate this parameter.
//    */
//   void init (const SystemNonBondedInteraction & sysNbInter,
// 	     const MDSystem & sys,
// 	     const ScalorType & rlist,
// 	     const IndexType & NTread,
// 	     const ScalorType & DeviceNeighborListExpansion = 5,
// 	     const RectangularBoxGeometry::BoxDirection_t & bdir = 7,
// 	     const IndexType & cellListDivideLevel = 1);
//   void reinit (const MDSystem & sys,
// 	       const ScalorType & rlist,
// 	       const IndexType & NTread,
// 	       const ScalorType & DeviceNeighborListExpansion = 5,
// 	       const RectangularBoxGeometry::BoxDirection_t & bdir = 7,
// 	       const IndexType & cellListDivideLevel = 1);
//   /** 
//    * Build neighbor list
//    * 
//    * @param sys System on which the neighbor list will be built.
//    * @param timer A timer monitors the time consumption.
//    */
//   void build (const MDSystem & sys,
// 	      MDTimer *timer = NULL);
//   /** 
//    * Rebuild neighbor list
//    * 
//    * @param sys System on which the neighbor list will be built.
//    * @param timer A timer monitors the time consumption.
//    */
//   void reBuild (const MDSystem & sys,
// 		MDTimer *timer = NULL);
//   /** 
//    * Judge if the neighbor should be rebuild. This function will not actually
//    * rebuild the list.
//    * 
//    * @param sys System on which the neighbor list will be built.
//    * @param diffTol If the displacement of any atom hits this tolrance, the
//    * the program will suggest you to rebuild the neighbor list.
//    * @param timer A timer monitors the time consumption.
//    * 
//    * @return ture if it suggest you to rebuild the neighbor list. false otherwise.
//    */
//   bool judgeRebuild (const MDSystem & sys,
// 		     const ScalorType &diffTol,
// 		     MDTimer *timer=NULL);
//   virtual void reshuffle (const IndexType * indexTable,
// 			  const IndexType & numAtom,
// 			  MDTimer *timer=NULL);
//   void reshuffleCell (const IndexType * indexTable,
// 		      const IndexType & numAtom,
// 		      MDTimer *timer=NULL);
//   void buildCellList (const MDSystem & sys,
// 		      MDTimer * timer = NULL);
//   void reBuildCellList (const MDSystem & sys,
// 			MDTimer * timer = NULL);
// };

	     

#endif
