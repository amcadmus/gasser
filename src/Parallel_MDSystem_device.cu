#define DEVICE_CODE
#include "Parallel_Timer.h"
#include "Parallel_MDSystem.h"
#include "Parallel_Interface.h"
#include "Parallel_BondList.h"

#include "compile_error_mixcode.h"

using namespace Parallel::Timer;

Parallel::MDSystem::
MDSystem ()
    : atomName(NULL), atomIndex(NULL),
      resdName(NULL), resdIndex(NULL),
      numFreedom (0)
{
}

Parallel::MDSystem::
~MDSystem()
{
  freeAPointer ((void**)&atomName);
  freeAPointer ((void**)&atomIndex);
  freeAPointer ((void**)&resdName);
  freeAPointer ((void**)&resdIndex);
}


void Parallel::MDSystem::
reallocGroProperty (const IndexType & memSize__)
{
  if (memSize__ == 0) return ;
  IndexType memSize_ = memSize__;

  size_t sizec = memSize_ * sizeof(char) * StringSize;
  size_t sizeidx = memSize_ * sizeof(IndexType);
  
  atomName = (char *) realloc (atomName, sizec);
  if (atomName == NULL) throw MDExcptFailedReallocOnHost("atomName", sizec);
  atomIndex = (IndexType *) realloc (atomIndex, sizeidx);
  if (atomIndex == NULL) throw MDExcptFailedReallocOnHost("atomIndex", sizeidx);
  resdName = (char *) realloc (resdName, sizec);
  if (resdName == NULL) throw MDExcptFailedReallocOnHost("resdName", sizec);
  resdIndex = (IndexType *) realloc (resdIndex, sizeidx); 
  if (resdIndex == NULL) throw MDExcptFailedReallocOnHost("resdIndex", sizeidx);
}



void Parallel::MDSystem::
init (const char * confFileName,
      const Topology::System & sysTop,
      const ScalorType & cellSize,
      const IndexType & divideLevel)
{
  numFreedom = sysTop.getNumFreedom();
  MDDataItemMask_t maskConfig = (MDDataItemMask_Coordinate |
  				 MDDataItemMask_CoordinateNoi |
  				 MDDataItemMask_Velocity |
  				 MDDataItemMask_GlobalIndex);
  if (Parallel::Interface::myRank() == 0){
    globalNumAtom = globalHostData.numAtomInGroFile (confFileName);
    // globalHostData.reallocCoordinate (globalNumAtom);
    // globalHostData.reallocVelocity   (globalNumAtom);
    // globalHostData.reallocGroProperty(globalNumAtom);
    // globalHostData.reallocGlobalIndex(globalNumAtom);
    globalHostData.easyMalloc (globalNumAtom);
    reallocGroProperty (globalNumAtom);
    globalHostData.initConf_GroFile (confFileName,
				     atomName, atomIndex,
				     resdName, resdIndex);
    for (int i = 0; i < globalNumAtom; ++i){
      hostMoveParticleToBox (globalHostData.getGlobalBox(),
			     & ((globalHostData.cptr_coordinate()[i]).x), 
			     & ((globalHostData.cptr_coordinate()[i]).y), 
			     & ((globalHostData.cptr_coordinate()[i]).z),
			     & (globalHostData.cptr_coordinateNoi()[i].x),
			     & (globalHostData.cptr_coordinateNoi()[i].y),
			     & (globalHostData.cptr_coordinateNoi()[i].z));
    }
    // globalHostData.initTopology_property (sysTop);
    for (int i = 0; i < globalNumAtom; ++i){
      globalHostData.cptr_coordinate()[i].w = globalHostData.cptr_type()[i];
    }
  }
  
  distributeGlobalMDData (globalHostData, localHostData);

  DeviceMDData & ddata (deviceData);
  ddata.copyFromHost (localHostData, MDDataItemMask_All);

  deviceData.initCellStructure (cellSize, divideLevel);
  // printf ("ncell: %d\n", deviceData.getNumCell().x);
  cellRelation.rebuild (deviceData);
  
  // for (IndexType i = 0; i < deviceData.numData(); ++i){
  //   deviceData.dptr_coordinate()[i].x += 2;
  //   deviceData.dptr_coordinate()[i].y += 2;
  //   deviceData.dptr_coordinate()[i].z += 2;
  // }  
  // deviceData.rebuild ();
  // printf ("rank %d\n", Parallel::Interface::myRank());

  deviceData.copyToHost (localHostData, maskConfig);
  SystemBondedInteraction sysBdInter (sysTop);
  localHostData.initTopology (sysTop, sysBdInter);
  deviceData.copyFromHost (localHostData, MDDataItemMask_All);

  // MDDataItemMask_t maskWithBond = maskConfig | MDDataItemMask_Bond;
  // MDDataItemMask_t maskWithBondAngle = maskConfig | MDDataItemMask_Bond | MDDataItemMask_Angle;
  // HostMDData ht1 ;
  // HostMDData ht2 ;
  // DeviceMDData dt1;  
  // ht1.easyMalloc (16, 7, 0, 0);
  // ht1.numData() = 3;
  // dt1.easyMalloc (14, 5, 0, 0);
  // dt1.numData() = 5;
  // dt1.copyToHost (ht1, maskWithBond);
  // // ht2.easyMalloc (20, 1, 1, 1);
  // // ht2.numData() = 10;
  // // ht1.copy (ht2, maskWithBond);
  // // ht2.copy (ht1, maskWithBond);

  // deviceData.dptr_coordinate()[124].x = 3;
  // deviceData.dptr_coordinate()[126].x = 4;
  // SubCellList subList;
  // localHostData.buildSubList (1, 2, 1, 2, 1, 3, subList);
  // HostSubCellList hsub ;
  // localHostData.buildSubList (1, 2, 1, 2, 1, 2, hsub);
  // hsub.setHostData (localHostData);
  // HostCellListedMDData hdata (localHostData);
  // localHostData.add (hdata);

  // HostTransferPackage hpkg;
  // hpkg.reinit (subList);
  // DeviceTransferPackage dpkg;
  // dpkg.reinit (subList);
  // hpkg.pack (localHostData);
  // dpkg.copyFromHost (hpkg);
  // dpkg.unpack_add (deviceData);
  // hpkg.unpack_replace (localHostData);

  // DeviceCellRelation relation;
  // relation.build (deviceData);
  // DeviceBondList dbdlist (deviceData);
  // buildDeviceBondList (deviceData, relation, dbdlist);
  // HostBondList hdblist;
  // dbdlist.copyToHost (hdblist);
  
  // HostCellListedMDData hdata1;
  // hdata1.copy (localHostData,
  // 	       MDDataItemMask_All ^
  // 	       MDDataItemMask_Bond ^
  // 	       MDDataItemMask_Angle);
  // hpkg.unpack_replace (hdata1);
  // DeviceCellListedMDData ddata1;
  // ddata1.copyFromDevice (deviceData,
  // 			 MDDataItemMask_All ^
  // 			 MDDataItemMask_Bond ^
  // 			 MDDataItemMask_Angle);
  // dpkg.unpack_replace (ddata1);
  
  // hsub.add (hsub);
  // DeviceTransferPackage dpkg;
  // dpkg.reinit (subList);
  // dpkg.pack (deviceData);
  // HostTransferPackage hpkg;
  // hpkg.reinit (subList);
  // dpkg.copyToHost (hpkg);
  // hpkg.unpack_replace (localHostData);

  hostBuff.copy (localHostData, MDDataItemMask_All);
  hostBuff.clearData ();
  
  redistribtransUtil  .setHostData (localHostData, hostBuff);
  redistribcopyUtil   .setData     (localHostData, deviceData);
  transCoordstransUtil.setHostData (localHostData, hostBuff);
  transCoordscopyUtil .setData     (localHostData, deviceData);
  
  collectUtil.setHostData (localHostData, hostBuff, globalHostData);
  
  return;
}


bool Parallel::MDSystem::
reinitCellStructure (const ScalorType & cellSize,
		     const IndexType & divideLevel)
{
  bool reinited = deviceData.reinitCellStructure (cellSize, divideLevel);
  if (reinited){
    cellRelation.rebuild (deviceData);

    deviceData.mallocToHost   (localHostData);
    hostBuff  .mallocFromHost (localHostData);
  
    redistribtransUtil  .setHostData (localHostData, hostBuff);
    redistribcopyUtil   .setData     (localHostData, deviceData);
    transCoordstransUtil.setHostData (localHostData, hostBuff);
    transCoordscopyUtil .setData     (localHostData, deviceData);
  
    collectUtil.setHostData (localHostData, hostBuff, globalHostData);
  }
    
  return reinited;
}



void Parallel::MDSystem::
redistribute ()
{
  DeviceTimer::tic (item_ApplyBondaryCondition);
  deviceData.applyPeriodicBondaryConditionOnGhostCells ();
  DeviceTimer::toc (item_ApplyBondaryCondition);
  HostTimer::tic(item_Redistribute_DHCopy);
  localHostData.clearData ();
  redistribcopyUtil .copyToHost ();
  redistribcopyUtil .clearDeviceSent ();
  HostTimer::toc(item_Redistribute_DHCopy);
  redistribtransUtil.redistributeHost ();
  HostTimer::tic(item_Redistribute_DHCopy);
  redistribcopyUtil .copyFromHost ();
  HostTimer::toc(item_Redistribute_DHCopy);
}

void Parallel::MDSystem::
transferGhost ()
{
  HostTimer::tic (item_TransferGhost_DHCopy);
  transCoordscopyUtil .copyToHost ();
  HostTimer::toc (item_TransferGhost_DHCopy);
  HostTimer::tic (item_TransferGhost_Tranfer);
  transCoordstransUtil.transCoords ();
  HostTimer::toc (item_TransferGhost_Tranfer);
  HostTimer::tic (item_TransferGhost_DHCopy);
  transCoordscopyUtil .copyFromHost ();
  HostTimer::toc (item_TransferGhost_DHCopy);
}

void Parallel::MDSystem::
clearGhost ()
{
  transCoordscopyUtil.clearGhost();
}

void Parallel::MDSystem::
collectLocalData ()
{
  collectUtil.collect();
}

void Parallel::MDSystem::
writeGlobalData_GroFile (const char * filename)
{
  if (Parallel::Interface::myRank() == 0){
    globalHostData.writeData_GroFile (filename,
				      atomName, atomIndex,
				      resdName, resdIndex);
  }
}

void Parallel::MDSystem::
finalize ()
{
  redistribtransUtil.clear ();
  transCoordstransUtil.clear ();
}

Parallel::SystemRedistributeCopyUtil::
SystemRedistributeCopyUtil ()
    : hpkgInner (Parallel::HostAllocator::hostMallocPageLocked),
      hpkgOuter (Parallel::HostAllocator::hostMallocPageLocked)
{
  mask = MDDataItemMask_All;
}

void Parallel::SystemRedistributeCopyUtil::
setData (HostCellListedMDData & hdata,
	 DeviceCellListedMDData & ddata)
{
  SubCellList tmp;
  ptr_hdata = &hdata;
  IndexType devideLevel = ptr_hdata->getDevideLevel ();
  HostIntVectorType numCell = ptr_hdata->getNumCell ();

  ptr_hdata->buildSubListGhostCell (hostSubOuter);
  ptr_hdata->buildSubListRealCell (hostSubInner);
  ptr_hdata->buildSubList (2 * devideLevel, numCell.x - 2 * devideLevel,
			   2 * devideLevel, numCell.y - 2 * devideLevel,
			   2 * devideLevel, numCell.z - 2 * devideLevel,
			   tmp);
  hostSubInner.sub (tmp);

  ptr_ddata = &ddata;
  ptr_ddata->buildSubListGhostCell (deviceSubOuter);
  ptr_ddata->buildSubListRealCell (deviceSubInner);
  ptr_ddata->buildSubList (2 * devideLevel, numCell.x - 2 * devideLevel,
			   2 * devideLevel, numCell.y - 2 * devideLevel,
			   2 * devideLevel, numCell.z - 2 * devideLevel,
			   tmp);
  deviceSubInner.sub (tmp);

  hpkgInner.reinit (hostSubInner);
  hpkgOuter.reinit (hostSubOuter);
  dpkgInner.reinit (deviceSubInner);
  dpkgOuter.reinit (deviceSubOuter);
}


Parallel::SystemTransCoordsCopyUtil::
SystemTransCoordsCopyUtil ()
    : mask (MDDataItemMask_Coordinate |
	    MDDataItemMask_GlobalIndex),
      hpkgInner (Parallel::HostAllocator::hostMallocPageLocked),
      hpkgOuter (Parallel::HostAllocator::hostMallocPageLocked)
{
}

void Parallel::SystemTransCoordsCopyUtil::
setData (HostCellListedMDData & hdata,
	 DeviceCellListedMDData & ddata)
{
  SubCellList tmp;
  ptr_hdata = &hdata;
  IndexType devideLevel = ptr_hdata->getDevideLevel ();
  HostIntVectorType numCell = ptr_hdata->getNumCell ();

  ptr_hdata->buildSubListGhostCell (hostSubOuter);
  ptr_hdata->buildSubListRealCell (hostSubInner);
  ptr_hdata->buildSubList (2 * devideLevel, numCell.x - 2 * devideLevel,
			   2 * devideLevel, numCell.y - 2 * devideLevel,
			   2 * devideLevel, numCell.z - 2 * devideLevel,
			   tmp);
  hostSubInner.sub (tmp);

  ptr_ddata = &ddata;
  ptr_ddata->buildSubListGhostCell (deviceSubOuter);
  ptr_ddata->buildSubListRealCell (deviceSubInner);
  ptr_ddata->buildSubList (2 * devideLevel, numCell.x - 2 * devideLevel,
			   2 * devideLevel, numCell.y - 2 * devideLevel,
			   2 * devideLevel, numCell.z - 2 * devideLevel,
			   tmp);
  deviceSubInner.sub (tmp);

  hpkgInner.reinit (hostSubInner);
  hpkgOuter.reinit (hostSubOuter);
  dpkgInner.reinit (deviceSubInner);
  dpkgOuter.reinit (deviceSubOuter);
}




