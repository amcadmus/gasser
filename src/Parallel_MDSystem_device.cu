#define DEVICE_CODE
#include "Parallel_MDSystem.h"
#include "Parallel_Interface.h"

#include "compile_error_mixcode.h"


Parallel::MDSystem::
MDSystem ()
    : atomName(NULL), atomIndex(NULL),
      resdName(NULL), resdIndex(NULL)
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
      const Topology::System & sysTop)
{
  if (Parallel::Interface::myRank() == 0){
    globalNumAtom = globalHostData.numAtomInGroFile (confFileName);
    // globalHostData.reallocCoordinate (globalNumAtom);
    // globalHostData.reallocVelocity   (globalNumAtom);
    // globalHostData.reallocGroProperty(globalNumAtom);
    // globalHostData.reallocGlobalIndex(globalNumAtom);
    globalHostData.easyRealloc (globalNumAtom);
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
    globalHostData.initTopology (sysTop);
    for (int i = 0; i < globalNumAtom; ++i){
      globalHostData.cptr_coordinate()[i].w = globalHostData.cptr_type()[i];
    }
  }
  
  distributeGlobalMDData (globalHostData, localHostData);

  DeviceMDData & ddata (deviceData);
  ddata.copyFromHost (localHostData, MDDataItemMask_All);

  deviceData.initCellStructure (3.2, 2);
  // printf ("ncell: %d\n", deviceData.getNumCell().x);

  cellRelation.build (deviceData);
  
  // for (IndexType i = 0; i < deviceData.numData(); ++i){
  //   deviceData.dptr_coordinate()[i].x += 2;
  //   deviceData.dptr_coordinate()[i].y += 2;
  //   deviceData.dptr_coordinate()[i].z += 2;
  // }  
  // deviceData.rebuild ();
  // printf ("rank %d\n", Parallel::Interface::myRank());

  deviceData.copyToHost (localHostData, MDDataItemMask_All);
  hostBuff.copy (localHostData, MDDataItemMask_All);
  hostBuff.clearData ();

  redistribtransUtil  .setHostData (localHostData, hostBuff);
  redistribcopyUtil   .setData     (localHostData, deviceData);
  transCoordstransUtil.setHostData (localHostData, hostBuff);
  transCoordscopyUtil .setData     (localHostData, deviceData);
  
  collectUtil.setHostData (localHostData, hostBuff, globalHostData);
  
  return;
}


void Parallel::MDSystem::
redistribute ()
{
  localHostData.clearData ();
  redistribcopyUtil .copyToHost ();
  redistribcopyUtil .clearDeviceSent ();
  redistribtransUtil.redistributeHost ();
  redistribcopyUtil .copyFromHost ();
}

void Parallel::MDSystem::
transferGhost ()
{
  transCoordscopyUtil .copyToHost ();
  transCoordstransUtil.transCoords ();
  transCoordscopyUtil .copyFromHost ();
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


Parallel::SystemRedistributeCopyUtil::
SystemRedistributeCopyUtil ()
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
    : mask (MDDataItemMask_Coordinate | MDDataItemMask_GlobalIndex)
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




