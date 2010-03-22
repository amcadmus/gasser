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
			     & (globalHostData.cptr_coordinateNoiX()[i]),
			     & (globalHostData.cptr_coordinateNoiY()[i]),
			     & (globalHostData.cptr_coordinateNoiZ()[i]));
    }
    globalHostData.initTopology (sysTop);
    for (int i = 0; i < globalNumAtom; ++i){
      globalHostData.cptr_coordinate()[i].w = globalHostData.cptr_type()[i];
    }
  }
  
  distributeGlobalMDData (globalHostData, localHostData);

  // deviceData.mallocAll (localHostData.memSize());
  DeviceMDData & ddata (deviceData);
  ddata.copyFromHost (localHostData, MDDataItemMask_All);
  // deviceData.copyToHost   (localHostData);  

  deviceData.initCellStructure (2);

  // deviceData.rebuild ();

  // deviceData.dptr_coordinate()[124].x = 7;
  // deviceData.dptr_coordinate()[124].y = 7;
  // deviceData.dptr_coordinate()[124].z = 7;

  // deviceData.coord[124].x = -0.5;
  // deviceData.coord[128].x = -0.5;
  // deviceData.coord[132].x = -0.5;
  // deviceData.dptr_mass()[124] = 7;
  // deviceData.dptr_mass()[128] = 7;
  // deviceData.dptr_mass()[132] = 7;
  
  for (IndexType i = 0; i < deviceData.numData(); ++i){
    deviceData.dptr_coordinate()[i].x += 2;
    deviceData.dptr_coordinate()[i].y += 2;
    deviceData.dptr_coordinate()[i].z += 2;
  }
  
  // printf ("rank %d\n", Parallel::Interface::myRank());
  deviceData.rebuild ();
  printf ("rank %d\n", Parallel::Interface::myRank());

  deviceData.copyToHost (localHostData, MDDataItemMask_All);
  hostBuff.copy (localHostData, MDDataItemMask_All);
  hostBuff.clearData ();
  
  // hostBuff.copy (localHostData, MDDataItemMask_All);  
  // localHostData.copy (hostBuff);
  
  // HostSubCellList sub0, sub1;
  // localHostData.buildSubListRealCell (sub0);
  // hostBuff.buildSubListRealCell (sub1);
  // sub1.add (sub0, MDDataItemMask_AllExceptForce);
  
  // // DeviceCellListedMDData ddata1;
  // // ddata1.copyFromHost (localHostData);
  // // DeviceCellListedMDData ddata2;
  // // ddata2.copyFromDevice (deviceData);

  // SubCellList subList;
  // deviceData.buildSubListAllCell (subList);
  // //  for (IndexType i = 0; i < subList.size(); ++i){
  // //   IndexType ix, iy, iz;
  // //   deviceData.D1toD3 (subList[i], ix, iy, iz);
  // //   printf ("%d %d %d\n", ix, iy, iz);
  // // }
 
  // deviceData.buildSubListRealCell (subList);
  // // for (IndexType i = 0; i < subList.size(); ++i){
  // //   IndexType ix, iy, iz;
  // //   deviceData.D1toD3 (subList[i], ix, iy, iz);
  // //   printf ("%d %d %d\n", ix, iy, iz);
  // // }

  // deviceData.buildSubListGhostCell (subList);

  // Parallel::DeviceTransferPackage dpkg ;
  // MDDataItemMask_t mask = MDDataItemMask_AllExceptForce;
  // dpkg.reinit (subList);
  // dpkg.pack (deviceData, mask);

  // Parallel::HostTransferPackage hpkg;
  // dpkg.copyToHost (hpkg);

  // hpkg.cptr_forceX()[0] = 1.2;
  // hpkg.cptr_mass()[0] = 7;
  // hpkg.cptr_mass()[1] = 7;
  // hpkg.cptr_mass()[2] = 7;

  // dpkg.copyFromHost (hpkg);
  // dpkg.unpack_add (deviceData);
  // hpkg.unpack_add (localHostData);
  
  // int i = 1;

  redistribtransUtil.setHostData (localHostData, hostBuff);
  redistribcopyUtil .setData     (localHostData, deviceData);
  
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
  mask = MDDataItemMask_AllExceptForce;
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
  ptr_hdata->buildSubListAllCell (hostSubInner);
  ptr_hdata->buildSubList (2 * devideLevel, numCell.x - 2 * devideLevel,
			   2 * devideLevel, numCell.y - 2 * devideLevel,
			   2 * devideLevel, numCell.z - 2 * devideLevel,
			   tmp);
  hostSubInner.sub (tmp);

  ptr_ddata = &ddata;
  ptr_ddata->buildSubListGhostCell (deviceSubOuter);
  ptr_ddata->buildSubListAllCell (deviceSubInner);
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

