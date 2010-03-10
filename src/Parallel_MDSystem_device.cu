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
    globalHostData.reallocAll (globalNumAtom);
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
  deviceData.copyFromHost (localHostData);
  // deviceData.copyToHost   (localHostData);  

  deviceData.initCellStructure (2);

  // deviceData.rebuild ();

  deviceData.coord[124].z = 3.1;
  deviceData.coord[128].z = 0.9;
  deviceData.coord[132].z = 1.1;
  
  deviceData.rebuild ();
}



