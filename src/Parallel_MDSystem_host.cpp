#define HOST_CODE
#define CPP_FILE

#include "Parallel_MDSystem.h"
#include "Parallel_TransferEngine.h"

Parallel::MDSystem::
MDSystem (const Environment & env_)
    : atomName(NULL), atomIndex(NULL),
      resdName(NULL), resdIndex(NULL)
{
  env = &env_;
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
  Parallel::TransferEngine transSend (*env);
  Parallel::TransferEngine transRecv (*env);
  
  if (env->myRank() == 0){
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
  }
  
  int Nx, Ny, Nz;
  env->numProcDim (Nx, Ny, Nz);
  ScalorType gbx, gby, gbz;
  ScalorType hx(0), hy(0), hz(0);
  
  if (env->myRank() == 0){
    IndexType naiveLocalMemSize  = (globalNumAtom / env->numProc() * 2);
    if (naiveLocalMemSize < 100)  naiveLocalMemSize = 100 ;
    sendHostData.reallocAll (naiveLocalMemSize);
    
    hx = (gbx = globalHostData.getGlobalBox().size.x) / Nx;
    hy = (gby = globalHostData.getGlobalBox().size.y) / Ny;
    hz = (gbz = globalHostData.getGlobalBox().size.z) / Nz;
  }
  
  for (int ix = 0; ix < Nx; ++ix){
    for (int iy = 0; iy < Ny; ++iy){
      for (int iz = 0; iz < Nz; ++iz){

	if (env->myRank() == 0){
	  sendHostData.clearData();
	  ScalorType lx = ix * hx;
	  ScalorType ux = lx + hx;
	  ScalorType ly = iy * hy;
	  ScalorType uy = ly + hy;
	  ScalorType lz = iz * hz;
	  ScalorType uz = lz + hz;
	  for (IndexType i = 0; i < globalHostData.numAtom(); ++i){
	    if (globalHostData.cptr_coordinate()[i].x >= lx &&
		globalHostData.cptr_coordinate()[i].x <  ux &&
		globalHostData.cptr_coordinate()[i].y >= ly &&
		globalHostData.cptr_coordinate()[i].y <  uy &&
		globalHostData.cptr_coordinate()[i].z >= lz &&
		globalHostData.cptr_coordinate()[i].z <  uz){
	      sendHostData.pushBackAtom (
		  globalHostData.cptr_coordinate()[i],
		  globalHostData.cptr_coordinateNoiX()[i],
		  globalHostData.cptr_coordinateNoiY()[i],
		  globalHostData.cptr_coordinateNoiZ()[i],
		  globalHostData.cptr_velocityX()[i],
		  globalHostData.cptr_velocityY()[i],
		  globalHostData.cptr_velocityZ()[i],
		  globalHostData.cptr_globalIndex()[i],
		  globalHostData.cptr_type()[i],
		  globalHostData.cptr_mass()[i],
		  globalHostData.cptr_charge()[i]);
	    }
	  }
	}

	int workingRank;
	env->cartCoordToRank (ix, iy, iz, workingRank);
	int recvNumAtom, recvMemSize;
	int sendNumAtom, sendMemSize;
	
	if (env->myRank() == 0){
	  sendNumAtom = sendHostData.numAtom();
	  sendMemSize = sendHostData.memSize();
	  transSend.clearRegistered ();
	  transSend.registerBuff (&sendNumAtom, sizeof(IndexType));
	  transSend.registerBuff (&sendMemSize, sizeof(IndexType));
	  transSend.build ();
	  transSend.Isend (workingRank, 0);
	  // transSend.wait ();
	}
	if (env->myRank() == workingRank){
	  transRecv.clearRegistered();
	  transRecv.registerBuff (&recvNumAtom, sizeof(IndexType));
	  transRecv.registerBuff (&recvMemSize, sizeof(IndexType));
	  transRecv.build ();
	  transRecv.Irecv (0, 0);
	  transRecv.wait ();
	  localHostData.reallocAll (recvMemSize);
	}
	if (env->myRank () == 0){
	  transSend.wait ();
	}
	if (env->myRank () == 0){
	  transSend.clearRegistered ();
	  transSend.registerBuff (sendHostData.cptr_coordinate(),
				  sizeof(HostCoordType) * sendHostData.numAtom());
	  transSend.registerBuff (sendHostData.cptr_coordinateNoiX(),
				  sizeof(IntScalorType) * sendHostData.numAtom());
	  transSend.registerBuff (sendHostData.cptr_coordinateNoiY(),
				  sizeof(IntScalorType) * sendHostData.numAtom());
	  transSend.registerBuff (sendHostData.cptr_coordinateNoiZ(),
				  sizeof(IntScalorType) * sendHostData.numAtom());
	  transSend.registerBuff (sendHostData.cptr_velocityX(),
				  sizeof(ScalorType) * sendHostData.numAtom());
	  transSend.registerBuff (sendHostData.cptr_velocityY(),
				  sizeof(ScalorType) * sendHostData.numAtom());
	  transSend.registerBuff (sendHostData.cptr_velocityZ(),
				  sizeof(ScalorType) * sendHostData.numAtom());
	  transSend.registerBuff (sendHostData.cptr_globalIndex(),
				  sizeof(IndexType) * sendHostData.numAtom());
	  transSend.registerBuff (sendHostData.cptr_type(),
				  sizeof(TypeType) * sendHostData.numAtom());
	  transSend.registerBuff (sendHostData.cptr_mass(),
				  sizeof(ScalorType) * sendHostData.numAtom());
	  transSend.registerBuff (sendHostData.cptr_charge(),
				  sizeof(ScalorType) * sendHostData.numAtom());
	  ScalorType tmpx, tmpy, tmpz;
	  tmpx = sendHostData.getGlobalBox().size.x;
	  tmpy = sendHostData.getGlobalBox().size.y;
	  tmpz = sendHostData.getGlobalBox().size.z;
	  transSend.registerBuff (&tmpx,
				  sizeof(ScalorType) * 1);
	  transSend.registerBuff (&tmpy,
				  sizeof(ScalorType) * 1);
	  transSend.registerBuff (&tmpz,
				  sizeof(ScalorType) * 1);
	  transSend.build ();
	  transSend.Isend (workingRank, 0);
	}
	if (env->myRank () == workingRank){
	  transRecv.clearRegistered();
	  transRecv.registerBuff (localHostData.cptr_coordinate(),
				  sizeof(HostCoordType) * recvNumAtom);
	  transRecv.registerBuff (localHostData.cptr_coordinateNoiX(),
				  sizeof(IntScalorType) * recvNumAtom);
	  transRecv.registerBuff (localHostData.cptr_coordinateNoiY(),
				  sizeof(IntScalorType) * recvNumAtom);
	  transRecv.registerBuff (localHostData.cptr_coordinateNoiZ(),
				  sizeof(IntScalorType) * recvNumAtom);
	  transRecv.registerBuff (localHostData.cptr_velocityX(),
				  sizeof(ScalorType) * recvNumAtom);
	  transRecv.registerBuff (localHostData.cptr_velocityY(),
				  sizeof(ScalorType) * recvNumAtom);
	  transRecv.registerBuff (localHostData.cptr_velocityZ(),
				  sizeof(ScalorType) * recvNumAtom);
	  transRecv.registerBuff (localHostData.cptr_globalIndex(),
				  sizeof(IndexType) * recvNumAtom);
	  transRecv.registerBuff (localHostData.cptr_type(),
				  sizeof(TypeType) * recvNumAtom);
	  transRecv.registerBuff (localHostData.cptr_mass(),
				  sizeof(ScalorType) * recvNumAtom);
	  transRecv.registerBuff (localHostData.cptr_charge(),
				  sizeof(ScalorType) * recvNumAtom);
	  ScalorType tmpx, tmpy, tmpz;
	  transRecv.registerBuff (&tmpx,
				  sizeof(ScalorType) * 1);
	  transRecv.registerBuff (&tmpy,
				  sizeof(ScalorType) * 1);
	  transRecv.registerBuff (&tmpz,
				  sizeof(ScalorType) * 1);
	  transRecv.build ();
	  transRecv.Irecv (0, 0);
	  transRecv.wait ();
	  localHostData.setGlobalBox (tmpx, tmpy, tmpz);
	  localHostData.numAtom () = recvNumAtom;
	}
	if (env->myRank () == 0){
	  transSend.wait ();
	}
      }
    }
  }


  // MPI_Barrier (env->communicator());
  // SummationEngine sume (*env);
  // IndexType *sumNum;
  // sume.sumIndexAll (&localHostData.numAtom(), 1, &sumNum);
  // printf ("myrank is %d, my num atom is %d\n", env->myRank(), localHostData.numAtom());
  // fflush (stdout);
  // if (env->myRank() == 0){
  //   printf ("sumTotal is %d, global is %d\n", sumNum[0], globalHostData.numAtom());
  //   fflush (stdout);
  // }
  
}
  



