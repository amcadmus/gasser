// #define HOST_CODE
// #define CPP_FILE

#define DEVICE_CODE

#include "Parallel_MDSystem.h"
#include "NonBondedInteraction.h"
#include "GPU_Environment.h"
#include "Parallel_Interface.h"
#include "Parallel_Statistic.h"
#include "Parallel_InteractionEngine.h"
#include "Parallel_Integrator.h"
#include "Parallel_BondList.h"
#include "SystemNonBondedInteraction.h"
#include "SystemBondedInteraction.h"
#include "Parallel_Timer.h"
#include "BondInteraction.h"
#include "AngleInteraction.h"

#include "compile_error_mixcode.h"

using namespace Parallel::Timer;

int main(int argc, char * argv[])
{
  Parallel::Interface::initEnvironment (&argc, &argv);
  // int div[3];
  // div[2] = env.numProc();
  // div[1] = div[0] = 1;
  // // div[0] = 2;
  // // div[1] = 2;
  // // div[2] = 2;
  // env.init (div);
  if (Parallel::Interface::numProc() == 8){
    Parallel::Interface::initCart (2,2,2);
  }
  else{
    Parallel::Interface::initCart (1,
				   1,
				   Parallel::Interface::numProc());
  }

  
  GPU::Environment genv;
  if (Parallel::Interface::numProc() == 1){
    genv.setDeviceId (atoi(argv[3]));
  }
  else if (Parallel::Interface::myRank() == 0){
    genv.setDeviceId (0);
  }
  else{
    genv.setDeviceId(2);
  }
  
  HostTimer::init();
  DeviceTimer::init ();
  HostTimer::tic (item_Total);
  
  IndexType nstep = 20;
  char * filename;  
  if (argc != 4){
    printf ("Usage:\n%s conf.gro nstep device\n", argv[0]);
    return 1;
  }
  if (argc != 1){
    nstep = atoi(argv[2]);
    filename = argv[1];
  }

  Parallel::MDSystem sys;

  Topology::System sysTop;
  LennardJones6_12Parameter ljparam;
  ljparam.reinit (1.f, 1.f, 0.f, 3.2f);
  sysTop.addNonBondedInteraction (Topology::NonBondedInteraction(0, 0, ljparam));
  Topology::Molecule mol;
  mol.pushAtom (Topology::Atom (1.0, 0.0, 0));
  mol.pushAtom (Topology::Atom (1.0, 0.0, 0));

  HarmonicSpringParameter hsparam;
  FENEParameter feneparam;
  hsparam.reinit (10.f, 1.f);
  feneparam.reinit (30.f, 2.f);
  mol.addBond (Topology::Bond (0, 1, hsparam));
  mol.addBond (Topology::Bond (0, 1, feneparam));
  sysTop.addMolecules (mol, sys.numAtomInGroFile(filename) / 2);
  
  sys.init (filename, sysTop);

  sys.redistribute ();
  sys.deviceData.applyPeriodicBondaryCondition ();
  SystemNonBondedInteraction sysNbInter (sysTop);
  SystemBondedInteraction    sysBdInter (sysTop);

  Parallel::InteractionEngine interEng (sys.deviceData);
  interEng.registNonBondedInteraction (sysNbInter);
  interEng.registBondedInteraction    (sysBdInter);
  interEng.clearInteraction (sys.deviceData);
  Parallel::DeviceCellRelation relation;
  relation.build (sys.deviceData);
  Parallel::DeviceCellRelation relation_buildBdList;
  Parallel::SubCellList ghost, innerShell;
  sys.deviceData.buildSubListGhostCell (ghost);
  sys.deviceData.buildSubListInnerShell (innerShell);
  relation_buildBdList.build (sys.deviceData, innerShell, ghost);
  Parallel::DeviceBondList dbdlist;
  dbdlist.reinit (sys.deviceData);
  buildDeviceBondList (sys.deviceData, relation, dbdlist);
  Parallel::HostStatistic hst;
  Parallel::DeviceStatistic dst;
  hst.reinit (sys.localHostData);
  dst.reinit (sys.deviceData);
  
  Parallel::Integrator::VelocityVerlet vv (sys.deviceData);
  Parallel::TranslationalFreedomRemover trRemover (sys.deviceData);
  ScalorType dt = 0.0001;

  sys.globalHostData.initWriteData_xtcFile ("traj.xtc");

  IndexType stFeq = 1;
  for (IndexType i = 0; i < nstep; ++i){
    if ((i)%10 == 0){
      DeviceTimer::tic (item_RemoveTransFreedom);
      trRemover.remove (sys.deviceData);
      DeviceTimer::toc (item_RemoveTransFreedom);
    }
    dst.clearData ();
    DeviceTimer::tic (item_Integrate);
    vv.step1 (sys.deviceData, dt);
    DeviceTimer::toc (item_Integrate);
    DeviceTimer::tic (item_NonBondedInterStatistic);
    interEng.clearInteraction (sys.deviceData);
    DeviceTimer::toc (item_NonBondedInterStatistic);
    HostTimer::tic (item_TransferGhost);
    sys.transferGhost ();
    HostTimer::toc (item_TransferGhost);
    
    DeviceTimer::tic (item_BuildBondList);
    buildDeviceBondList (sys.deviceData, relation_buildBdList, dbdlist);
    DeviceTimer::toc (item_BuildBondList);

    if ((i+1) % stFeq == 0){
      DeviceTimer::tic (item_NonBondedInterStatistic);
      interEng.applyNonBondedInteraction (sys.deviceData, relation, dst);
      DeviceTimer::toc (item_NonBondedInterStatistic);
      DeviceTimer::tic (item_BondedInterStatistic);
      interEng.applyBondedInteraction (sys.deviceData, dbdlist, dst);
      DeviceTimer::toc (item_BondedInterStatistic);
    }
    else {
      DeviceTimer::tic (item_NonBondedInteraction);
      interEng.applyNonBondedInteraction (sys.deviceData, relation);
      DeviceTimer::toc (item_NonBondedInteraction);
      DeviceTimer::tic (item_BondedInteraction);
      interEng.applyBondedInteraction (sys.deviceData, dbdlist, dst);
      DeviceTimer::toc (item_BondedInteraction);
    }
    HostTimer::tic (item_TransferGhost);
    sys.clearGhost ();
    HostTimer::toc (item_TransferGhost);
    DeviceTimer::tic (item_Integrate);
    vv.step2 (sys.deviceData, dt, dst);
    DeviceTimer::toc (item_Integrate);

    // for (IndexType j = 0; j < sys.deviceData.DeviceMDData::memSize(); ++j){
    //   sys.deviceData.dptr_coordinate()[j].z -= 1.;
    // }

    DeviceTimer::tic (item_BuildCellList);
    sys.deviceData.rebuild (dbdlist);
    DeviceTimer::toc (item_BuildCellList);
    DeviceTimer::tic (item_ApplyBondaryCondition);
    sys.deviceData.applyPeriodicBondaryCondition ();
    DeviceTimer::toc (item_ApplyBondaryCondition);
    HostTimer::tic (item_Redistribute);
    sys.redistribute ();
    HostTimer::toc (item_Redistribute);
    
    dst.copyToHost (hst);
    hst.collectData ();

    if ((i+1)%stFeq == 0 && Parallel::Interface::myRank() == 0){
      printf ("%09d %07e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e\n",
  	      (i+1),  
  	      (i+1) * dt, 
  	      hst.NonBondedEnergy(),
  	      hst.BondedEnergy(),
  	      hst.kineticEnergy(),
  	      hst.NonBondedEnergy() +
  	      hst.BondedEnergy() +
  	      hst.kineticEnergy(),
  	      hst.pressureXX(),
  	      hst.pressureYY(),
  	      hst.pressureZZ(),
  	      hst.pressure());
      fflush (stdout);
    }
    
    if ((i+1)%100 == 0){
      DeviceTimer::tic (item_DataIO);
      sys.updateHost ();  
      sys.collectLocalData ();
      sys.globalHostData.writeData_xtcFile (i, dt*i);
      DeviceTimer::toc (item_DataIO);
    }
  }
  
  // // int myRank = Parallel::Interface::myRank();
  // // printf ("rank: %d, %f\n", myRank, hst.NonBondedEnergy());

  
  // // sys.transferCoordinate ();  
  
  // // for (IndexType i = 0; i < sys.deviceData.numData(); ++i){
  // //   sys.deviceData.dptr_coordinate()[i].x += 2;
  // //   sys.deviceData.dptr_coordinate()[i].y += 2;
  // //   sys.deviceData.dptr_coordinate()[i].z += 2;
  // // }
  // // sys.deviceData.rebuild ();
  // // sys.redistribute ();
  // // sys.deviceData.applyPeriodicBondaryCondition ();
  
  // // for (IndexType i = 0; i < sys.deviceData.numData(); ++i){
  // //   sys.deviceData.dptr_coordinate()[i].x += 2;
  // //   sys.deviceData.dptr_coordinate()[i].y += 2;
  // //   sys.deviceData.dptr_coordinate()[i].z += 2;
  // // }
  // // sys.deviceData.rebuild ();  
  // // sys.redistribute ();
  // // sys.deviceData.applyPeriodicBondaryCondition ();


  // // for (IndexType k = 0; k < Parallel::Interface::numProc(); ++k){
  // //   if (Parallel::Interface::myRank() == k){
  // //     IndexType totalNumCell = sys.localHostData.getNumCell().x *
  // // 	  sys.localHostData.getNumCell().y *
  // // 	  sys.localHostData.getNumCell().z ;
  // //     for (IndexType i = 0; i < totalNumCell; ++i){
  // // 	IndexType ix, iy, iz;
  // // 	IndexType numThreadsInCell = Parallel::Interface::numThreadsInCell();
  // // 	IndexType cellShift = numThreadsInCell * i ;
  // // 	sys.localHostData.D1toD3 (i, ix, iy, iz);
  // // 	for (IndexType j = 0; j < sys.localHostData.cptr_numAtomInCell()[i]; ++j){
  // // 	  printf ("%d %d %d   %f %f %f\n",
  // // 		  ix, iy, iz,
  // // 		  sys.localHostData.cptr_coordinate()[j+cellShift].x,
  // // 		  sys.localHostData.cptr_coordinate()[j+cellShift].y,
  // // 		  sys.localHostData.cptr_coordinate()[j+cellShift].z);
  // // 	}
  // //     }
  // //   }
  // //   Parallel::Interface::barrier();
  // // }

  sys.updateHost ();
  
  sys.collectLocalData ();

  sys.writeGlobalData_GroFile ("confout.gro");
  
  char name[1024];
  sprintf (name, "id%d.coord", Parallel::Interface::myRank());
  sys.writeLocalData_SimpleFile (name);


  sys.globalHostData.endWriteData_xtcFile ();
  
  
  HostTimer::toc (item_Total);
  
  printRecord (stderr);

  DeviceTimer::finalize ();
  HostTimer::finalize ();

  Parallel::Interface::finalizeEnvironment ();
  
  return 0;
}

