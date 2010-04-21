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
#include "SystemNonBondedInteraction.h"
#include "Parallel_Timer.h"

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
  genv.setDeviceId (atoi(argv[3]));
  
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
  Topology::Molecule mol;
  mol.pushAtom (Topology::Atom (1.0, 0.0, 0));
  LennardJones6_12Parameter ljparam;
  ljparam.reinit (1.f, 1.f, 0.f, 3.2f);
  sysTop.addNonBondedInteraction (Topology::NonBondedInteraction(0, 0, ljparam));
  sysTop.addMolecules (mol, sys.numAtomInGroFile(filename));
  
  sys.init (filename, sysTop);
  sys.redistribute ();
  sys.deviceData.applyPeriodicBondaryCondition ();
  SystemNonBondedInteraction sysNbInter (sysTop);

  Parallel::InteractionEngine interEng (sys.deviceData);
  interEng.registNonBondedInteraction (sysNbInter);
  interEng.clearInteraction (sys.deviceData);
  Parallel::DeviceCellRelation relation;
  relation.build (sys.deviceData);

  Parallel::HostStatistic hst;
  Parallel::DeviceStatistic dst;
  hst.reinit (sys.localHostData);
  dst.reinit (sys.deviceData);
  
  Parallel::Integrator::VelocityVerlet vv (sys.deviceData);
  Parallel::TranslationalFreedomRemover trRemover (sys.deviceData);
  ScalorType dt = 0.001;

  sys.globalHostData.initWriteData_xtcFile ("traj.xtc");

  IndexType stFeq = 10;
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
    if ((i+1) % stFeq == 0){
      DeviceTimer::tic (item_NonBondedInterStatistic);
      interEng.applyNonBondedInteraction (sys.deviceData, relation, dst);
      DeviceTimer::toc (item_NonBondedInterStatistic);
    }
    else {
      DeviceTimer::tic (item_NonBondedInteraction);
      interEng.applyNonBondedInteraction (sys.deviceData, relation);
      DeviceTimer::toc (item_NonBondedInteraction);
    }
    HostTimer::tic (item_TransferGhost);
    sys.clearGhost ();
    HostTimer::toc (item_TransferGhost);
    DeviceTimer::tic (item_Integrate);
    vv.step2 (sys.deviceData, dt, dst);
    DeviceTimer::toc (item_Integrate);
    DeviceTimer::tic (item_BuildCellList);
    sys.deviceData.rebuild ();
    DeviceTimer::toc (item_BuildCellList);
    DeviceTimer::tic (item_ApplyBondaryCondition);
    sys.deviceData.applyPeriodicBondaryCondition ();
    DeviceTimer::toc (item_ApplyBondaryCondition);
    HostTimer::tic (item_Redistribute);
    sys.redistribute ();
    HostTimer::toc (item_Redistribute);
  
    dst.copyToHost (hst);
    hst.collectData ();

    if ((i+1)%stFeq == 0){
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
      sys.updateHost ();  
      sys.collectLocalData ();
      sys.globalHostData.writeData_xtcFile (i, dt*i);
    }
  }
  
  // int myRank = Parallel::Interface::myRank();
  // printf ("rank: %d, %f\n", myRank, hst.NonBondedEnergy());

  
  // sys.transferCoordinate ();  
  
  // for (IndexType i = 0; i < sys.deviceData.numData(); ++i){
  //   sys.deviceData.dptr_coordinate()[i].x += 2;
  //   sys.deviceData.dptr_coordinate()[i].y += 2;
  //   sys.deviceData.dptr_coordinate()[i].z += 2;
  // }
  // sys.deviceData.rebuild ();
  // sys.redistribute ();
  // sys.deviceData.applyPeriodicBondaryCondition ();
  
  // for (IndexType i = 0; i < sys.deviceData.numData(); ++i){
  //   sys.deviceData.dptr_coordinate()[i].x += 2;
  //   sys.deviceData.dptr_coordinate()[i].y += 2;
  //   sys.deviceData.dptr_coordinate()[i].z += 2;
  // }
  // sys.deviceData.rebuild ();  
  // sys.redistribute ();
  // sys.deviceData.applyPeriodicBondaryCondition ();


  // for (IndexType k = 0; k < Parallel::Interface::numProc(); ++k){
  //   if (Parallel::Interface::myRank() == k){
  //     IndexType totalNumCell = sys.localHostData.getNumCell().x *
  // 	  sys.localHostData.getNumCell().y *
  // 	  sys.localHostData.getNumCell().z ;
  //     for (IndexType i = 0; i < totalNumCell; ++i){
  // 	IndexType ix, iy, iz;
  // 	IndexType numThreadsInCell = Parallel::Interface::numThreadsInCell();
  // 	IndexType cellShift = numThreadsInCell * i ;
  // 	sys.localHostData.D1toD3 (i, ix, iy, iz);
  // 	for (IndexType j = 0; j < sys.localHostData.cptr_numAtomInCell()[i]; ++j){
  // 	  printf ("%d %d %d   %f %f %f\n",
  // 		  ix, iy, iz,
  // 		  sys.localHostData.cptr_coordinate()[j+cellShift].x,
  // 		  sys.localHostData.cptr_coordinate()[j+cellShift].y,
  // 		  sys.localHostData.cptr_coordinate()[j+cellShift].z);
  // 	}
  //     }
  //   }
  //   Parallel::Interface::barrier();
  // }

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

