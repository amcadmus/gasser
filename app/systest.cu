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

int prog (int argc, char * argv[]);

int main(int argc, char * argv[])
{
  Parallel::Interface::initMPI (&argc, &argv);
  // Parallel::Interface::initEnvironment ("Device Emulation (CPU)");
  Parallel::Interface::initEnvironment ();

  int flag = 0;
  if (Parallel::Interface::isActive()){
    flag = prog (argc, argv);
  }
  
  Parallel::Interface::finalizeEnvironment ();

  return flag;
}

int prog (int argc, char * argv[])
{  
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
  LennardJones6_12CapParameter ljparamCap;
  ljparamCap.reinit (1.f, 1.f, 0.f, 3.2f, 10000.f);
  LennardJones6_12Parameter ljparam;
  ljparam.reinit (1.f, 1.f, 0.f, 3.2f);
  sysTop.addNonBondedInteraction (Topology::NonBondedInteraction(0, 0, ljparam));
  sysTop.addMolecules (mol, sys.numAtomInGroFile(filename));

  // Topology::System sysTop;
  // LennardJones6_12Parameter ljparam;
  // ljparam.reinit (1.f, 1.f, 0.f, 3.2f);
  // sysTop.addNonBondedInteraction (Topology::NonBondedInteraction(0, 0, ljparam));
  // Topology::Molecule mol;
  // mol.pushAtom (Topology::Atom (1.0, 0.0, 0));
  // mol.pushAtom (Topology::Atom (1.0, 0.0, 0));
  // HarmonicSpringParameter hsparam;
  // hsparam.reinit (10.f, 1.f);
  // mol.addBond (Topology::Bond (0, 1, hsparam));
  // mol.addBond (Topology::Bond (0, 1, hsparam));
  // Topology::Molecule mol1;
  // mol.pushAtom (Topology::Atom (1.0, 0.0, 0));  
  // sysTop.addMolecules (mol, 10);
  // sysTop.addMolecules (mol1, sys.numAtomInGroFile(filename) - 20);

  // // conf.bond2.gro
  // Topology::System sysTop;
  // LennardJones6_12Parameter ljparam;
  // ljparam.reinit (0.f, 1.f, 0.f, 3.2f);
  // sysTop.addNonBondedInteraction (Topology::NonBondedInteraction(0, 0, ljparam));
  // Topology::Molecule mol;
  // mol.pushAtom (Topology::Atom (1.0, 0.0, 0));
  // mol.pushAtom (Topology::Atom (1.0, 0.0, 0));
  // mol.pushAtom (Topology::Atom (1.0, 0.0, 0));
  // mol.pushAtom (Topology::Atom (1.0, 0.0, 0));
  // mol.pushAtom (Topology::Atom (1.0, 0.0, 0));
  // mol.pushAtom (Topology::Atom (1.0, 0.0, 0));
  // HarmonicSpringParameter hsparam;
  // hsparam.reinit (10.f, 1.2f);
  // mol.addBond (Topology::Bond (0, 1, hsparam));
  // hsparam.reinit (11.f, 1.f);
  // mol.addBond (Topology::Bond (0, 1, hsparam));
  // mol.addBond (Topology::Bond (1, 2, hsparam));
  // mol.addBond (Topology::Bond (2, 3, hsparam));
  // mol.addBond (Topology::Bond (3, 4, hsparam));
  // mol.addBond (Topology::Bond (4, 5, hsparam));
  // AngleHarmonicParameter angleparam;
  // angleparam.reinit (10, M_PI / 2.);
  // mol.addAngle (Topology::Angle (0, 1, 2, angleparam));
  // angleparam.reinit (11, M_PI / 2.);
  // mol.addAngle (Topology::Angle (0, 1, 2, angleparam));
  // mol.addAngle (Topology::Angle (1, 2, 3, angleparam));
  // mol.addAngle (Topology::Angle (2, 3, 4, angleparam));
  // mol.addAngle (Topology::Angle (3, 4, 5, angleparam));
  // sysTop.addMolecules (mol, 2);

  // Topology::System sysTop;
  // Topology::Molecule mol;
  // mol.pushAtom (Topology::Atom (1.0, 0.0, 0));
  // mol.pushAtom (Topology::Atom (1.0, 0.0, 0));
  // LennardJones6_12Parameter ljparam;
  // ljparam.reinit (0.f, 1.f, 0.f, 3.2f);
  // sysTop.addNonBondedInteraction (Topology::NonBondedInteraction(0, 0, ljparam));
  // HarmonicSpringParameter hsparam;
  // hsparam.reinit (10.f, .7f);
  // mol.addBond (Topology::Bond (0, 1, hsparam));
  // // mol.addBond (Topology::Bond (0, 1, hsparam));
  // sysTop.addMolecules (mol, sys.numAtomInGroFile(filename) / 2);
  
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
  // Parallel::DeviceBondList dbdlist;
  // dbdlist.reinit (sys.deviceData);
  // buildDeviceBondList (sys.deviceData, relation, dbdlist);
  Parallel::HostStatistic hst;
  Parallel::DeviceStatistic dst;
  hst.reinit (sys.localHostData);
  dst.reinit (sys.deviceData);
  
  Parallel::Integrator::VelocityVerlet vv (sys.deviceData);
  Parallel::TranslationalFreedomRemover trRemover (sys.deviceData);
  ScalorType dt = 0.0001;

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
    DeviceTimer::tic (item_ClearInteraction);
    interEng.clearInteraction (sys.deviceData);
    DeviceTimer::toc (item_ClearInteraction);
    HostTimer::tic (item_TransferGhost);
    sys.transferGhost ();
    HostTimer::toc (item_TransferGhost);

    // buildDeviceBondList (sys.deviceData, relation, dbdlist);

    if ((i+1) % stFeq == 0){
      DeviceTimer::tic (item_NonBondedInterStatistic);
      interEng.applyNonBondedInteraction (sys.deviceData, relation, dst);
      DeviceTimer::toc (item_NonBondedInterStatistic);
      // interEng.applyBondedInteraction (sys.deviceData, dbdlist, dst);
    }
    else {
      DeviceTimer::tic (item_NonBondedInteraction);
      interEng.applyNonBondedInteraction (sys.deviceData, relation);
      DeviceTimer::toc (item_NonBondedInteraction);
      // interEng.applyBondedInteraction (sys.deviceData, dbdlist, dst);
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
    sys.deviceData.rebuild ();
    DeviceTimer::toc (item_BuildCellList);
    // DeviceTimer::tic (item_ApplyBondaryCondition);
    // sys.deviceData.applyPeriodicBondaryCondition ();
    // DeviceTimer::toc (item_ApplyBondaryCondition);
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
    
    if ((i+1)%1000 == 0){
      sys.updateHost ();  
      sys.collectLocalData ();
      sys.globalHostData.writeData_xtcFile (i, dt*i);
      // printf ("%e %e %e\n",
      // 	      sys.globalHostData.cptr_coordinate()[1].x,
      // 	      sys.globalHostData.cptr_coordinate()[1].y,
      // 	      sys.globalHostData.cptr_coordinate()[1].z);
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
  
  Parallel::Interface::barrier();
  sys.finalize();
  
  HostTimer::toc (item_Total);
  
  printRecord (stderr);

  DeviceTimer::finalize ();
  HostTimer::finalize ();

  Parallel::Interface::finalizeEnvironment ();
  
  return 0;
}

