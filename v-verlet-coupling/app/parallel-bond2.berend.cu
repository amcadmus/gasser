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
  Parallel::Interface::initEnvironment (64, "Tesla C1060");
  
  int flag = 0;
  if (Parallel::Interface::isActive()){
    try{
      flag = prog (argc, argv);
    }
    catch (MDExcptCuda & e){
      fprintf (stderr, "%s\n", e.what());	
    }
    catch (MDException &e){
      fprintf (stderr, "%s\n", e.what());
    }
  }

  Parallel::Interface::finalizeEnvironment ();
  
  return flag;
}


int prog (int argc, char *argv[])
{
    
  HostTimer::init();
  DeviceTimer::init ();
  HostTimer::tic (item_Total);
  
  IndexType nstep = 20;
  char * filename;  
  if (argc != 3){
    printf ("Usage:\n%s conf.gro nstep \n", argv[0]);
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

  SystemNonBondedInteraction sysNbInter (sysTop);
  SystemBondedInteraction    sysBdInter (sysTop);
  // if (Parallel::Interface::myRank() == 0) sysBdInter.printEverything ();
  
  sys.init (filename, sysTop, sysNbInter.maxRcut() + .1f);
  // sys.init (filename, sysTop, 3.34);
  sys.redistribute ();
  sys.globalHostData.initWriteData_xtcFile ("traj.xtc");

  Parallel::InteractionEngine interEng;
  interEng.reinit(sys.deviceData);
  interEng.registNonBondedInteraction (sysNbInter);
  interEng.registBondedInteraction    (sysBdInter);
  interEng.clearInteraction (sys.deviceData);
  Parallel::DeviceCellRelation relation;
  relation.rebuild (sys.deviceData);
  Parallel::DeviceBondList dbdlist;
  dbdlist.reinit (sys.deviceData);
  buildDeviceBondList (sys.deviceData, relation, dbdlist);
  Parallel::HostStatistic hst;
  Parallel::DeviceStatistic dst;
  
  ScalorType dt = 0.001;
  Parallel::TranslationalFreedomRemover trRemover;
  trRemover.reinit (sys.deviceData);
  Parallel::Integrator::BerendsenLeapFrog blf;
  blf.reinit (sys, dt, &interEng, &trRemover, 10, &dbdlist);

  blf.TCouple (1.f, 1);
  blf.addPcoupleGroup (mdRectBoxDirectionX |
  		       mdRectBoxDirectionY |
  		       mdRectBoxDirectionZ,
  		       1, 1, 1);
  // blf.addPcoupleGroup (mdRectBoxDirectionX |
  // 		       mdRectBoxDirectionY,
  // 		       0, 1, 10);
  // blf.addPcoupleGroup (mdRectBoxDirectionZ,
  // 		       0, 1, 1);

  if (Parallel::Interface::myRank() == 0){
    printf ("#       1            2              3             4             5             6             7           8           9          10          11\n");
    printf ("#       n         time           nb E           b E           k E             T       total E         Pxx         Pyy         Pzz           P\n");
  }
  
  IndexType stFeq = 10;
  for (IndexType i = 0; i < nstep; ++i){
    if ((i+1)%stFeq == 0){
      dst.clearData ();
      blf.oneStep (sys, dst);
      dst.copyToHost (hst);
      hst.collectData ();
    }
    else {
      blf.oneStep (sys);
    }    
    if ((i+1)%stFeq == 0 && Parallel::Interface::myRank() == 0){
      printf ("%09d %07e %.7e %.7e %.7e %.7e %.7e %.5e %.5e %.5e %.5e %.3f %.3f %.3f\n",
  	      (i+1),  
  	      (i+1) * dt, 
  	      hst.NonBondedEnergy(),
  	      hst.BondedEnergy(),
  	      hst.kineticEnergy(),
  	      hst.kineticEnergy() * 2. / sys.getNumFreedom(),
  	      hst.NonBondedEnergy() +
  	      hst.BondedEnergy() +
  	      hst.kineticEnergy(),
  	      hst.pressureXX(sys.deviceData.getGlobalBox()),
  	      hst.pressureYY(sys.deviceData.getGlobalBox()),
  	      hst.pressureZZ(sys.deviceData.getGlobalBox()),
  	      hst.pressure(sys.deviceData.getGlobalBox()),
	      sys.deviceData.getGlobalBoxSize().x,
	      sys.deviceData.getGlobalBoxSize().y,
	      sys.deviceData.getGlobalBoxSize().z
	  );
      fflush (stdout);
    }
    
    if ((i+1)%1000 == 0){
      DeviceTimer::tic (item_DataIO);
      sys.updateHost ();  
      sys.collectLocalData ();
      sys.globalHostData.writeData_xtcFile (i, dt*i);
      DeviceTimer::toc (item_DataIO);
    }
  }

  sys.updateHost ();
  sys.collectLocalData ();
  sys.writeGlobalData_GroFile ("confout.gro");
  sys.globalHostData.endWriteData_xtcFile ();
  
  // char name[1024];
  // sprintf (name, "id%d.coord", Parallel::Interface::myRank());
  // sys.writeLocalData_SimpleFile (name);

  Parallel::Interface::barrier();
  sys.finalize();
  
  HostTimer::toc (item_Total);
  
  printRecord (stderr);

  DeviceTimer::finalize ();
  HostTimer::finalize ();

  return 0;
}
