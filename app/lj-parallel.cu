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
#include "Parallel_HostAllocator.h"

#include "compile_error_mixcode.h"

using namespace Parallel::Timer;

int prog (int argc, char * argv[]);

int main(int argc, char * argv[])
{
  Parallel::Interface::initMPI (&argc, &argv);
  Parallel::Interface::initEnvironment (20, "Device Emulation (CPU)");
  // Parallel::Interface::initEnvironment (64, "Tesla C1060");
  
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

int prog (int argc, char * argv[])
{  
  HostTimer::init();
  DeviceTimer::init ();
  HostTimer::tic (item_Total);
  
  IndexType nstep = 20;
  char * filename;  
  if (argc != 3){
    printf ("Usage:\n%s conf.gro nstep\n", argv[0]);
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

  SystemNonBondedInteraction sysNbInter (sysTop);
  SystemBondedInteraction    sysBdInter (sysTop);

  sys.init (filename, sysTop, sysNbInter.maxRcut(), 1);
  sys.redistribute ();

  Parallel::InteractionEngine interEng (sys.deviceData);
  interEng.registNonBondedInteraction (sysNbInter);
  interEng.registBondedInteraction    (sysBdInter);
  interEng.clearInteraction (sys.deviceData);
  Parallel::DeviceCellRelation relation;
  relation.rebuild (sys.deviceData);
  Parallel::HostStatistic hst;
  Parallel::DeviceStatistic dst;
  hst.reinit (sys.localHostData);
  dst.reinit (sys.deviceData);
  
  Parallel::Integrator::VelocityVerlet vv (sys.deviceData);
  Parallel::TranslationalFreedomRemover trRemover (sys.deviceData);
  ScalorType dt = 0.001;

  sys.globalHostData.initWriteData_xtcFile ("traj.xtc");

  IndexType stFeq = 100;
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

