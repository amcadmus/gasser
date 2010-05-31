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
  Parallel::Interface::initEnvironment (4, "Device Emulation (CPU)");
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
  ljparam.reinit (1.f, 1.f, 0.f, 3.19f);
  sysTop.addNonBondedInteraction (Topology::NonBondedInteraction(0, 0, ljparam));
  sysTop.addMolecules (mol, sys.numAtomInGroFile(filename));

  SystemNonBondedInteraction sysNbInter (sysTop);
  SystemBondedInteraction    sysBdInter (sysTop);

  // sys.init (filename, sysTop, sysNbInter.maxRcut(), 1);
  sys.init (filename, sysTop, 4., 1);
  sys.redistribute ();
  IndexType count = 0;
  IntVectorType numCell = sys.deviceData.getNumCell();
  for (IndexType i = 0;
       i < numCell.x * numCell.y * numCell.z; ++i){
    count += sys.deviceData.dptr_numAtomInCell()[i];
  }
  printf ("num atom in sys is %d\n", count);

  sys.globalHostData.initWriteData_xtcFile ("traj.xtc");
  
  sys.reinitCellStructure (3.21);
  count = 0;
  numCell = sys.deviceData.getNumCell();
  for (IndexType i = 0;
       i < numCell.x * numCell.y * numCell.z; ++i){
    count += sys.deviceData.dptr_numAtomInCell()[i];
  }
  sys.redistribute();
  printf ("num atom in sys is %d\n", count);

  // sys.reinitCellStructure (4);
  // count = 0;
  // numCell = sys.deviceData.getNumCell();
  // for (IndexType i = 0;
  //      i < numCell.x * numCell.y * numCell.z; ++i){
  //   count += sys.deviceData.dptr_numAtomInCell()[i];
  // }
  // sys.redistribute();
  // printf ("num atom in sys is %d\n", count);

  // return 0;
  
  Parallel::InteractionEngine interEng;
  interEng.reinit (sys.deviceData);
  interEng.registNonBondedInteraction (sysNbInter);
  interEng.registBondedInteraction    (sysBdInter);
  interEng.clearInteraction (sys.deviceData);
  Parallel::DeviceCellRelation relation;
  relation.rebuild (sys.deviceData);
  Parallel::HostStatistic hst;
  Parallel::DeviceStatistic dst;
  // hst.reinit (sys.localHostData);
  // dst.reinit (sys.deviceData);
  
  ScalorType dt = 0.001;
  Parallel::Integrator::VelocityVerlet vv;
  vv.reinit(sys.deviceData);
  Parallel::Integrator::LeapFrog lf;
  lf.reinit(sys.deviceData);
  Parallel::TranslationalFreedomRemover trRemover;
  trRemover.reinit (sys.deviceData);
  Parallel::Integrator::BerendsenLeapFrog blf;
  blf.reinit (sys, dt, &interEng, &trRemover);

  blf.TCouple (1.f, 0.1);
  blf.addPcoupleGroup (mdRectBoxDirectionX |
  		       mdRectBoxDirectionY |
  		       mdRectBoxDirectionZ,
  		       -4, 1, 1);
      
  IndexType stFeq = 1;

  printf ("#       1            2              3             4             5             6             7           8           9          10          11     12     13     14\n");
  printf ("#       n         time           nb E           b E           k E             T       total E         Pxx         Pyy         Pzz           P     Lx     Ly     Lz\n");
  for (IndexType i = 0; i < nstep; ++i){
    dst.clearData ();

    blf.oneStep (sys, dst);
    
    dst.copyToHost (hst);
    hst.collectData ();

    if ((i+1)%stFeq == 0 && Parallel::Interface::myRank() == 0){
      printf ("%09d %07e %.7e %.7e %.7e %.7e %.7e %.5e %.5e %.5e %.5e %.3f %.3f %.3f\n",
  	      (i+1),  
  	      (i+1) * dt, 
  	      hst.NonBondedEnergy(),
  	      hst.BondedEnergy(),
  	      hst.kineticEnergy(),
  	      hst.NonBondedEnergy() +
  	      hst.BondedEnergy() +
  	      hst.kineticEnergy(),
  	      hst.kineticEnergy() * 2. / sys.getNumFreedom(),
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

