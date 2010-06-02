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
  Parallel::Interface::initEnvironment (32, "Tesla C1060");
  
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
  if (Parallel::Interface::myRank() == 0) sysBdInter.printEverything ();
  
  sys.init (filename, sysTop, sysNbInter.maxRcut());
  sys.redistribute ();

  Parallel::InteractionEngine interEng (sys.deviceData);
  interEng.registNonBondedInteraction (sysNbInter);
  interEng.registBondedInteraction    (sysBdInter);
  interEng.clearInteraction (sys.deviceData);
  Parallel::DeviceCellRelation relation;
  relation.rebuild (sys.deviceData);
  Parallel::DeviceCellRelation relation_buildBdList;
  Parallel::SubCellList ghost, innerShell;
  sys.deviceData.buildSubListGhostCell (ghost);
  sys.deviceData.buildSubListInnerShell (innerShell);
  ghost.add (innerShell);
  relation_buildBdList.rebuild (sys.deviceData, innerShell, ghost);
  Parallel::DeviceBondList dbdlist;
  dbdlist.reinit (sys.deviceData);
  buildDeviceBondList (sys.deviceData, relation, dbdlist);
  Parallel::HostStatistic hst;
  Parallel::DeviceStatistic dst;
  // hst.reinit (sys.localHostData);
  // dst.reinit (sys.deviceData);
  
  Parallel::Integrator::VelocityVerlet vv (sys.deviceData);
  Parallel::TranslationalFreedomRemover trRemover (sys.deviceData);
  ScalorType dt = 0.0001;

  sys.globalHostData.initWriteData_xtcFile ("traj.xtc");
  // printf ("# numCell: %d %d %d\n",
  // 	  sys.deviceData.getNumCell().x,
  // 	  sys.deviceData.getNumCell().y,
  // 	  sys.deviceData.getNumCell().z);
  // printf ("# frameUp: %f %f %f\n",
  // 	  sys.deviceData.getFrameUp().x,
  // 	  sys.deviceData.getFrameUp().y,
  // 	  sys.deviceData.getFrameUp().z);
  // printf ("# frameLow: %f %f %f\n",
  // 	  sys.deviceData.getFrameLow().x,
  // 	  sys.deviceData.getFrameLow().y,
  // 	  sys.deviceData.getFrameLow().z);

  if (Parallel::Interface::myRank() == 0){
    printf ("#       1            2              3             4             5             6             7           8           9          10          11\n");
    printf ("#       n         time           nb E           b E           k E             T       total E         Pxx         Pyy         Pzz           P\n");
  }
  
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
    DeviceTimer::tic (item_ClearInteraction);
    interEng.clearInteraction (sys.deviceData);
    DeviceTimer::toc (item_ClearInteraction);
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
      interEng.applyBondedInteraction (sys.deviceData, dbdlist);
      DeviceTimer::toc (item_BondedInteraction);
    }
    HostTimer::tic (item_TransferGhost);
    sys.clearGhost ();
    HostTimer::toc (item_TransferGhost);
    DeviceTimer::tic (item_Integrate);
    vv.step2 (sys.deviceData, dt, dst);
    DeviceTimer::toc (item_Integrate);

    DeviceTimer::tic (item_BuildCellList);
    sys.deviceData.rebuild (dbdlist);
    DeviceTimer::toc (item_BuildCellList);
    HostTimer::tic (item_Redistribute);
    sys.redistribute ();
    HostTimer::toc (item_Redistribute);
    
    dst.copyToHost (hst);
    hst.collectData ();

    if ((i+1)%stFeq == 0 && Parallel::Interface::myRank() == 0){
      printf ("%09d %07e %.7e %.7e %.7e %.7e %.7e %.5e %.5e %.5e %.5e\n",
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
  	      hst.pressure(sys.deviceData.getGlobalBox()));
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
