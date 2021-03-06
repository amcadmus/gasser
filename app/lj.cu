#include <stdio.h>
#include "MDSystem_interface.h"
#include "common.h"
#include "BoxGeometry.h"
#include "MDSystem.h"
#include "RandomGenerator.h"
#include "Auxiliary.h"
#include "NeighborList_interface.h"
#include"Statistic.h"
#include "Integrator_interface.h"
#include "InteractionEngine_interface.h"
#include "Reshuffle_interface.h"
#include "Displacement_interface.h"

#include "Topology.h"
#include "SystemBondedInteraction.h"

#include "BondInteraction.h"
#include "NonBondedInteraction.h"


#define NThreadsPerBlockCell	128
#define NThreadsPerBlockAtom	96

int main(int argc, char * argv[])
{
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
  printf ("# setting device to %d\n", atoi(argv[3]));
  cudaSetDevice (atoi(argv[3]));
  checkCUDAError ("set device");

  MDSystem sys;
  sys.initConfig(filename);

  Topology::System sysTop;
  Topology::Molecule mol;
  mol.pushAtom (Topology::Atom (1.0, 0.0, 0));
  LennardJones6_12Parameter ljparam;
  ScalorType rcut = 3.2f;
  ljparam.reinit (1.f, 1.f, 0.f, rcut);
  sysTop.addNonBondedInteraction (Topology::NonBondedInteraction(0, 0, ljparam));
  sysTop.addMolecules (mol, sys.hdata.numAtom);

  sys.initTopology (sysTop);
  sys.initDeviceData ();
  
  SystemNonBondedInteraction sysNbInter;
  sysNbInter.reinit (sysTop);
  
  ScalorType maxrcut = sysNbInter.maxRcut();
  ScalorType nlistExten = 0.3;
  ScalorType rlist = maxrcut + nlistExten;
  CellList clist (sys, rlist, NThreadsPerBlockCell, NThreadsPerBlockAtom);
  NeighborList nlist (sysNbInter, sys, rlist, NThreadsPerBlockAtom, 10.f);
  sys.normalizeDeviceData ();
  clist.rebuild (sys, NULL);
  nlist.rebuild (sys, clist, NULL);
  Displacement_max disp (sys, NThreadsPerBlockAtom);
  disp.recordCoord (sys);
  
  MDStatistic st(sys);
  VelocityVerlet inte_vv (sys, NThreadsPerBlockAtom);
  ScalorType refT = 1.;
  VelocityRescale inte_vr (sys, NThreadsPerBlockAtom, refT, 0.1);
  TranslationalFreedomRemover tfremover (sys, NThreadsPerBlockAtom);
  InteractionEngine inter (sys, NThreadsPerBlockAtom);
  inter.registNonBondedInteraction (sysNbInter);
  
  MDTimer timer;
  unsigned i;
  ScalorType dt = 0.001;
  ScalorType seed = 1;
  RandomGenerator_MT19937::init_genrand (seed);

  Reshuffle resh (sys);
  
  timer.tic(mdTimeTotal);
  if (resh.calIndexTable (clist, &timer)){
    sys.reshuffle   (resh.indexTable, sys.hdata.numAtom, &timer);
    clist.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);  
    nlist.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);  
    disp.reshuffle  (resh.indexTable, sys.hdata.numAtom, &timer);  
  }
  
  printf ("# prepare ok, start to run\n");
  // sys.recoverDeviceData (&timer);
  // sys.updateHostFromRecovered (&timer);
  // sys.writeHostDataGro ("confstart.gro", 0, 0.f, &timer);
  printf ("# prepare ok, start to run\n");
  printf ("#*     1     2           3         4            5       6          7-9\n");
  printf ("#* nstep  time  nonBondedE  kineticE  temperature  totalE  pressurexyz\n");
  try{
    // sys.initWriteXtc ("traj.xtc");
    // sys.recoverDeviceData (&timer);
    // sys.updateHostFromRecovered (&timer);
    // sys.writeHostDataXtc (0, 0*dt, &timer);
    for (i = 0; i < nstep; ++i){
      if (i%10 == 0){
	tfremover.remove (sys, &timer);
      }
      if ((i+1) % 100 == 0){
	st.clearDevice();
	inte_vv.step1 (sys, dt, &timer);
	inter.clearInteraction (sys);
	ScalorType maxdr = disp.calMaxDisplacemant (sys, &timer);
	if (maxdr > nlistExten * 0.5){
	  // printf ("# Rebuild at step %09i ... ", i+1);
	  // fflush(stdout);
	  // rebuild
	  sys.normalizeDeviceData (&timer);
	  disp.recordCoord (sys);
	  clist.rebuild (sys, &timer);
	  inter.applyNonBondedInteraction (sys, clist, rcut, st, &timer);
	  nlist.rebuild (sys, clist, &timer);
	  // printf ("done\n");
	  // fflush(stdout);
	}
	else{
	  inter.applyNonBondedInteraction (sys, nlist, st, NULL, &timer);
	  // inter.applyNonBondedInteraction (sys, rcut, st, &timer);
	}
	inte_vv.step2 (sys, dt, st, &timer);
	st.updateHost();
	printf ("%09d %07e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e \n",
		(i+1),  
		(i+1) * dt, 
		st.nonBondedEnergy(),
		st.kineticEnergy(),
		st.kineticEnergy() / (sys.ddata.numAtom - 1) * 2./3.,
		st.nonBondedEnergy() + st.bondedEnergy() + st.kineticEnergy(),
		st.pressureXX(sys.box),
		st.pressureYY(sys.box),
		st.pressureZZ(sys.box),
		st.pressure(sys.box));
	fflush(stdout);
      }
      else {
	inte_vv.step1 (sys, dt, &timer);
	inter.clearInteraction (sys);
	ScalorType maxdr = disp.calMaxDisplacemant (sys, &timer);
	if (maxdr > nlistExten * 0.5){
	  // printf ("# Rebuild at step %09i ... ", i+1);
	  // fflush(stdout);
	  // rebuild
	  sys.normalizeDeviceData (&timer);
	  disp.recordCoord (sys);
	  clist.rebuild (sys, &timer);
	  inter.applyNonBondedInteraction (sys, clist, rcut, &timer);
	  nlist.rebuild (sys, clist, &timer);
	  // printf ("done\n");
	  // fflush(stdout);
	}
	else{
	  inter.applyNonBondedInteraction (sys, nlist, NULL, &timer);
	  // inter.applyNonBondedInteraction (sys, rcut, &timer);
	}
	inte_vv.step2 (sys, dt, &timer);
      }
      // if (maxdr > nlistExten * 0.5){
      // 	printf ("# Rebuild at step %09i ... ", i+1);
      // 	fflush(stdout);
      // 	// rebuild
      // 	sys.normalizeDeviceData ();
      // 	disp.recordCoord (sys);
      // 	clist.rebuild (sys, &timer);
      // 	nlist.rebuild (sys, clist, &timer);
      // 	printf ("done\n");
      // 	fflush(stdout);
      // }
      // if ((i+1) % 1000 == 0){
      // 	sys.recoverDeviceData (&timer);
      // 	sys.updateHostFromRecovered (&timer);
      // 	sys.writeHostDataXtc (i+1, (i+1)*dt, &timer);
      // }
      if ((i+1) % 100 == 0){
	if (resh.calIndexTable (clist, &timer)){
	  sys.reshuffle   (resh.indexTable, sys.hdata.numAtom, &timer);
	  clist.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);  
	  nlist.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);  
	  disp.reshuffle  (resh.indexTable, sys.hdata.numAtom, &timer);  
	}
      }
    }
    // sys.endWriteXtc();
    // sys.recoverDeviceData (&timer);
    // sys.updateHostFromRecovered (&timer);
    // sys.writeHostDataGro ("confout.gro", nstep, nstep*dt, &timer);
    timer.toc(mdTimeTotal);
    timer.printRecord (stderr);
  }
  catch (MDExcptCuda & e){
    // resh.recoverMDDataToHost (sys, &timer);
    // sys.writeHostDataXtc (i+1, (i+1)*dt, &timer);
    timer.toc(mdTimeTotal);
    timer.printRecord (stderr);
    return 1;
  }
  catch (MDException &e){
    fprintf (stderr, "%s\n", e.what());
    return 1;
  }
  
  
  return 0;
}

  
