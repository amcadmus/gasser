#include <stdio.h>
#include "common.h"
#include "BoxGeometry.h"
#include "MDSystem.h"
#include "RandomGenerator.h"
#include "Auxiliary.h"
#include "NeighborList_interface.h"
#include"Statistic.h"
#include "Integrator_interface.h"
#include "InteractionEngine_interface.h"
#include "tmp.h"
#include "Reshuffle_interface.h"

#include "MDSystem_interface.h"

#include "Topology.h"
#include "SystemBondedInteraction.h"

#include "BondInteraction.h"
#include "NonBondedInteraction.h"


#define NThreadsPerBlockCell	64
#define NThreadsPerBlockAtom	64

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
  ljparam.reinit (1.f, 1.f, 0.f, 3.2f);
  mol.addNonBondedInteraction (Topology::NonBondedInteraction(0, 0, ljparam));
  sysTop.addMolecules (mol, sys.hdata.numAtom);

  sys.initTopology (sysTop);
  sys.initDeviceData ();

  SystemNonBondedInteraction sysNbInter;
  sysNbInter.reinit (sysTop);
  
  ScalorType maxrcut = sysNbInter.maxRcut();
  ScalorType nlistExten = 0.0;
  ScalorType rlist = maxrcut + nlistExten;
  NeighborList nlist (sysNbInter, sys, rlist, NThreadsPerBlockCell, 5,
		      RectangularBoxGeometry::mdRectBoxDirectionX |
		      RectangularBoxGeometry::mdRectBoxDirectionY |
		      RectangularBoxGeometry::mdRectBoxDirectionZ,
		      1);
  nlist.buildCellList(sys);

  MDStatistic st(sys);
  VelocityVerlet inte_vv (sys, NThreadsPerBlockAtom);
  ScalorType refT = 0.9977411970749;
  VelocityRescale inte_vr (sys, NThreadsPerBlockAtom, refT, 0.1);
  TranslationalFreedomRemover tfremover (sys, NThreadsPerBlockAtom);
  InteractionEngine_interface inter (sys, NThreadsPerBlockAtom);
  inter.registNonBondedInteraction (sysNbInter);
  
  MDTimer timer;
  unsigned i;
  ScalorType dt = 0.001;
  ScalorType seed = 1;
  RandomGenerator_MT19937::init_genrand (seed);

  // Reshuffle resh (sys, nlist, NThreadsPerBlockCell);
  
  timer.tic(mdTimeTotal);
  // resh.calIndexTable (nlist, &timer);
  // sys.reshuffle   (resh.getIndexTable(), sys.hdata.numAtom, &timer);
  // nlist.reshuffle (resh.getIndexTable(), sys.hdata.numAtom, &timer);  
  
  printf ("# prepare ok, start to run\n");
  // sys.recoverDeviceData (&timer);
  // sys.updateHostFromRecovered (&timer);
  // sys.writeHostDataGro ("confstart.gro", 0, 0.f, &timer);
  try{
    // sys.initWriteXtc ("traj.xtc");
    // sys.recoverDeviceData (&timer);
    // sys.updateHostFromRecovered (&timer);
    // sys.writeHostDataXtc (0, 0*dt, &timer);
    for (i = 0; i < nstep; ++i){
      if (i%10 == 0){
	tfremover.remove (sys, &timer);
      }
      if ((i+1) % 10 == 0){
	st.clearDevice();
	inte_vv.step1 (sys, dt, &timer);
	inter.clearInteraction (sys);
	inter.applyNonBondedInteractionCell (sys, nlist, st, &timer);
	inte_vv.step2 (sys, dt, st, &timer);
	st.updateHost();
	printf ("%09d %07e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e\n",
		(i+1),  
		(i+1) * dt, 
		st.getStatistic(mdStatisticNonBondedPotential),
		st.getStatistic(mdStatisticBondedPotential),
		st.kineticEnergy(),
		st.getStatistic(mdStatisticNonBondedPotential) +
		st.getStatistic(mdStatisticBondedPotential) +
		st.kineticEnergy(),
		st.pressureXX(),
		st.pressureYY(),
		st.pressureZZ(),
		st.pressure());
	fflush(stdout);
      }
      else {
	inte_vv.step1 (sys, dt, &timer);
	inter.clearInteraction (sys);
	inter.applyNonBondedInteractionCell (sys, nlist, &timer);
	inte_vv.step2 (sys, dt, &timer);
      }
      nlist.reBuildCellList(sys, &timer);

      // if (nlist.judgeRebuild(sys, 0.5 * nlistExten, &timer)){
      // 	// printf ("# Rebuild at step %09i ... ", i+1);
      // 	// fflush(stdout);
      // 	nlist.reBuildCellList(sys, &timer);
      // 	// printf ("done\n");
      // 	// fflush(stdout);
      // }
      // if ((i+1) % 1000 == 0){
      // 	sys.recoverDeviceData (&timer);
      // 	sys.updateHostFromRecovered (&timer);
      // 	sys.writeHostDataXtc (i+1, (i+1)*dt, &timer);
      // }
      if ((i+1) % 100 == 0){
	// resh.calIndexTable (nlist, &timer);
	// sys.reshuffle   (resh.getIndexTable(), sys.hdata.numAtom, &timer);
	// nlist.reshuffle (resh.getIndexTable(), sys.hdata.numAtom, &timer);  
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

  
