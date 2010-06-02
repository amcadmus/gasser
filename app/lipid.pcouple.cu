#include <stdio.h>
#include "MDSystem_interface.h"
#include "common.h"
#include "BoxGeometry.h"
#include "MDSystem.h"
#include "RandomGenerator.h"
#include "Auxiliary.h"
#include "NeighborList_interface.h"
#include "Statistic.h"
#include "Integrator_interface.h"
#include "InteractionEngine_interface.h"
#include "tmp.h"
#include "Reshuffle_interface.h"

#include "Topology.h"
#include "SystemBondedInteraction.h"

#include "BondInteraction.h"
#include "NonBondedInteraction.h"


#define NThreadsPerBlockCell	192
#define NThreadsPerBlockAtom	16

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
  mol.pushAtom (Topology::Atom (1.0, 0.0, 1));
  mol.pushAtom (Topology::Atom (1.0, 0.0, 1));
  // nonbonded interactions
  CosTailParameter cosparam;
  cosparam.reinit (1.f, 0.95f, 0.f);
  sysTop.addNonBondedInteraction (Topology::NonBondedInteraction(0, 0, cosparam));  
  cosparam.reinit (1.f, 0.975f, 0.f);
  sysTop.addNonBondedInteraction (Topology::NonBondedInteraction(0, 1, cosparam));  
  cosparam.reinit (1.f, 1.f, 1.6f);
  sysTop.addNonBondedInteraction (Topology::NonBondedInteraction(1, 1, cosparam));
  // bonded interactions
  HarmonicSpringParameter hsparam;
  FENEParameter feneparam;
  hsparam.reinit (10.f, 4.f);
  feneparam.reinit (30.f, 1.5f);
  mol.addBond (Topology::Bond (0, 1, feneparam));
  mol.addBond (Topology::Bond (1, 2, feneparam));
  mol.addBond (Topology::Bond (0, 2, hsparam));
  sysTop.addMolecules (mol, sys.hdata.numAtom/3);

  sys.initTopology (sysTop);
  sys.initDeviceData ();
  
  SystemBondedInteraction sysBdInter;
  sysBdInter.reinit (sysTop);

  BondedInteractionList bdInterList;
  bdInterList.reinit (sys, sysTop, sysBdInter);

  SystemNonBondedInteraction sysNbInter;
  sysNbInter.reinit (sysTop);
  
  ScalorType maxrcut = sysNbInter.maxRcut();
  ScalorType nlistExten = 0.5;
  ScalorType rlist = maxrcut + nlistExten;
  ScalorType rebuildThreshold = 0.5 * nlistExten;
  NeighborList nlist (sysNbInter, sys, rlist, NThreadsPerBlockCell, 200,
		      RectangularBoxGeometry::mdRectBoxDirectionX |
		      RectangularBoxGeometry::mdRectBoxDirectionY);
  nlist.build(sys);
  MDStatistic st(sys);
  TranslationalFreedomRemover tfremover (sys, NThreadsPerBlockAtom);
  InteractionEngine_interface inter (sys, NThreadsPerBlockAtom);
  inter.registNonBondedInteraction (sysNbInter);
  inter.registBondedInteraction    (sysBdInter);
  
  MDTimer timer;
  unsigned i;
  ScalorType dt = 0.005;
  ScalorType seed = 1;
  RandomGenerator_MT19937::init_genrand (seed);

  BerendsenLeapFrog blpf (sys, NThreadsPerBlockAtom, dt,
			  inter,
			  nlist,
			  rebuildThreshold,
			  &bdInterList);
  blpf.TCouple (0.9977411970749, 0.1);
  blpf.addPcoupleGroup (PCoupleX | PCoupleY,
  			0., 1, 1);
  
  Reshuffle resh (sys, nlist, NThreadsPerBlockCell);

  timer.tic(mdTimeTotal);
  if (resh.calIndexTable (nlist, &timer)){
    sys.reshuffle   (resh.getIndexTable(), sys.hdata.numAtom, &timer);
    nlist.reshuffle (resh.getIndexTable(), sys.hdata.numAtom, &timer);
    bdInterList.reshuffle (resh.getIndexTable(), sys.hdata.numAtom, &timer);
  }
  
  printf ("# prepare ok, start to run\n");
  try{
    timer.tic(mdTimeTotal);
    sys.initWriteXtc ("traj.xtc");
    sys.recoverDeviceData (&timer);
    sys.updateHostFromRecovered (&timer);
    sys.writeHostDataXtc (0, 0*dt, &timer);
    for (i = 0; i < nstep; ++i){
      if (i%100 == 0){
	tfremover.remove (sys, &timer);
      }	  
      if ((i+1) % 10 == 0){
	st.clearDevice();
	blpf.oneStep (sys, st, &timer);
	st.updateHost();
	printf ("%09d %07e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.3f %.3f %.3f\n",
		(i+1),  
		(i+1) * dt, 
		st.getStatistic(mdStatisticNonBondedPotential),
		st.getStatistic(mdStatisticBondedPotential),
		st.kineticEnergy(),
		st.getStatistic(mdStatisticNonBondedPotential) +
		st.getStatistic(mdStatisticBondedPotential) +
		st.kineticEnergy(),
		st.pressureXX(sys.box),
		st.pressureYY(sys.box),
		st.pressureZZ(sys.box),
		st.pressure(sys.box),
		sys.box.size.x,
		sys.box.size.y,
		sys.box.size.z);
	fflush(stdout);
      }
      else {
	blpf.oneStep (sys, &timer);      
      }
      if ((i+1) % 1000 == 0){
      	sys.recoverDeviceData (&timer);
      	sys.updateHostFromRecovered (&timer);
      	sys.writeHostDataXtc (i+1, (i+1)*dt, &timer);
      }
      if ((i+1) % 100 == 0){
	if (resh.calIndexTable (nlist, &timer)){
	  sys.reshuffle   (resh.getIndexTable(), sys.hdata.numAtom, &timer);
	  nlist.reshuffle (resh.getIndexTable(), sys.hdata.numAtom, &timer);
	  bdInterList.reshuffle (resh.getIndexTable(), sys.hdata.numAtom, &timer);
	}
      }
    }
    sys.endWriteXtc();
    sys.recoverDeviceData (&timer);
    sys.updateHostFromRecovered (&timer);
    sys.writeHostDataGro ("confout.gro", nstep, nstep*dt, &timer);
    timer.toc(mdTimeTotal);
    timer.printRecord (stderr);
  }
  catch (MDExcptCuda & e){
    sys.writeHostDataXtc (i+1, (i+1)*dt, &timer);
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

