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
#include "tmp.h"
#include "Reshuffle_interface.h"


#include "Topology.h"
#include "SystemBondedInteraction.h"

#include "BondInteraction.h"
#include "NonBondedInteraction.h"


#define NThreadsPerBlockCell	128
#define NThreadsPerBlockAtom	96

// #define NThreadsPerBlockCell	16
// #define NThreadsPerBlockAtom	16


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
  ljparam.reinit (1.f, 1.f, 0.f, 2.5f);
  sysTop.addNonBondedInteraction (Topology::NonBondedInteraction(0, 0, ljparam));
  sysTop.addMolecules (mol, sys.hdata.numAtom);

  sys.initTopology (sysTop);
  sys.initDeviceData ();

  SystemNonBondedInteraction sysNbInter;
  sysNbInter.reinit (sysTop);
  
  ScalorType maxrcut = sysNbInter.maxRcut();
  ScalorType nlistExten = 0.4;
  ScalorType rlist = maxrcut + nlistExten;
  ScalorType rebuildThreshold = 0.5 * nlistExten;
  NeighborList nlist (sysNbInter, sys, rlist, NThreadsPerBlockCell, 5,
		      RectangularBoxGeometry::mdRectBoxDirectionX |
		      RectangularBoxGeometry::mdRectBoxDirectionY |
		      RectangularBoxGeometry::mdRectBoxDirectionZ);
  nlist.build(sys);
  MDStatistic st(sys);
  TranslationalFreedomRemover tfremover (sys, NThreadsPerBlockAtom);
  InteractionEngine_interface inter (sys, NThreadsPerBlockAtom);
  inter.registNonBondedInteraction (sysNbInter);

  MDTimer timer;
  unsigned i;
  ScalorType dt = 0.002;
  ScalorType seed = 1;
  RandomGenerator_MT19937::init_genrand (seed);

  // LeapFrog_TPCouple lpfrog (sys, NThreadsPerBlockAtom);
  // LeapFrog_TPCouple_Rescale blpf (sys,
  // 				  NThreadsPerBlockAtom,
  // 				  dt,
  // 				  inter,
  // 				  nlist,
  // 				  rebuildThreshold);
  LeapFrog_TPCouple_VCouple blpf (sys,
				  NThreadsPerBlockAtom,
				  dt,
				  inter,
				  nlist,
				  rebuildThreshold);
  // blpf.TCouple (1, 0.1);
  Thermostat_NoseHoover thermostat;
  // ScalorType refT = 0.9977411970749;
  ScalorType refT = 1.3;
  thermostat.reinit (refT, dt, 1, sys.ddata.numAtom * 3 - 3);
  blpf.addThermostat (thermostat);
  Barostat_ParrinelloRahman barostat;
  barostat.reinit (dt, 1, sys.box);
  // barostat.assignGroup (mdRectBoxDirectionX |
  // 			mdRectBoxDirectionY |
  // 			mdRectBoxDirectionZ,
  // 			0.07, 1);
  barostat.assignGroup (mdRectBoxDirectionX |
  			mdRectBoxDirectionZ,
  			0.07, 1);
  blpf.addBarostat (barostat);
  
  // blpf.addPcoupleGroup (PCoupleX | PCoupleY | PCoupleZ,
  // 			0, 1, 10);
  // blpf.addPcoupleGroup (PCoupleZ,
  // 			4., 10, 1);
  
  // blpf.addPcoupleGroup (PCoupleX | PCoupleY ,
  // 			0, 1, 10);
  // blpf.addPcoupleGroup (PCoupleZ,
  // 			0, 1, 1);

  Reshuffle resh (sys, nlist, NThreadsPerBlockCell);
  
  timer.tic(mdTimeTotal);
  resh.calIndexTable (nlist, &timer);
  sys.reshuffle   (resh.getIndexTable(), sys.hdata.numAtom, &timer);
  nlist.reshuffle (resh.getIndexTable(), sys.hdata.numAtom, &timer);  
  
  printf ("# prepare ok, start to run\n");
  // sys.recoverDeviceData (&timer);
  // sys.updateHostFromRecovered (&timer);
  // sys.writeHostDataGro ("confstart.gro", 0, 0.f, &timer);
  printf ("# prepare ok, start to run\n");
  printf ("#*     1     2           3         4            5       6          7-9       10   11-13      14\n");
  printf ("#* nstep  time  nonBondedE  kineticE  temperature  totalE  pressurexyz pressure  boxxyz  volume\n");
  try{
    sys.initWriteXtc ("traj.xtc");
    sys.recoverDeviceData (&timer);
    sys.updateHostFromRecovered (&timer);
    sys.writeHostDataXtc (0, 0*dt, &timer);
    for (i = 0; i < nstep; ++i){
      if (i%100 == 0){
	tfremover.remove (sys, &timer);
      }
      if ((i+1) % 100 == 0){
	st.clearDevice();
	blpf.oneStep (sys, st, &timer);
	st.updateHost();
	printf ("%09d %07e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.3f %.3f %.3f %.7e\n",
		(i+1),  
		(i+1) * dt, 
		st.getStatistic(mdStatisticNonBondedPotential),
		st.kineticEnergy(),
		st.kineticEnergy() / (sys.ddata.numAtom - 1) * 2./3.,
		st.getStatistic(mdStatisticNonBondedPotential) +
		st.getStatistic(mdStatisticBondedPotential) +
		st.kineticEnergy(),
		st.pressureXX(sys.box),
		st.pressureYY(sys.box),
		st.pressureZZ(sys.box),
		st.pressure(sys.box),
		sys.box.size.x,
		sys.box.size.y,
		sys.box.size.z,
		sys.box.size.x * sys.box.size.y * sys.box.size.z);
	fflush(stdout);
      }
      else {
	blpf.oneStep (sys, &timer);      
      }
      if ((i+1) % 100000 == 0){
      	sys.recoverDeviceData (&timer);
      	sys.updateHostFromRecovered (&timer);
      	sys.writeHostDataXtc (i+1, (i+1)*dt, &timer);
      }
      if ((i+1) % 100 == 0){
	resh.calIndexTable (nlist, &timer);
	sys.reshuffle   (resh.getIndexTable(), sys.hdata.numAtom, &timer);
	nlist.reshuffle (resh.getIndexTable(), sys.hdata.numAtom, &timer);  
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