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


#define NThreadsPerBlockCell	512
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
  ScalorType rcut = 5.f;
  ljparam.reinit (1.f, 1.f, 0.f, rcut);
  // NonBondedInteractionParameter * nbp (&ljparam);
  Topology::NonBondedInteraction nb00 (0, 0, ljparam);
  // printf ("# %f %f %f\n",
  // 	  ljparam.energyCorrection(1.8),
  // 	  nbp->energyCorrection(1.8),
  // 	  nb00.energyCorrection(1.8));
  sysTop.addNonBondedInteraction (nb00);
  sysTop.addMolecules (mol, sys.hdata.numAtom);

  sys.initTopology (sysTop);
  sys.initDeviceData ();

  SystemNonBondedInteraction sysNbInter;
  sysNbInter.reinit (sysTop);
  
  ScalorType maxrcut = sysNbInter.maxRcut();
  ScalorType nlistExten = 0.2;
  ScalorType rlist = maxrcut + nlistExten;
  NeighborList nlist (sysNbInter, sys, rlist, NThreadsPerBlockCell, 2,
		      RectangularBoxGeometry::mdRectBoxDirectionX |
		      RectangularBoxGeometry::mdRectBoxDirectionY |
		      RectangularBoxGeometry::mdRectBoxDirectionZ);
  nlist.build(sys);
  MDStatistic st(sys);
  VelocityVerlet inte_vv (sys, NThreadsPerBlockAtom);
  // ScalorType refT = 1.45;
  ScalorType refT = 1.45;
  VelocityRescale inte_vr (sys, NThreadsPerBlockAtom, refT, 0.1);
  TranslationalFreedomRemover tfremover (sys, NThreadsPerBlockAtom);
  InteractionEngine_interface inter (sys, NThreadsPerBlockAtom);
  inter.registNonBondedInteraction (sysNbInter);

  WidomTestParticleInsertion_NVT2 widom;
  widom.reinit (refT, 2, sys.box, 0, sysNbInter);
  // widom.reinit (refT, 2, 0, sysNbInter);
		
  MDTimer timer;
  unsigned i;
  ScalorType dt = 0.002;
  ScalorType seed = 1;
  RandomGenerator_MT19937::init_genrand (seed);

  Reshuffle resh (sys, nlist, NThreadsPerBlockCell);
  
  timer.tic(mdTimeTotal);
  if (resh.calIndexTable (nlist, &timer)){
    sys.reshuffle   (resh.getIndexTable(), sys.hdata.numAtom, &timer);
    nlist.reshuffle (resh.getIndexTable(), sys.hdata.numAtom, &timer);  
  }
  
  printf ("# prepare ok, start to run\n");
  // sys.recoverDeviceData (&timer);
  // sys.updateHostFromRecovered (&timer);
  // sys.writeHostDataGro ("confstart.gro", 0, 0.f, &timer);
  printf ("# prepare ok, start to run\n");
  printf ("#*     1     2           3         4            5       6          7-9        10  11\n");
  printf ("#* nstep  time  nonBondedE  kineticE  temperature  totalE  pressurexyz  pressure  mu\n");
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
	widom.generateTestCoords (sys);
	st.clearDevice();
	inte_vr.step1 (sys, dt, &timer);
	inter.clearInteraction (sys);
	inter.applyNonBondedInteraction (sys, nlist, st, &timer);
	inter.calculateWidomDeltaEnergy (sys, nlist, widom, &timer);
	inte_vr.step2 (sys, dt, st, &timer);
	st.updateHost();
	printf ("%09d %07e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e\n",
		(i+1),  
		(i+1) * dt, 
		st.NonBondedEnergy(),
		st.kineticEnergy(),
		st.kineticEnergy() / (sys.ddata.numAtom - 1) * 2./3.,
		st.NonBondedEnergy() + st.kineticEnergy(),
		st.pressureXX(sys.box),
		st.pressureYY(sys.box),
		st.pressureZZ(sys.box),
		st.pressure(sys.box),
		widom.expMu());	
	fflush(stdout);
	
      }
      else {
	inte_vr.step1 (sys, dt, &timer);
	inter.clearInteraction (sys);
	inter.applyNonBondedInteraction (sys, nlist, &timer);
	inte_vr.step2 (sys, dt, &timer);
      }
      if (nlist.judgeRebuild(sys, 0.5 * nlistExten, &timer)){
	// printf ("# Rebuild at step %09i ... ", i+1);
	// fflush(stdout);
	nlist.reBuild(sys, &timer);
	// printf ("done\n");
	// fflush(stdout);
      }
      // if ((i+1) % 1000 == 0){
      // 	sys.recoverDeviceData (&timer);
      // 	sys.updateHostFromRecovered (&timer);
      // 	sys.writeHostDataXtc (i+1, (i+1)*dt, &timer);
      // }
      if ((i+1) % 100 == 0){
      	if(resh.calIndexTable (nlist, &timer)){
	  sys.reshuffle   (resh.getIndexTable(), sys.hdata.numAtom, &timer);
	  nlist.reshuffle (resh.getIndexTable(), sys.hdata.numAtom, &timer);
	}
      }
    }
    // sys.endWriteXtc();
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

  
