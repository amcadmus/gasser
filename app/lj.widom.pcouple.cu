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


#define NThreadsPerBlockCell	160
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
  ScalorType rcut = 4.6;
  // ljparam.reinit (1.f, 1.f, 1.0f, rcut);
  ljparam.reinit (1.f, 1.f, .0f, rcut);

  printf ("# shift is %f\n", ljparam.shiftAtCut());
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
  ScalorType nlistExten = 0.3;
  ScalorType rebuildThreshold = 0.5 * nlistExten;
  ScalorType rlist = maxrcut + nlistExten;
  NeighborList nlist (sysNbInter, sys, rlist, NThreadsPerBlockCell, 2,
		      RectangularBoxGeometry::mdRectBoxDirectionX |
		      RectangularBoxGeometry::mdRectBoxDirectionY |
		      RectangularBoxGeometry::mdRectBoxDirectionZ, 1);
  nlist.build(sys);
  MDStatistic st(sys);
  TranslationalFreedomRemover tfremover (sys, NThreadsPerBlockAtom);
  InteractionEngine_interface inter (sys, NThreadsPerBlockAtom);
  inter.registNonBondedInteraction (sysNbInter);
		
  MDTimer timer;
  unsigned i;
  ScalorType dt = 0.002;
  IndexType energyFeq = 10;
  ScalorType seed = 1;
  ScalorType refT = 1.3;
  ScalorType refP = 0.12;
  ScalorType tauT = 0.1;
  ScalorType tauP = 1.;
  ScalorType betaP = 1.;
  RandomGenerator_MT19937::init_genrand (seed);

  LeapFrog_TPCouple_VCouple blpf (sys,
				  NThreadsPerBlockAtom,
				  dt,
				  inter,
				  nlist,
				  rebuildThreshold);
  Thermostat_NoseHoover thermostat;
  thermostat.reinit (refT, dt, tauT, sys.ddata.numAtom * 3 - 3);
  Barostat_ParrinelloRahman barostat;
  barostat.reinit (dt, tauP, sys.box);
  barostat.assignGroup (mdRectBoxDirectionX |
                        mdRectBoxDirectionY |
                        mdRectBoxDirectionZ,
                        refP, betaP);
  blpf.addThermostat (thermostat);
  blpf.addBarostat   (barostat);

  WidomTestParticleInsertion_NPT widom;
  IndexType ntest = 20;
  widom.reinit (refT, refP, ntest, 0, sysNbInter);
  
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
  printf ("#*     1     2           3         4            5       6       7         8    9  10\n");
  printf ("#* nstep  time  nonBondedE  kineticE  temperature  totalE  virial  pressure  box  mu\n");
  try{
    // sys.initWriteXtc ("traj.xtc");
    // sys.recoverDeviceData (&timer);
    // sys.updateHostFromRecovered (&timer);
    // sys.writeHostDataXtc (0, 0*dt, &timer);
    for (i = 0; i < nstep; ++i){
      if (i%10 == 0){
	tfremover.remove (sys, &timer);
      }
      if ((i+1) % energyFeq == 0){
	widom.generateTestCoords (sys);
	st.clearDevice();
	blpf.oneStep (sys, st, &timer);
	st.updateHost();
	inter.calculateWidomDeltaEnergy (sys, nlist, widom, &timer);
	printf ("%09d %07e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e\n",
		(i+1),  
		(i+1) * dt, 
		st.NonBondedEnergy() / sys.ddata.numAtom,
		st.kineticEnergy() / sys.ddata.numAtom,
		st.kineticEnergy() / (sys.ddata.numAtom - 1) * 2./3.,
		(st.NonBondedEnergy() + st.kineticEnergy()) / sys.ddata.numAtom,
		- st.virial() / sys.ddata.numAtom,
		st.pressure(sys.box),
		sys.box.size.x,
		widom.expMu()
	    );	
	fflush(stdout);
      }
      else {
	blpf.oneStep (sys, &timer);      
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

  
