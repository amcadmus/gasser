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


#define NThreadsPerBlockCell	160
#define NThreadsPerBlockAtom	96

// #define NThreadsPerBlockCell	160
// #define NThreadsPerBlockAtom	96


int main(int argc, char * argv[])
{
  IndexType nstep = 100000;
  IndexType confFeq = 100000;
  IndexType thermoFeq = 100;
  ScalorType rcut = 4.0;
  ScalorType nlistExten = 0.3;
  ScalorType refT = 1.30;
  ScalorType tauT = 1.;
  ScalorType refP = 0.12;
  ScalorType tauP = 1.;
  ScalorType beta = 1.;

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
  ljparam.reinit (1.f, 1.f, 0.f, rcut);
  // ScalorType shift = ljparam.calShiftAtCut ();
  // ljparam.reinit (1.f, 1.f, shift, rcut);
  sysTop.addNonBondedInteraction (Topology::NonBondedInteraction(0, 0, ljparam));
  sysTop.addMolecules (mol, sys.hdata.numAtom);

  sys.initTopology (sysTop);
  sys.initDeviceData ();

  SystemNonBondedInteraction sysNbInter;
  sysNbInter.reinit (sysTop);
  
  ScalorType maxrcut = sysNbInter.maxRcut();
  ScalorType rlist = maxrcut + nlistExten;
  ScalorType rebuildThreshold = 0.5 * nlistExten;
  NeighborList nlist (sysNbInter, sys, rlist, NThreadsPerBlockCell, 10,
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
  ScalorType seed = 1286812148;
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
			refP, beta);
  blpf.addThermostat (thermostat);
  blpf.addBarostat   (barostat);  

  Reshuffle resh (sys, nlist, NThreadsPerBlockCell);
  
  timer.tic(mdTimeTotal);
  resh.calIndexTable (nlist, &timer);
  sys.reshuffle   (resh.getIndexTable(), sys.hdata.numAtom, &timer);
  nlist.reshuffle (resh.getIndexTable(), sys.hdata.numAtom, &timer);  
  
  printf ("# prepare ok, start to run\n");
  printf ("#*     1     2           3         4            5       6         7    8       9  10  11\n");
  printf ("#* nstep  time  nonBondedE  kineticE  temperature  totalE  pressure  box  volume  h0  h1\n");

  // sys.recoverDeviceData (&timer);
  // sys.updateHostFromRecovered (&timer);
  // sys.writeHostDataGro ("confstart.gro", 0, 0.f, &timer);
  try{
    sys.initWriteXtc ("traj.xtc");
    sys.recoverDeviceData (&timer);
    sys.updateHostFromRecovered (&timer);
    sys.writeHostDataXtc (0, 0*dt, &timer);
    for (i = 0; i < nstep; ++i){
      if (i%100 == 0){
	tfremover.remove (sys, &timer);
      }
      if ((i+1) % thermoFeq == 0){
	st.clearDevice();
	blpf.oneStep (sys, st, &timer);
	st.updateHost();
	ScalorType e = st.NonBondedEnergy () + st.kineticEnergy();	
	ScalorType v = sys.box.size.x * sys.box.size.y * sys.box.size.z;
	ScalorType p0 = st.pressure(sys.box);
	ScalorType h0 = v * p0 + e;
	ScalorType h1 = v * refP + e;
	printf ("%09d %07e %.7e %.7e %.7e %.7e %.7e %.7f %.7e %.7e %.7e\n",
		(i+1),  
		(i+1) * dt, 
		st.NonBondedEnergy (),
		st.kineticEnergy(),
		st.kineticEnergy() * 2. / (3. * (double (sys.hdata.numAtom) - 1.)),
		e,
		p0,
		sys.box.size.x,
		v,
		h0,
		h1
	    );
	fflush(stdout);
      }
      else {
	blpf.oneStep (sys, &timer);      
      }
      if ((i+1) % confFeq == 0){
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
