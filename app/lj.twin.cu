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
#include "Displacement_interface.h"

#include "Topology.h"
#include "SystemBondedInteraction.h"

#include "BondInteraction.h"
#include "NonBondedInteraction.h"


#define NThreadsPerBlockCell	32
#define NThreadsPerBlockAtom	4

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
  ScalorType rcut1 = 2.f;
  ScalorType rcut2 = 6.f;
  ScalorType rThreshold = 0.3f;
  ljparam.reinit (1.f, 1.f, 0.f, rcut2);
  sysTop.addNonBondedInteraction (Topology::NonBondedInteraction(0, 0, ljparam));
  sysTop.addMolecules (mol, sys.hdata.numAtom);

  sys.initTopology (sysTop);
  sys.initDeviceData ();
  
  SystemNonBondedInteraction sysNbInter;
  sysNbInter.reinit (sysTop);
  CellList clist1 (sys, rcut1, NThreadsPerBlockCell, NThreadsPerBlockAtom);
  CellList clist2 (sys, rcut2, NThreadsPerBlockCell, NThreadsPerBlockAtom);
  NeighborList nlist (sysNbInter, sys, rcut1, NThreadsPerBlockAtom, 10.f);
  TwinRangeCorrectionRecorder tcrec (sys, NThreadsPerBlockAtom);
  
  sys.normalizeDeviceData ();
  clist1.rebuild (sys);
  clist2.rebuild (sys);
  nlist.rebuild (sys, clist1);
  Displacement_mean disp (sys, NThreadsPerBlockAtom);
  disp.recordCoord (sys);
  Reshuffle resh (sys);  
  
  MDStatistic st (sys);
  MDStatistic last_st (sys);
  TranslationalFreedomRemover tfremover (sys, NThreadsPerBlockAtom);
  InteractionEngine inter (sys, NThreadsPerBlockAtom);
  inter.registNonBondedInteraction (sysNbInter);
  inter.clearInteraction (sys);
  inter.applyNonBondedInteraction (sys, nlist, st);
  inter.calTwinRangeCorrection (sys, clist2, rcut1, rcut2, tcrec);
  tcrec.correct (sys, st);

  MDTimer timer;
  unsigned i;
  ScalorType dt = 0.001;
  ScalorType seed = 1;
  RandomGenerator_MT19937::init_genrand (seed);
  ScalorType refT = 1.34;
  ScalorType tauT = 1.;
  ScalorType refP = 0.1420;
  ScalorType tauP = 1.;
  ScalorType betaP = 1.;
  Thermostat_NoseHoover thermostat;
  thermostat.reinit (refT, dt, tauT, sys.ddata.numAtom * 3 - 3);
  Barostat_ParrinelloRahman barostat;
  barostat.reinit (dt, tauP, sys.box);
  barostat.assignGroup (mdRectBoxDirectionX |
  			mdRectBoxDirectionY |
  			mdRectBoxDirectionZ,
  			refP, betaP);
  LeapFrog_TPCouple_VCouple blpf (sys, NThreadsPerBlockAtom);
  blpf.addThermostat (thermostat);
  blpf.addBarostat   (barostat);

  timer.tic(mdTimeTotal);

  if (resh.calIndexTable (clist1)){
    sys.reshuffle   (resh.indexTable, sys.hdata.numAtom, &timer);
    clist1.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);
    clist2.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);  
    nlist.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);  
    disp.reshuffle  (resh.indexTable, sys.hdata.numAtom, &timer);
    tcrec.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);
  }
  
  printf ("# prepare ok, start to run\n");
  // sys.recoverDeviceData (&timer);
  // sys.updateHostFromRecovered (&timer);
  // sys.writeHostDataGro ("confstart.gro", 0, 0.f, &timer);
  printf ("# prepare ok, start to run\n");
  printf ("#*     1     2           3         4            5       6         7     8       9  10  11\n");
  printf ("#* nstep  time  nonBondedE  kineticE  temperature  totalE  pressure  boxl  volume  h0  h1\n");
  try{
    // sys.initWriteXtc ("traj.xtc");
    // sys.recoverDeviceData (&timer);
    // sys.updateHostFromRecovered (&timer);
    // sys.writeHostDataXtc (0, 0*dt, &timer);
    for (i = 0; i < nstep; ++i){
      last_st.deviceCopy (st);
      if (i%10 == 0){
	tfremover.remove (sys, &timer);
      }
      st.clearDevice();
      blpf.oneStep (sys, dt, last_st, st, &timer);
      ScalorType meandr = disp.calMeanDisplacemant (sys, &timer);
      if (meandr > rThreshold * 0.5){
	// // printf ("# Rebuild at step %09i ... \n", i+1);
	// // fflush(stdout);
	// // rebuild
	sys.normalizeDeviceData (&timer);
	disp.recordCoord (sys, &timer);
	clist1.rebuild (sys, &timer);
	nlist.rebuild (sys, clist1, &timer);
	clist2.rebuild (sys, &timer);
	inter.calTwinRangeCorrection (sys, clist2, rcut1, rcut2, tcrec, &timer);
	// // printf ("done\n");
	// // fflush(stdout);
      }
      inter.clearInteraction (sys);
      inter.applyNonBondedInteraction (sys, nlist, st, &timer);
      tcrec.correct (sys, st, &timer);

      timer.tic (mdTimeDataIO);
      if ((i+1) % 1 == 0){
      	st.updateHost();
      	ScalorType ep = st.nonBondedEnergy ();
      	ScalorType ek = st.kineticEnergy();
      	ScalorType e = ep + ek;
      	ScalorType v = sys.box.size.x * sys.box.size.y * sys.box.size.z;
      	ScalorType p0 = st.pressure(sys.box);
      	ScalorType h0 = v * p0 + e;
      	ScalorType h1 = v * refP + e;
      	printf ("%09d %07e %.7e %.7e %.7e %.7e %.7e %.7f %.7e %.7e %.7e\n",
      		(i+1),  
      		(i+1) * dt, 
      		ep,
      		ek, 
      		ek * 2. / (3. * (double (sys.hdata.numAtom) - 1.)),
      		e,
      		p0,
      		sys.box.size.x,
      		v,
      		h0,
      		h1
      	    );
      	fflush(stdout);
      }
      timer.toc (mdTimeDataIO);
      if ((i+1) % 100 == 0){
      	if (resh.calIndexTable (clist1, &timer)){
      	  sys.reshuffle   (resh.indexTable, sys.hdata.numAtom, &timer);
      	  clist1.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);  
      	  clist2.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);  
      	  nlist.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);  
      	  disp.reshuffle  (resh.indexTable, sys.hdata.numAtom, &timer);  
	  tcrec.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);
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

  
