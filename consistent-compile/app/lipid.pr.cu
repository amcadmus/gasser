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

#define NThreadsPerBlockCell	256
#define NThreadsPerBlockAtom	96

int main(int argc, char * argv[])
{
  IndexType nstep = 1000000;
  IndexType startStep = 0;
  
  char * filename;

  if (argc < 4){
    printf ("Usage:\n%s conf.gro nstep device [startStep]\n", argv[0]);
    return 1;
  }
  if (argc != 1){
    nstep = atoi(argv[2]);
    filename = argv[1];
  }
  if (argc == 5){
    startStep = atoi (argv[4]);
  }
  printf ("# setting device to %d\n", atoi(argv[3]));
  cudaSetDevice (atoi(argv[3]));
  checkCUDAError ("set device");

  ScalorType dt = 0.005;
  ScalorType nlistExten = 0.5;
  ScalorType rebuildThreshold = 0.5 * nlistExten;
  ScalorType refT = 1.08088629683116564358;
  ScalorType tauT = .1;
  ScalorType refP = 5e-5;
  ScalorType tauP = 1.;
  ScalorType betaP = 23.;
  IndexType energy_feq = 10;
  IndexType config_feq = 10000;

  ScalorType headDiff = 0.1f;
  ScalorType head0d = 0.95f - 0.5 * headDiff;
  ScalorType head1d = 0.95f + 0.5 * headDiff;
  ScalorType taild = 1.0f;
  ScalorType strongw = 1.6f;  
  ScalorType weekw = 1.48f;
  
  MDSystem sys;
  sys.initConfig(filename);

  Topology::System sysTop;
  Topology::Molecule mol0;
  mol0.pushAtom (Topology::Atom (1.0, 0.0, 0));
  mol0.pushAtom (Topology::Atom (1.0, 0.0, 1));
  mol0.pushAtom (Topology::Atom (1.0, 0.0, 2));
  Topology::Molecule mol1;
  mol1.pushAtom (Topology::Atom (1.0, 0.0, 3));
  mol1.pushAtom (Topology::Atom (1.0, 0.0, 4));
  mol1.pushAtom (Topology::Atom (1.0, 0.0, 2));
  Topology::Molecule mol2;
  mol2.pushAtom (Topology::Atom (1.0, 0.0, 3));
  mol2.pushAtom (Topology::Atom (1.0, 0.0, 5));
  mol2.pushAtom (Topology::Atom (1.0, 0.0, 6));
  Topology::Molecule mol3;
  mol3.pushAtom (Topology::Atom (1.0, 0.0, 0));
  mol3.pushAtom (Topology::Atom (1.0, 0.0, 7));
  mol3.pushAtom (Topology::Atom (1.0, 0.0, 6));

  CosTailParameter cosparam_h0h0, cosparam_h0h1, cosparam_h1h1;
  CosTailParameter cosparam_h0t,  cosparam_h1t;
  CosTailParameter cosparam_Stail,  cosparam_Wtail, cosparam_Rtail;
  
  cosparam_h0h0.reinit (1.f, head0d, 0.f);
  cosparam_h1h1.reinit (1.f, head1d, 0.f);
  cosparam_h0h1.reinit (1.f, 0.5*(head0d+head1d), 0.f);
  cosparam_h0t .reinit (1.f, 0.5*(head0d+taild), 0.f);
  cosparam_h1t .reinit (1.f, 0.5*(head1d+taild), 0.f);
  cosparam_Stail.reinit (1.f, taild, strongw);
  cosparam_Wtail.reinit (1.f, taild, weekw);
  cosparam_Rtail.reinit (1.f, taild, 0.f);
  
  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(0, 0, cosparam_h0h0));
  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(3, 3, cosparam_h1h1));
  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(0, 3, cosparam_h0h1));
  
  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(0, 1, cosparam_h0t));
  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(0, 2, cosparam_h0t));
  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(0, 4, cosparam_h0t));
  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(0, 5, cosparam_h0t));
  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(0, 6, cosparam_h0t));
  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(0, 7, cosparam_h0t));

  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(3, 1, cosparam_h1t));
  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(3, 2, cosparam_h1t));
  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(3, 4, cosparam_h1t));
  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(3, 5, cosparam_h1t));
  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(3, 6, cosparam_h1t));
  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(3, 7, cosparam_h1t));

  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(1, 1, cosparam_Stail));
  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(2, 2, cosparam_Stail));
  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(4, 4, cosparam_Stail));
  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(5, 5, cosparam_Stail));
  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(6, 6, cosparam_Stail));
  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(7, 7, cosparam_Stail));

  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(1, 2, cosparam_Stail));
  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(2, 4, cosparam_Stail));
  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(5, 6, cosparam_Stail));
  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(6, 7, cosparam_Stail));

  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(1, 5, cosparam_Wtail));
  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(2, 6, cosparam_Wtail));
  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(2, 5, cosparam_Wtail));
  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(1, 6, cosparam_Wtail));
  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(4, 7, cosparam_Wtail));
  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(4, 6, cosparam_Wtail));
  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(2, 7, cosparam_Wtail));

  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(1, 4, cosparam_Rtail));
  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(1, 7, cosparam_Rtail));
  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(4, 5, cosparam_Rtail));
  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(5, 7, cosparam_Rtail));

  HarmonicSpringParameter hsparam;
  FENEParameter feneparam;
  hsparam.reinit (10.f, 4.f);
  feneparam.reinit (30.f, 1.5f);

  mol0.addBond (Topology::Bond (0, 1, feneparam));
  mol0.addBond (Topology::Bond (1, 2, feneparam));
  mol0.addBond (Topology::Bond (0, 2, hsparam));
  mol1.addBond (Topology::Bond (0, 1, feneparam));
  mol1.addBond (Topology::Bond (1, 2, feneparam));
  mol1.addBond (Topology::Bond (0, 2, hsparam));
  mol2.addBond (Topology::Bond (0, 1, feneparam));
  mol2.addBond (Topology::Bond (1, 2, feneparam));
  mol2.addBond (Topology::Bond (0, 2, hsparam));
  mol3.addBond (Topology::Bond (0, 1, feneparam));
  mol3.addBond (Topology::Bond (1, 2, feneparam));
  mol3.addBond (Topology::Bond (0, 2, hsparam));

  sysTop.addMolecules (mol0, sys.hdata.numAtom/12);
  sysTop.addMolecules (mol1, sys.hdata.numAtom/12);
  sysTop.addMolecules (mol2, sys.hdata.numAtom/12);
  sysTop.addMolecules (mol3, sys.hdata.numAtom/12);

  sys.initTopology (sysTop);
  sys.initDeviceData ();
  
  SystemBondedInteraction sysBdInter;
  sysBdInter.reinit (sysTop);

  BondedInteractionList bdInterList;
  bdInterList.reinit (sys, sysTop, sysBdInter);

  SystemNonBondedInteraction sysNbInter;
  sysNbInter.reinit (sysTop);
  
  ScalorType maxrcut = sysNbInter.maxRcut();
  ScalorType rlist = maxrcut + nlistExten;
  CellList clist (sys, rlist, NThreadsPerBlockCell, NThreadsPerBlockAtom,
		  RectangularBoxGeometry::mdRectBoxDirectionX |
		  RectangularBoxGeometry::mdRectBoxDirectionY);
  NeighborList nlist (sysNbInter, sys, rlist, NThreadsPerBlockAtom, 100.f);
  sys.normalizeDeviceData ();
  clist.rebuild (sys);
  nlist.rebuild (sys, clist);
  Displacement_max disp (sys, NThreadsPerBlockAtom);
  disp.recordCoord (sys);
  
  MDStatistic st (sys);
  MDStatistic last_st (sys);
  TranslationalFreedomRemover tfremover (sys, NThreadsPerBlockAtom);
  InteractionEngine inter (sys, NThreadsPerBlockAtom);
  inter.registNonBondedInteraction (sysNbInter);
  inter.registBondedInteraction    (sysBdInter);
  inter.clearInteraction (sys);
  inter.applyNonBondedInteraction (sys, nlist, st, NULL);
  inter.applyBondedInteraction    (sys, bdInterList, st);

  MDTimer timer;
  unsigned i;
  ScalorType seed = 1289028167;
  RandomGenerator_MT19937::init_genrand (seed);
  
  Thermostat_NoseHoover thermostat;
  thermostat.reinit (refT, dt, tauT, sys.ddata.numAtom * 3 - 3);
  Barostat_ParrinelloRahman barostat;
  barostat.reinit (dt, tauP, sys.box);
  barostat.assignGroup (mdRectBoxDirectionX, refP, betaP);
  LeapFrog_TPCouple_VCouple blpf (sys, NThreadsPerBlockAtom);
  blpf.addThermostat (thermostat);
  blpf.addBarostat   (barostat);

  Reshuffle resh (sys);
  
  timer.tic(mdTimeTotal);
  if (resh.calIndexTable (clist, &timer)){
    sys.reshuffle   (resh.indexTable, sys.hdata.numAtom, &timer);
    clist.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);  
    nlist.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);  
    disp.reshuffle  (resh.indexTable, sys.hdata.numAtom, &timer);
    bdInterList.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);
  }
  
  printf ("# prepare ok, start to run\n");
  printf ("#*     1     2           3        4         5       6           7           8           9        10    11    12    13\n");
  printf ("#* nstep  time  nonBondedE  bondedE  kineticE  totalE  pressurexx  pressureyy  pressurezz  pressure  boxx  boxy  boxz\n");
  // sys.recoverDeviceData (&timer);
  // sys.updateHostFromRecovered (&timer);
  // sys.writeHostDataGro ("confstart.gro", 0, 0.f, &timer);
  try{
    sys.initWriteXtc ("traj.xtc");
    sys.recoverDeviceData (&timer);
    sys.updateHostFromRecovered (&timer);
    sys.writeHostDataXtc (0, 0*dt, &timer);
    for (i = startStep; i < nstep + startStep; ++i){
      last_st.deviceCopy (st);
      if (i%10 == 0){
	tfremover.remove (sys, &timer);
      }
      st.clearDevice();
      blpf.oneStep (sys, dt, last_st, st, &timer);
      ScalorType maxdr = disp.calMaxDisplacemant (sys, &timer);
      if (maxdr > rebuildThreshold){
	// printf ("# Rebuild at step %09i ... ", i+1);
	// fflush(stdout);
	// rebuild
	sys.normalizeDeviceData (&timer);
	disp.recordCoord (sys);
	clist.rebuild (sys, &timer);
	nlist.rebuild (sys, clist, &timer);
	// printf ("done\n");
	// fflush(stdout);
      }
      inter.clearInteraction (sys);
      inter.applyNonBondedInteraction (sys, nlist,  st, NULL, &timer);
      inter.applyBondedInteraction    (sys, bdInterList, st, &timer);
      if ((i+1) % energy_feq == 0){
	st.updateHost();
	printf ("%09d %07e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.3f %.3f %.3f\n",
		(i+1),  
		(i+1) * dt, 
		st.nonBondedEnergy(),
		st.bondedEnergy(),
		st.kineticEnergy(),
		st.kineticEnergy() + st.nonBondedEnergy() + st.bondedEnergy(),
		st.pressureXX(sys.box),
		st.pressureYY(sys.box),
		st.pressureZZ(sys.box),
		st.pressure(sys.box),
		sys.box.size.x,
		sys.box.size.y,
		sys.box.size.z);
	fflush(stdout);
      }
      if ((i+1) % config_feq == 0){
      	sys.recoverDeviceData (&timer);
      	sys.updateHostFromRecovered (&timer);
      	sys.writeHostDataXtc (i+1, (i+1)*dt, &timer);
      }
      if ((i+1) % 100 == 0){
	if (resh.calIndexTable (clist, &timer)){
	  sys.reshuffle   (resh.indexTable, sys.hdata.numAtom, &timer);
	  clist.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);  
	  nlist.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);  
	  disp.reshuffle  (resh.indexTable, sys.hdata.numAtom, &timer);  
	  bdInterList.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);
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

  
