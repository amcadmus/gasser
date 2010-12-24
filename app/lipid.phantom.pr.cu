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

#define NThreadsPerBlockCell	64
#define NThreadsPerBlockAtom	64

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

  ScalorType dt = 0.001;
  ScalorType nlistExten = 0.5;
  ScalorType rebuildThreshold = 0.5 * nlistExten;
  ScalorType refT = 1.08088629683116564358;
  ScalorType tauT = .1;
  ScalorType refP = 5e-5;
  ScalorType tauP = 1.;
  ScalorType betaP = 23.;
  IndexType energy_feq = 100;
  IndexType config_feq = 10000;
  
  MDSystem sys;
  sys.initConfig(filename);

  IndexType idx = 1;
  for (; idx < sys.hdata.numAtom; ++idx){
    if (strcmp ("HD", &sys.hdata.resdName[idx*StringSize]) == 0){
      break;
    }
  }
  IndexType ntail = idx - 1;
  for (; idx < sys.hdata.numAtom; ++idx){
    if (strcmp ("Sol", &sys.hdata.resdName[idx*StringSize]) == 0){
      break;
    }
  }
  IndexType nLipidAtom = idx;
  IndexType nLipidMol = nLipidAtom / (1 + ntail);
  IndexType nSolAtom = sys.hdata.numAtom - nLipidAtom;
  
  Topology::System sysTop;
  Topology::Molecule molLipid;
  Topology::Molecule molSol;
  molLipid.pushAtom (Topology::Atom (1.0, 0.0, 0));
  for (IndexType i = 0; i < ntail; ++i){
    molLipid.pushAtom (Topology::Atom (1.0, 0.0, 1));
  }
  molSol.pushAtom (Topology::Atom (1.0, 0.0, 2));  

  ScalorType sigma_tt = 1.f;
  ScalorType sigma_hh = 1.1f;
  ScalorType sigma_ht = 0.5 * (sigma_hh + sigma_tt);
  ScalorType rcut_tt  = 2.f;
  ScalorType rcut_ht  = sigma_ht;
  ScalorType rcut_hh  = sigma_hh;

  ScalorType sigma_s  = sigma_hh;
  ScalorType sigma_hs = 0.5 * (sigma_hh + sigma_s);
  ScalorType sigma_ts = 0.5 * (sigma_tt + sigma_s);
  ScalorType rcut_hs = sigma_hs;
  ScalorType rcut_ts = sigma_ts;
  ScalorType capValue = 100.f;
  
  LennardJones6_12BCapParameter p_tt_tmp, p_tt;
  p_tt_tmp.reinit (1.f, 1.f, 0.f, rcut_tt, capValue);
  p_tt.reinit (1.f, 1.f, -p_tt_tmp.shiftAtCut(), rcut_tt, capValue);
  LennardJones6_12BCapParameter p_ht_tmp, p_ht;
  p_ht_tmp.reinit (1.f, 1.f, 0.f, rcut_ht, capValue);
  p_ht.reinit (1.f, 1.f, -p_ht_tmp.shiftAtCut(), rcut_ht, capValue);
  LennardJones6_12BCapParameter p_hh_tmp, p_hh;
  p_hh_tmp.reinit (1.f, 1.f, 0.f, rcut_hh, capValue);
  p_hh.reinit (1.f, 1.f, -p_hh_tmp.shiftAtCut(), rcut_hh, capValue);

  LennardJones6_12BCapParameter p_hs_tmp, p_hs;
  p_hs_tmp.reinit (1.f, 1.f, 0.f, rcut_hs, capValue);
  p_hs.reinit (1.f, 1.f, -p_hs_tmp.shiftAtCut(), rcut_hs, capValue);
  LennardJones6_12BCapParameter p_ts_tmp, p_ts;
  p_ts_tmp.reinit (1.f, 1.f, 0.f, rcut_ts, capValue);
  p_ts.reinit (1.f, 1.f, -p_ts_tmp.shiftAtCut(), rcut_ts, capValue);

  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(0, 0, p_hh));
  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(0, 1, p_ht));
  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(1, 1, p_tt));
  
  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(2, 0, p_hs));
  sysTop.addNonBondedInteraction(Topology::NonBondedInteraction(2, 1, p_ts));  

  FENE2Parameter p_fene;
  p_fene.reinit (100.f, 0.2f * sigma_tt, 0.7f * sigma_tt);

  for (IndexType i = 0; i < ntail; ++i){
    molLipid.addBond (Topology::Bond (i, i+1, p_fene));
    molLipid.addExclusion (Topology::Exclusion (i, i+1));
  }
  
  sysTop.addMolecules (molLipid, nLipidMol);
  sysTop.addMolecules (molSol, nSolAtom);

  sys.initTopology (sysTop);
  sys.initDeviceData ();
  
  SystemBondedInteraction sysBdInter;
  sysBdInter.reinit (sysTop);

  BondedInteractionList bdInterList;
  bdInterList.reinit (sys, sysTop, sysBdInter);

  SystemNonBondedInteraction sysNbInter;
  sysNbInter.reinit (sysTop);

  ExclusionList exclList;
  exclList.reinit (sys, sysTop, sysNbInter);
  
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
  inter.applyNonBondedInteraction (sys, nlist, st, &exclList);
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
  // blpf.addThermostat (thermostat);
  // blpf.addBarostat   (barostat);

  Reshuffle resh (sys);
  
  timer.tic(mdTimeTotal);
  if (resh.calIndexTable (clist, &timer)){
    sys.reshuffle   (resh.indexTable, sys.hdata.numAtom, &timer);
    clist.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);  
    nlist.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);  
    disp.reshuffle  (resh.indexTable, sys.hdata.numAtom, &timer);
    bdInterList.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);
    exclList.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);
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
      inter.applyNonBondedInteraction (sys, nlist,  st, &exclList, &timer);
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
      if ((i+1) % 50 == 0){
      	if (resh.calIndexTable (clist, &timer)){
      	  sys.reshuffle   (resh.indexTable, sys.hdata.numAtom, &timer);
      	  clist.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);  
      	  nlist.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);  
      	  disp.reshuffle  (resh.indexTable, sys.hdata.numAtom, &timer);  
      	  bdInterList.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);
      	  exclList.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);
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

  
