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


// #define NThreadsPerBlockCell	32
// #define NThreadsPerBlockAtom	4

#define NThreadsPerBlockCell	160
#define NThreadsPerBlockAtom	96

int main(int argc, char * argv[])
{
  IndexType nstep = 100000;
  IndexType confFeq = 100000;
  IndexType thermoFeq = 100;
  ScalorType rcut = 4.0;
  ScalorType nlistExten = 0.3;
  ScalorType refT = 1.50;
  ScalorType tauT = 1.;
  ScalorType refP = 0.2;
  ScalorType tauP = 1.;
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
  sysTop.addNonBondedInteraction (Topology::NonBondedInteraction(0, 0, ljparam));
  sysTop.addMolecules (mol, sys.hdata.numAtom);

  sys.initTopology (sysTop);
  sys.initDeviceData ();
  
  SystemNonBondedInteraction sysNbInter;
  sysNbInter.reinit (sysTop);
  ScalorType energyCorr = sysNbInter.energyCorrection ();
  ScalorType pressureCorr = sysNbInter.pressureCorrection ();
  
  ScalorType maxrcut = sysNbInter.maxRcut();
  ScalorType rlist = maxrcut + nlistExten;
  CellList clist (sys, rlist, NThreadsPerBlockCell, NThreadsPerBlockAtom);
  NeighborList nlist (sysNbInter, sys, rlist, NThreadsPerBlockAtom, 10.f);
  sys.normalizeDeviceData ();
  clist.rebuild (sys, NULL);
  nlist.rebuild (sys, clist, NULL);
  Displacement_max disp (sys, NThreadsPerBlockAtom);
  disp.recordCoord (sys);
  
  MDStatistic st(sys);
  MDStatistic old_st(sys), tmp_st(sys);
  old_st.setEnergyCorr (energyCorr);
  old_st.setPressureCorr (pressureCorr);
  tmp_st.setEnergyCorr (energyCorr);
  tmp_st.setPressureCorr (pressureCorr);
  TranslationalFreedomRemover tfremover (sys, NThreadsPerBlockAtom);
  InteractionEngine inter (sys, NThreadsPerBlockAtom);
  inter.registNonBondedInteraction (sysNbInter);
  
  MDTimer timer;
  unsigned i;
  ScalorType dt = 0.002;
  ScalorType seed = 1;
  RandomGenerator_MT19937::init_genrand (seed);

  VelocityVerlet inte_vv (sys, NThreadsPerBlockAtom);
  VelocityRescale inte_vr (sys, NThreadsPerBlockAtom, refT, 0.1);
  NoseHoover_Chains2_Isobaric nhcp;
  nhcp.reinit (sys, NThreadsPerBlockAtom, refT, tauT, refP, tauP, pressureCorr);

  Reshuffle resh (sys);
  
  timer.tic(mdTimeTotal);
  if (resh.calIndexTable (clist, &timer)){
    sys.reshuffle   (resh.indexTable, sys.hdata.numAtom, &timer);
    clist.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);  
    nlist.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);  
    disp.reshuffle  (resh.indexTable, sys.hdata.numAtom, &timer);  
  }
  
  printf ("# prepare ok, start to run\n");
  sys.recoverDeviceData (&timer);
  sys.updateHostFromRecovered (&timer);
  sys.writeHostDataGro ("confstart.gro", 0, 0.f, &timer);
  printf ("# prepare ok, start to run\n");
  printf ("#*     1     2           3         4            5       6         7    8    9     10  11     12   13   14    15   16    17\n");
  printf ("#* nstep  time  nonBondedE  kineticE  temperature  totalE  pressure  box  volume  h0  h1  NHC_H  vep  xi1  vxi1  xi2  vxi2\n");

  try{
    sys.initWriteXtc ("traj.xtc");
    sys.recoverDeviceData (&timer);
    sys.updateHostFromRecovered (&timer);
    sys.writeHostDataXtc (0, 0*dt, &timer);
    for (i = 0; i < nstep; ++i){
      if (i%10 == 0){
	tfremover.remove (sys, &timer);
      }
      old_st.deviceCopy (st);
      old_st.updateHost ();
      nhcp.operator_L_CP (0.5 * dt, sys, old_st, &timer);
      
      nhcp.operator_L_v  (0.5 * dt, sys, &timer);

      nhcp.operator_L_r  (dt, sys, &timer);
      nhcp.operator_L_box (dt, sys.box);
      
      st.clearDevice();
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
      }
      
      tmp_st.clearDevice();
      nhcp.operator_L_v  (0.5 * dt, sys, tmp_st, &timer);
      mdStatisticItem_t tmp_array[3];
      tmp_array[0] = mdStatisticVirialXX;
      tmp_array[1] = mdStatisticVirialYY;
      tmp_array[2] = mdStatisticVirialZZ;
      tmp_st.add (st, 3, tmp_array);
      
      nhcp.operator_L_CP (0.5 * dt, sys, tmp_st, st, &timer);

      if ((i+1) % thermoFeq == 0){
	st.updateHost ();
	ScalorType e = st.nonBondedEnergy () + st.kineticEnergy();	
	ScalorType v = sys.box.size.x * sys.box.size.y * sys.box.size.z;
	ScalorType p0 = st.pressure(sys.box);
	ScalorType h0 = v * p0 + e;
	ScalorType h1 = v * refP + e;

	printf ("%09d %07e %.7e %.7e %.7e %.7e %.7e %.7f %.7e %.7e %.7e %.5e %.5e %.5e %.5e\n",
		(i+1),  
		(i+1) * dt, 
		st.nonBondedEnergy (),
		st.kineticEnergy(),
		st.kineticEnergy() * 2. / (3. * (double (sys.hdata.numAtom) - 1.)),
		e,
		p0,
		sys.box.size.x,
		v,
		h0,
		h1,
		e + nhcp.HamiltonianContribution (sys.box),
		nhcp.vep,
		nhcp.xi1,
		nhcp.vxi1);
	fflush(stdout);
      }

      if ((i+1) % confFeq == 0){
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

  
