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
  ScalorType lattice_k = 1000.f;
  char * filename;
  IndexType numAtom_A = 1000;
  IndexType numAtom_B = 1000;
  IndexType numAtom_a = 1000;
  IndexType numAtom_b = 1000;
  IndexType numAtom_c = 1000;
  
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
  if (sys.hdata.numAtom !=
      numAtom_a + numAtom_b + numAtom_c + numAtom_A + numAtom_B) {
    printf ("# inconsistent number of atom!\n");
    exit (1);
  }
  for (unsigned i = 0; i < sys.hdata.numAtom ; ++i){
    if (strcmp(&(sys.hdata.atomName[i*StringSize]), "lja") == 0) {
      continue;
    }
    else {
      numAtom_a = i;
      break;
    }
  }
  for (unsigned i = numAtom_a; i < sys.hdata.numAtom ; ++i){
    if (strcmp(&(sys.hdata.atomName[i*StringSize]), "ljb") == 0) {
      continue;
    }
    else {
      numAtom_b = i;
      break;
    }
  }
  numAtom_b -= numAtom_a;
  for (unsigned i = numAtom_b; i < sys.hdata.numAtom ; ++i){
    if (strcmp(&(sys.hdata.atomName[i*StringSize]), "ljc") == 0) {
      continue;
    }
    else {
      numAtom_c = i;
      break;
    }
  }
  numAtom_c -= numAtom_b;
  for (unsigned i = numAtom_c; i < sys.hdata.numAtom ; ++i){
    if (strcmp(&(sys.hdata.atomName[i*StringSize]), "ljA") == 0) {
      continue;
    }
    else {
      numAtom_A = i;
      break;
    }
  }
  numAtom_A -= numAtom_c;
  for (unsigned i = numAtom_A; i < sys.hdata.numAtom ; ++i){
    if (strcmp(&(sys.hdata.atomName[i*StringSize]), "ljB") == 0) {
      continue;
    }
    else {
      numAtom_B = i;
      break;
    }
  }
  numAtom_B -= numAtom_A;
  printf ("# numAtoms are: %d %d %d %d %d\n",
	  numAtom_a, numAtom_b, numAtom_c,
	  numAtom_A, numAtom_B);
  
  LennardJones6_12Parameter ljparam_attractive;
  LennardJones6_12Parameter ljparam_repulsive;
  LennardJones6_12Parameter ljparam_attractive0p8;
  ljparam_attractive.reinit    (1.f, 1.f, 1.f, 0.f, rcut);
  ljparam_attractive0p8.reinit (1.f, 1.f, .8f, 0.f, rcut);
  ljparam_repulsive .reinit (1.f, 1.f,-1.f, 0.f, rcut);

  Topology::System sysTop;
  Topology::Molecule mol_a;
  mol_a.pushAtom (Topology::Atom (1.0, 0.0, 0));
  Topology::Molecule mol_b;
  mol_b.pushAtom (Topology::Atom (1.0, 0.0, 1));
  Topology::Molecule mol_c;
  mol_c.pushAtom (Topology::Atom (1.0, 0.0, 4));
  Topology::Molecule mol_A;
  mol_A.pushAtom (Topology::Atom (1.0, 0.0, 2));
  Topology::Molecule mol_B;
  mol_B.pushAtom (Topology::Atom (1.0, 0.0, 3));

  sysTop.addMolecules (mol_a, numAtom_a);
  sysTop.addMolecules (mol_b, numAtom_b);
  sysTop.addMolecules (mol_c, numAtom_c);
  sysTop.addMolecules (mol_A, numAtom_A);
  sysTop.addMolecules (mol_B, numAtom_B);

  sysTop.addNonBondedInteraction (Topology::NonBondedInteraction(2, 2, ljparam_attractive));
  sysTop.addNonBondedInteraction (Topology::NonBondedInteraction(3, 3, ljparam_attractive));
  sysTop.addNonBondedInteraction (Topology::NonBondedInteraction(2, 3, ljparam_repulsive));

  sysTop.addNonBondedInteraction (Topology::NonBondedInteraction(4, 2, ljparam_attractive0p8));
  sysTop.addNonBondedInteraction (Topology::NonBondedInteraction(4, 3, ljparam_attractive0p8));

  sysTop.addNonBondedInteraction (Topology::NonBondedInteraction(0, 2, ljparam_attractive));
  sysTop.addNonBondedInteraction (Topology::NonBondedInteraction(1, 3, ljparam_attractive));
  sysTop.addNonBondedInteraction (Topology::NonBondedInteraction(0, 3, ljparam_repulsive));
  sysTop.addNonBondedInteraction (Topology::NonBondedInteraction(1, 2, ljparam_repulsive));  

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
  NoseHoover_Chains2 nhc;
  nhc.reinit (sys, NThreadsPerBlockAtom, refT, tauT);

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
  printf ("#*     1     2           3         4            5       6         8   9\n");
  printf ("#* nstep  time  nonBondedE  kineticE  temperature  totalE  pressure box\n");

  try{
    sys.initWriteXtc ("traj.xtc");
    sys.recoverDeviceData (&timer);
    sys.updateHostFromRecovered (&timer);
    sys.writeHostDataXtc (0, 0*dt, &timer);
    for (i = 0; i < nstep; ++i){
      if (i%10 == 0){
	tfremover.remove (sys, &timer);
      }
      
      nhc.operator_L (0.5 * dt, sys, &timer);
      inte_vv.step1 (sys, dt, &timer);

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
	nlist.rebuild (sys, clist, &timer);
	// printf ("done\n");
	// fflush(stdout);
      }
      inter.applyNonBondedInteraction (sys, nlist, st, NULL, &timer);
      inter.applyLatticeInteraction (sys, lattice_k, 0, 1, &timer);

      inte_vv.step2 (sys, dt, &timer);
      if ((i+1) % thermoFeq == 0){	
	nhc.operator_L (0.5 * dt, sys, st, &timer);
      }
      else {
	nhc.operator_L (0.5 * dt, sys, &timer);	
      }      

      if ((i+1) % thermoFeq == 0){
	st.updateHost ();
	printf ("%09d %07e %.7e %.7e %.7e %.7e %.7e %.7e\n",
		(i+1),  
		(i+1) * dt, 
		st.nonBondedEnergy(),
		st.kineticEnergy(),
		st.kineticEnergy() * 2. / 3. / (double (sys.hdata.numAtom) - 3.),
		st.nonBondedEnergy() +
		st.kineticEnergy(),
		st.pressure(sys.box),
		sys.box.size.x
	    );
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

  
