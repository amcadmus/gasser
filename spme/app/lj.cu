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
#include "EwaldSumRec.h"
#include "SPMERec.h"

#define NThreadsPerBlockCell	128
#define NThreadsPerBlockAtom	96

void printCoord (const MDSystem & sys)
{
  FILE * fp = fopen ("conf.save", "w");
  fprintf (fp, "%d\n", sys.hdata.numAtom);
  for (IndexType i = 0; i < sys.hdata.numAtom; ++i){
    fprintf (fp, "%e %e %e   %e\n",
	     sys.hdata.coord[i].x,
	     sys.hdata.coord[i].y,
	     sys.hdata.coord[i].z,
	     sys.hdata.charge[i]);
  }
  fprintf (fp, "%f\n", sys.box.size.x);
  fclose (fp);
}


int main(int argc, char * argv[])
{
  IndexType nstep = 20;
  char * filename;
  ScalorType rcut = 3.2f;
  ScalorType nlistExten = 0.3;
  ScalorType nlistExtenFactor = 10.f;
  ScalorType dt = 0.001;
  ScalorType beta = 1.3;
  IndexType order = 4;
  int Kvalue = 8;
  
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
  Topology::Molecule mol0;
  mol0.pushAtom (Topology::Atom (1.0,-1.0, 0));
  Topology::Molecule mol1;
  mol1.pushAtom (Topology::Atom (1.0, 1.0, 0));
  LennardJones6_12Parameter ljparam;
  ljparam.reinit (1.f, 1.f, 0.f, rcut);
  sysTop.addNonBondedInteraction (Topology::NonBondedInteraction(0, 0, ljparam));
  sysTop.addMolecules (mol0, sys.hdata.numAtom/2);
  sysTop.addMolecules (mol1, sys.hdata.numAtom - sys.hdata.numAtom/2);

  sys.initTopology (sysTop);
  sys.initDeviceData ();

  printCoord (sys);

  MDStatistic myst(sys);
  MatrixType vecA;
  vecA.xy = vecA.yx = vecA.yz = vecA.zy = vecA.xz = vecA.zx = 0.f;
  vecA.xx = sys.box.size.x;
  vecA.yy = sys.box.size.y;
  vecA.zz = sys.box.size.z;
  IntVectorType K;
  K.x = K.y = K.z = Kvalue;
  K.x += 2;
  K.y += 1;

  // EwaldSumRec ewald;
  // ewald.reinit (vecA, K, beta, sys.hdata.numAtom, 7, 13);
  // ewald.applyInteraction (sys, &myst);
  // sys.updateHostFromDevice ();
  // printf ("force: %f %f %f\n",
  // 	  sys.hdata.forcx[0],
  // 	  sys.hdata.forcy[0],
  // 	  sys.hdata.forcz[0]);
  // myst.updateHost ();
  // printf ("energy: %e, pressure: %e %e %e\n",
  // 	  myst.nonBondedEnergy(),
  // 	  myst.pressureXX(sys.box),
  // 	  myst.pressureYY(sys.box),
  // 	  myst.pressureZZ(sys.box));
  
  SystemNonBondedInteraction sysNbInter;
  sysNbInter.reinit (sysTop);
  
  ScalorType maxrcut = sysNbInter.maxRcut();
  ScalorType rlist = maxrcut + nlistExten;
  CellList clist (sys, rlist, NThreadsPerBlockCell, NThreadsPerBlockAtom);
  NeighborList nlist (sysNbInter,
		      sys,
		      rlist,
		      NThreadsPerBlockAtom,
		      nlistExtenFactor);
  sys.normalizeDeviceData ();
  clist.rebuild (sys, NULL);
  nlist.rebuild (sys, clist, NULL);
  Displacement_max disp (sys, NThreadsPerBlockAtom);
  disp.recordCoord (sys);
  
  MDStatistic st(sys);
  VelocityVerlet inte_vv (sys, NThreadsPerBlockAtom);
  ScalorType refT = 1.;
  VelocityRescale inte_vr (sys, NThreadsPerBlockAtom, refT, 0.1);
  TranslationalFreedomRemover tfremover (sys, NThreadsPerBlockAtom);
  InteractionEngine inter (sys, NThreadsPerBlockAtom);
  inter.registNonBondedInteraction (sysNbInter);
  
  MDTimer timer;
  unsigned i;
  ScalorType seed = 1;
  RandomGenerator_MT19937::init_genrand (seed);

  Reshuffle resh (sys);

  SPMERecIk spme;
  spme.reinit (vecA, K, order, beta, sys.hdata.numAtom, 7, 13);
  spme.applyInteraction (sys, NULL, NULL);
  FILE * fp = fopen ("myForce", "w");
  for (IndexType i = 0; i < sys.ddata.numAtom; ++i){
    fprintf (fp, "%.16e\t%.16e\t%.16e\n",
	     sys.ddata.forcx[i],
	     sys.ddata.forcy[i],
	     sys.ddata.forcz[i]);
  }
  fclose (fp);
  timer.tic(mdTimeTotal);
  if (resh.calIndexTable (clist, &timer)){
    sys.reshuffle   (resh.indexTable, sys.hdata.numAtom, &timer);
    clist.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);  
    nlist.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);  
    disp.reshuffle  (resh.indexTable, sys.hdata.numAtom, &timer);  
  }
  
  printf ("# prepare ok, start to run\n");
  // sys.recoverDeviceData (&timer);
  // sys.updateHostFromRecovered (&timer);
  // sys.writeHostDataGro ("confstart.gro", 0, 0.f, &timer);
  printf ("# prepare ok, start to run\n");
  printf ("#*     1     2           3         4            5       6          7-9\n");
  printf ("#* nstep  time  nonBondedE  kineticE  temperature  totalE  pressurexyz\n");
  try{
    // sys.initWriteXtc ("traj.xtc");
    // sys.recoverDeviceData (&timer);
    // sys.updateHostFromRecovered (&timer);
    // sys.writeHostDataXtc (0, 0*dt, &timer);
    for (i = 0; i < nstep; ++i){
      // if (i%10 == 0){
      // 	tfremover.remove (sys, &timer);
      // }
      // if ((i+1) % 100 == 0){
      // 	st.clearDevice();
      // 	inte_vv.step1 (sys, dt, &timer);
      // 	inter.clearInteraction (sys);
      // 	ScalorType maxdr = disp.calMaxDisplacemant (sys, &timer);
      // 	if (maxdr > nlistExten * 0.5){
      // 	  // printf ("# Rebuild at step %09i ... ", i+1);
      // 	  // fflush(stdout);
      // 	  // rebuild
      // 	  sys.normalizeDeviceData (&timer);
      // 	  disp.recordCoord (sys);
      // 	  clist.rebuild (sys, &timer);
      // 	  inter.applyNonBondedInteraction (sys, clist, rcut, st, &timer);
      // 	  nlist.rebuild (sys, clist, &timer);
      // 	  // printf ("done\n");
      // 	  // fflush(stdout);
      // 	}
      // 	else{
      // 	  inter.applyNonBondedInteraction (sys, nlist, st, NULL, &timer);
      // 	  // inter.applyNonBondedInteraction (sys, rcut, st, &timer);
      // 	}
      // 	inte_vv.step2 (sys, dt, st, &timer);
      // 	st.updateHost();
      // 	printf ("%09d %07e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e \n",
      // 		(i+1),  
      // 		(i+1) * dt, 
      // 		st.nonBondedEnergy(),
      // 		st.kineticEnergy(),
      // 		st.kineticEnergy() / (sys.ddata.numAtom - 1) * 2./3.,
      // 		st.nonBondedEnergy() + st.bondedEnergy() + st.kineticEnergy(),
      // 		st.pressureXX(sys.box),
      // 		st.pressureYY(sys.box),
      // 		st.pressureZZ(sys.box),
      // 		st.pressure(sys.box));
      // 	fflush(stdout);
      // }
      // else {
      // 	inte_vv.step1 (sys, dt, &timer);
      // 	inter.clearInteraction (sys);
      // 	ScalorType maxdr = disp.calMaxDisplacemant (sys, &timer);
      // 	if (maxdr > nlistExten * 0.5){
      // 	  // printf ("# Rebuild at step %09i ... ", i+1);
      // 	  // fflush(stdout);
      // 	  // rebuild
      // 	  sys.normalizeDeviceData (&timer);
      // 	  disp.recordCoord (sys);
      // 	  clist.rebuild (sys, &timer);
      // 	  inter.applyNonBondedInteraction (sys, clist, rcut, &timer);
      // 	  nlist.rebuild (sys, clist, &timer);
      // 	  // printf ("done\n");
      // 	  // fflush(stdout);
      // 	}
      // 	else{
      // 	  inter.applyNonBondedInteraction (sys, nlist, NULL, &timer);
      // 	  // inter.applyNonBondedInteraction (sys, rcut, &timer);
      // 	}
      // 	inte_vv.step2 (sys, dt, &timer);
      // }
      // // if (maxdr > nlistExten * 0.5){
      // // 	printf ("# Rebuild at step %09i ... ", i+1);
      // // 	fflush(stdout);
      // // 	// rebuild
      // // 	sys.normalizeDeviceData ();
      // // 	disp.recordCoord (sys);
      // // 	clist.rebuild (sys, &timer);
      // // 	nlist.rebuild (sys, clist, &timer);
      // // 	printf ("done\n");
      // // 	fflush(stdout);
      // // }
      // // if ((i+1) % 1000 == 0){
      // // 	sys.recoverDeviceData (&timer);
      // // 	sys.updateHostFromRecovered (&timer);
      // // 	sys.writeHostDataXtc (i+1, (i+1)*dt, &timer);
      // // }
      // if ((i+1) % 100 == 0){
      // 	if (resh.calIndexTable (clist, &timer)){
      // 	  sys.reshuffle   (resh.indexTable, sys.hdata.numAtom, &timer);
      // 	  clist.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);  
      // 	  nlist.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);  
      // 	  disp.reshuffle  (resh.indexTable, sys.hdata.numAtom, &timer);  
      // 	}
      // }
    }
    // sys.endWriteXtc();
    // sys.recoverDeviceData (&timer);
    // sys.updateHostFromRecovered (&timer);
    // sys.writeHostDataGro ("confout.gro", nstep, nstep*dt, &timer);
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

  
