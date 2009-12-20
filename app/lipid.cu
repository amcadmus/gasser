/**
 * @file   lipid.cu
 * @author Han Wang <han_wang@math.pku.edu.cn>
 * @date   Thu Nov 19 12:53:01 2009
 * 
 * @brief  the main program to test the membrane simulation.
 * 
 * 
 */
#include <stdio.h>
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
#include "BondList_interface.h"

#include "MDSystem_interface.h"


#define NThreadsPerBlockCell	192
#define NThreadsPerBlockAtom	128

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
  sys.initConfig(filename, "lipid.map");
  sys.initNBForce(2);
  ScalorType cosparam[mdForceNParamCosTail];  
  CosTail::initParameter (cosparam, 1.f, 0.95f, 0.f);
  sys.addNBForce (0, 0, mdForceCosTail, cosparam);
  CosTail::initParameter (cosparam, 1.f, 0.975f, 0.f);
  sys.addNBForce (0, 1, mdForceCosTail, cosparam);
  CosTail::initParameter (cosparam, 1.f, 1.0f, 1.6f);
  sys.addNBForce (1, 1, mdForceCosTail, cosparam);

  sys.initBond ();
  ScalorType hsparam[mdForceNParamHarmonicSpring] ;
  HarmonicSpring::initParameter (hsparam, 10.f, 4.f);
  ScalorType feneparam[mdForceNParamFENE];
  FENE::initParameter (feneparam, 30.f, 1.5f);
  
  for (unsigned i = 0; i < sys.hdata.numAtom; i+=3){
    sys.addBond (i, i+1, mdForceFENE, feneparam);
    sys.addBond (i+2, i+1, mdForceFENE, feneparam);
    sys.addBond (i, i+2, mdForceHarmonicSpring, hsparam);
  }   
  sys.buildBond();

  ScalorType tmpsum = 0.;
  for (IndexType i = 0; i < sys.hdata.numAtom; ++i){
    tmpsum += 0.5 * (sys.hdata.velox[i] * sys.hdata.velox[i] +
		     sys.hdata.veloy[i] * sys.hdata.veloy[i] +
		     sys.hdata.veloz[i] * sys.hdata.veloz[i] );
  }
  printf ("# tmpsum is %f\n", tmpsum);

  ScalorType maxrcut = sys.calMaxNBRcut ();
  printf ("# max rcut is %f\n", maxrcut);
  ScalorType nlistExten = 0.5f;
  ScalorType rlist = maxrcut + nlistExten;
  // NeighborList nlist(sys, rlist, NThreadsPerBlockCell, 20,
  // 		     RectangularBoxGeometry::mdRectBoxDirectionX |
  // 		     RectangularBoxGeometry::mdRectBoxDirectionY |
  // 		     RectangularBoxGeometry::mdRectBoxDirectionZ);;
  NeighborList nlist(sys, rlist, NThreadsPerBlockCell, 40,
  		     RectangularBoxGeometry::mdRectBoxDirectionX |
  		     RectangularBoxGeometry::mdRectBoxDirectionY);
  
  nlist.build (sys);
  
  Reshuffle resh (sys, nlist, NThreadsPerBlockCell);
  resh.shuffleSystem ( sys, nlist);

  MDStatistic st(sys);

  VelocityVerlet inte (sys, NThreadsPerBlockAtom);;

  ScalorType refT = 0.9977411970749;
  VelocityRescale inte_vr (sys, NThreadsPerBlockAtom, refT, 0.1);
// // printf ("%f %f\n", ddata.velox[0], ddata.velox[1]);
  // inte.removeTranslationalFreedom (ddata);
  // // printf ("%f %f\n", ddata.velox[0], ddata.velox[1]);

  TranslationalFreedomRemover tfremover (sys, NThreadsPerBlockAtom);

  InteractionEngine_interface interaction(sys, NThreadsPerBlockAtom);;

  MDTimer timer;
  unsigned i;
  ScalorType dt = 0.005;
  ScalorType seed = 1;
  RandomGenerator_MT19937::init_genrand (seed);
  
  printf ("# prepare ok, start to run\n");

  try{
    timer.tic(mdTimeTotal);
    sys.initWriteXtc ("traj.xtc");
    sys.writeHostDataXtc (0, 0*dt, &timer);
    for (i = 0; i < nstep; ++i){ 
      if (i%10 == 0){
	tfremover.remove (sys, &timer);
      }
      // if (i%1 == 0){
      //   if (i == 0) nlist.build (sys);
      //   else nlist.reBuild (sys);
      //   // resh.shuffleSystem (sys, nlist);
      // }   
	  
      if ((i+1) % 10 == 0){
	st.clearDevice();
	inte_vr.step1 (sys, dt, &timer);
	interaction.applyInteraction (sys, nlist, st, &timer);
	inte_vr.step2 (sys, dt, st, &timer);
	st.updateHost();
	printf ("%09d %07e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e\n",
		(i+1),  
		(i+1) * dt, 
		st.getStatistic(mdStatisticNonBondedPotential),
		st.getStatistic(mdStatisticBondedPotential),
		st.kineticEnergy(),
		st.getStatistic(mdStatisticNonBondedPotential) +
		st.getStatistic(mdStatisticBondedPotential) +
		st.kineticEnergy(),
		st.getStatistic(mdStatisticVirialXX)*0.5,
		st.getStatistic(mdStatisticVirialYY)*0.5,
		st.getStatistic(mdStatisticVirialZZ)*0.5,
		st.pressureXX(),
		st.pressureYY(),
		st.pressureZZ(),
		st.pressure());
	fflush(stdout);
      }
      else {
	inte_vr.step1 (sys, dt, &timer);
	interaction.applyInteraction (sys, nlist, &timer);
	inte_vr.step2 (sys, dt, &timer);
      }
      if (nlist.judgeRebuild(sys, 0.5 * nlistExten, &timer)){
	printf ("# Rebuild at step %09i ... ", i+1);
	fflush(stdout);
	nlist.reBuild(sys, &timer);
	printf ("done\n");
	fflush(stdout);
      }
      if ((i+1) % 1000 == 0){
	// sys.updateHost(&timer);
	resh.recoverMDDataToHost (sys, &timer);
	sys.writeHostDataXtc (i+1, (i+1)*dt, &timer);
      }
      if ((i+1) % 200 == 0){
	resh.shuffleSystem (sys, nlist, &timer);
      }
    }
    sys.endWriteXtc();
    timer.toc(mdTimeTotal);
    timer.printRecord (stderr);
  }
  catch (MDExcptCuda & e){
    resh.recoverMDDataToHost (sys, &timer);
    sys.writeHostDataXtc (i+1, (i+1)*dt, &timer);
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
