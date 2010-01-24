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

//   ScalorType cap = 100;
//   ScalorType cosparam[mdForceNParamCosTail_cap];  
//   CosTail_cap::initParameter (cosparam, 1.f, 0.95f, 0.f, cap);
//   sys.addNBForce (0, 0, mdForceCosTail_cap, cosparam);
//   CosTail_cap::initParameter (cosparam, 1.f, 0.975f, 0.f, cap);
//   sys.addNBForce (0, 1, mdForceCosTail_cap, cosparam);
//   CosTail_cap::initParameter (cosparam, 1.f, 1.0f, 1.6f, cap);
//   sys.addNBForce (1, 1, mdForceCosTail_cap, cosparam);
  
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
  sys.initAngle();
  sys.buildAngle();
  
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

  sys.writeHostDataGro ("confstart.gro", 0, 0.f, &timer);
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
    resh.recoverMDDataToHost (sys, &timer);
    sys.writeHostDataGro ("confout.gro", nstep, nstep*dt, &timer);
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


// typedef enum paramIndex {
//   epsilon		= 0,
//   bv		= 1,
//   rc		= 2,
//   wc		= 3,
//   wci		= 4,
//   rcut		= 5,
//   cap		= 6,
//   capR		= 7,
//   pcapR		= 8
// } paramIndex_t;




// void
// force (const ScalorType * param,
// 		    ScalorType diffx,
// 		    ScalorType diffy,
// 		    ScalorType diffz,
// 		    ScalorType *fx, 
// 		    ScalorType *fy, 
// 		    ScalorType *fz)
// {
//   ScalorType dr2 = diffx*diffx + diffy*diffy + diffz*diffz;
//   ScalorType fscalor = 0.f;
//   ScalorType ri = 1.f / sqrtf (dr2);

//   if (dr2 == 0.f){
//     *fx = - param[cap];
//     *fy = 0.f;
//     *fz = 0.f;
//     return;
//   }
//   else if (dr2 < param[capR]*param[capR]){
//     ScalorType s = 1.f / sqrtf(diffx*diffx + diffy*diffy + diffz*diffz);
//     fscalor = -s * param[cap];
//   }
//   else if (dr2 <= param[rc] * param[rc]){
//     ScalorType sri = param[bv] * ri;
//     ScalorType sri2 = sri * sri;
//     ScalorType sri6 = sri2*sri2*sri2;
//     fscalor = - 24.f * param[epsilon] * (2.f * sri6*sri6 - sri6) * ri * ri;
//   }
//   else if (dr2 < param[rcut] * param[rcut]) {
//     ScalorType tmp = M_PIF * param[wci];
//     ScalorType term = 0.5f *  (dr2 * ri - param[rc]) * tmp;
//     fscalor = param[epsilon] * tmp *
//     	cosf(term) * sinf(term) * ri;
//   }

//   *fx = diffx * fscalor;
//   *fy = diffy * fscalor;
//   *fz = diffz * fscalor;
// }


// ScalorType
// forcePoten (const ScalorType * param,
// 			 ScalorType diffx,
// 			 ScalorType diffy,
// 			 ScalorType diffz,
// 			 ScalorType *fx, 
// 			 ScalorType *fy, 
// 			 ScalorType *fz)
// {
//   ScalorType dr2 = diffx*diffx + diffy*diffy + diffz*diffz;
//   ScalorType rvalue = 0.f;
//   ScalorType fscalor = 0.f;
//   ScalorType ri = 1.f / sqrtf (dr2);

// //   printf ("capR pcarR is %f %f\n", param[capR], param[pcapR]);

//   if (dr2 == 0.f){
//     *fy = 0.f;
//     *fz = 0.f;
//     return 0.f;
//   }
//   else if (dr2 < param[capR]*param[capR]){
//     ScalorType s = 1.f / sqrtf(diffx*diffx + diffy*diffy + diffz*diffz);
//     fscalor = -s * param[cap];
//     rvalue = -param[cap] * (sqrtf(dr2) - param[capR]) + param[pcapR] ;
//   }
//   else if (dr2 <= param[rc] * param[rc]){
//     ScalorType sri = param[bv] * ri;
//     ScalorType sri2 = sri * sri;
//     ScalorType sri6 = sri2*sri2*sri2;
//     fscalor = - 24.f * param[epsilon] * (2.f * sri6*sri6 - sri6) * ri * ri;
//     rvalue = 4.f * param[epsilon] * (sri6*sri6 - sri6);
//     if (param[wc] == 0.f) rvalue += param[epsilon];
//     else rvalue += 0.f;
//   }
//   else if (dr2 < param[rcut] * param[rcut]) {
//     ScalorType term = 0.5f * M_PIF * (dr2 * ri - param[rc]) * param[wci];
//     ScalorType cost = cosf(term);
//     fscalor = param[epsilon] * M_PIF * param[wci] *
// 	cost * sinf(term) * ri;
//     rvalue = - param[epsilon] * cost * cost;
//   }

//   *fx = diffx * fscalor;
//   *fy = diffy * fscalor;
//   *fz = diffz * fscalor;

//   return rvalue;
// }

// void printfF (ScalorType * param, const char * name)
// {
//   FILE * fp = fopen (name, "w");
//   ScalorType fx, fy, fz;
//   ScalorType x;
//   ScalorType lower = 0.8, upper = 3;
//   IndexType N = 1000;
//   ScalorType h = (upper - lower) / N;

//   for (unsigned i = 0; i <= N; ++i){
//     x = lower + i * h;
//     ScalorType p ;
//     p = forcePoten (param, x, 0., 0., &fx, &fy, &fz);
//     fprintf (fp, "%f\t%f\t%f\n", x, p, fx);
//   }
//   fclose (fp);
// }	     
