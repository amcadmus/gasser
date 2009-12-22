/**
 * @file   lipid.headDiff.cu
 * @author Han Wang <han_wang@math.pku.edu.cn>
 * @date   Sun Nov 22 12:53:45 2009
 * 
 * @brief  For lipid with different head radius.
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


#define NThreadsPerBlockCell	256
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
  sys.initConfig(filename, "lipid.headDiff.map");
  sys.initNBForce(8);

  ScalorType headDiff = 0.2f;
  ScalorType head1d = 0.95f - 0.5 * headDiff;
  ScalorType head2d = 0.95f + 0.5 * headDiff;
  ScalorType taild = 1.0f;
  ScalorType cosparam[mdForceNParamCosTail];  
  ScalorType strongw = 1.6f;
  ScalorType weekw = 1.5f;
  
  CosTail::initParameter (cosparam, 1.f, head1d, 0.f);
  sys.addNBForce (0, 0, mdForceCosTail, cosparam);
  CosTail::initParameter (cosparam, 1.f, head2d, 0.f);
  sys.addNBForce (3, 3, mdForceCosTail, cosparam);
  CosTail::initParameter (cosparam, 1.f, 0.5*(head1d+head2d), 0.f);
  sys.addNBForce (0, 3, mdForceCosTail, cosparam);

  CosTail::initParameter (cosparam, 1.f, 0.5 *(head1d+taild), 0.f);
  sys.addNBForce (0, 1, mdForceCosTail, cosparam);
  sys.addNBForce (0, 2, mdForceCosTail, cosparam);
  sys.addNBForce (0, 4, mdForceCosTail, cosparam);
  sys.addNBForce (0, 5, mdForceCosTail, cosparam);
  sys.addNBForce (0, 6, mdForceCosTail, cosparam);
  sys.addNBForce (0, 7, mdForceCosTail, cosparam);
  CosTail::initParameter (cosparam, 1.f, 0.5 *(head2d+taild), 0.f);
  sys.addNBForce (3, 1, mdForceCosTail, cosparam);
  sys.addNBForce (3, 2, mdForceCosTail, cosparam);
  sys.addNBForce (3, 4, mdForceCosTail, cosparam);
  sys.addNBForce (3, 5, mdForceCosTail, cosparam);
  sys.addNBForce (3, 6, mdForceCosTail, cosparam);
  sys.addNBForce (3, 7, mdForceCosTail, cosparam);
  
  CosTail::initParameter (cosparam, 1.f, taild, strongw);
  sys.addNBForce (1, 1, mdForceCosTail, cosparam);
  sys.addNBForce (2, 2, mdForceCosTail, cosparam);
  sys.addNBForce (4, 4, mdForceCosTail, cosparam);
  sys.addNBForce (5, 5, mdForceCosTail, cosparam);
  sys.addNBForce (6, 6, mdForceCosTail, cosparam);
  sys.addNBForce (7, 7, mdForceCosTail, cosparam);

  sys.addNBForce (1, 2, mdForceCosTail, cosparam);
  sys.addNBForce (2, 4, mdForceCosTail, cosparam);
  sys.addNBForce (5, 6, mdForceCosTail, cosparam);
  sys.addNBForce (6, 7, mdForceCosTail, cosparam);

  CosTail::initParameter (cosparam, 1.f, taild, weekw);
  sys.addNBForce (1, 5, mdForceCosTail, cosparam);
  sys.addNBForce (2, 6, mdForceCosTail, cosparam);
  sys.addNBForce (2, 5, mdForceCosTail, cosparam);
  sys.addNBForce (1, 6, mdForceCosTail, cosparam);

  sys.addNBForce (4, 7, mdForceCosTail, cosparam);
  sys.addNBForce (4, 6, mdForceCosTail, cosparam);
  sys.addNBForce (2, 7, mdForceCosTail, cosparam);

  CosTail::initParameter (cosparam, 1.f, taild, 0.f);
  sys.addNBForce (1, 4, mdForceCosTail, cosparam);
  sys.addNBForce (1, 7, mdForceCosTail, cosparam);
  sys.addNBForce (5, 4, mdForceCosTail, cosparam);
  sys.addNBForce (5, 7, mdForceCosTail, cosparam);  

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
  ScalorType nlistExten = 0.6f;
  ScalorType rlist = maxrcut + nlistExten;
  // NeighborList nlist(sys, rlist, NThreadsPerBlockCell, 20,
  // 		     RectangularBoxGeometry::mdRectBoxDirectionX |
  // 		     RectangularBoxGeometry::mdRectBoxDirectionY |
  // 		     RectangularBoxGeometry::mdRectBoxDirectionZ);;
  NeighborList nlist(sys, rlist, NThreadsPerBlockCell, 100,
  		     RectangularBoxGeometry::mdRectBoxDirectionX |
  		     RectangularBoxGeometry::mdRectBoxDirectionY);
  
  nlist.build (sys);
  
  Reshuffle resh (sys, nlist, NThreadsPerBlockCell);
  resh.shuffleSystem ( sys, nlist);

  MDStatistic st(sys);

  InteractionEngine_interface interaction(sys, NThreadsPerBlockAtom);;

  VelocityVerlet inte (sys, NThreadsPerBlockAtom);;

  ScalorType dt = 0.005;
  ScalorType rebuildThreshold = 0.5 * nlistExten;
  ScalorType refT = 1.0808863f;
  BerendsenLeapFrog blpf (sys, NThreadsPerBlockAtom, dt,
			  interaction,
			  nlist, rebuildThreshold);
  blpf.TCouple (refT, 0.1);
  blpf.addPcoupleGroup (PCoupleX,
  			0., 20, 23);

  TranslationalFreedomRemover tfremover (sys, NThreadsPerBlockAtom);

  MDTimer timer;
  unsigned i;

  printf ("# prepare ok, start to run\n");
  printf ("#*     1     2           3        4         5       6           7           8           9        10    11    12    13\n");
  printf ("#* nstep  time  nonBondedE  bondedE  kineticE  totalE  pressurexx  pressureyy  pressurezz  pressure  boxx  boxy  boxz\n");
  
  try{
    timer.tic(mdTimeTotal);
    sys.initWriteXtc ("traj.xtc");
    sys.writeHostDataXtc (0, 0*dt, &timer);
    for (i = 0; i < nstep; ++i){ 
      if (i%10 == 0){
	tfremover.remove (sys, &timer);
      }
	  
      if ((i+1) % 200 == 0){
	st.clearDevice();
	blpf.oneStep (sys, st, &timer);
	st.updateHost();
	printf ("%09d %07e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e\n",
		(i+1),  
		(i+1) * dt, 
		st.getStatistic(mdStatisticNonBondedPotential),
		st.getStatistic(mdStatisticBondedPotential),
		st.kineticEnergy(),
		st.getStatistic(mdStatisticNonBondedPotential) +
		st.getStatistic(mdStatisticBondedPotential) +
		st.kineticEnergy(),
		st.pressureXX(),
		st.pressureYY(),
		st.pressureZZ(),
		st.pressure(),
		sys.box.size.x,
		sys.box.size.y,
		sys.box.size.z);
	fflush(stdout);
      }
      else {
	blpf.oneStep (sys, &timer);
      }
      if ((i+1) % 20000 == 0){
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
































// #include <stdio.h>
// #include "common.h"
// #include "BoxGeometry.h"
// #include "MDSystem.h"
// #include "RandomGenerator.h"
// #include "Auxiliary.h"
// #include "NeighborList_interface.h"
// #include"Statistic.h"
// #include "Integrator_interface.h"
// #include "InteractionEngine_interface.h"
// #include "tmp.h"
// #include "Reshuffle_interface.h"
// #include "BondList_interface.h"

// #include "MDSystem_interface.h"


// #define NThreadsPerBlock	32

// namespace tmpcos{
//     typedef enum paramIndex {
//       epsilon		= 0,
//       bv		= 1,
//       rc		= 2,
//       wc		= 3,
//       wci		= 4,
//       rcut		= 5
//     } paramIndex_t;
    
//     ScalorType tmpfp (const ScalorType * param,
// 		      ScalorType diffx,
// 		      ScalorType diffy,
// 		      ScalorType diffz,
// 		      ScalorType *fx, 
// 		      ScalorType *fy,
// 		      ScalorType *fz);
// };

// ScalorType tmpcos::tmpfp (const ScalorType * param,
// 			  ScalorType diffx,
// 			  ScalorType diffy,
// 			  ScalorType diffz,
// 			  ScalorType *fx, 
// 			  ScalorType *fy,
// 			  ScalorType *fz)
// {
//   ScalorType dr2 = diffx*diffx + diffy*diffy + diffz*diffz;
//   if (dr2 <= param[rc] * param[rc]){
//     ScalorType ri2 =1.f/(dr2);
//     ScalorType sri2 = param[bv] * param[bv] * ri2;
//     ScalorType sri6 = sri2*sri2*sri2;
//     ScalorType scalor = 4.f * param[epsilon] * (12.f * sri6*sri6 - 6.f * sri6) * ri2;
//     *fx = diffx * scalor;
//     *fy = diffy * scalor;
//     *fz = diffz * scalor;
//     ScalorType tmp = 4.f * param[epsilon] * (sri6*sri6 - sri6);
//     if (param[wc] == 0.f){
//       return tmp + param[epsilon];
//     }
//     else {
//       return tmp;
//     }
//   }
//   else if (dr2 < param[rcut] * param[rcut]) {
//     ScalorType r = sqrtf(dr2);
//     ScalorType ri = 1.f/(r);
//     ScalorType term = 0.5f * M_PI * (r - param[rc]) * param[wci];
//     ScalorType cost = cosf(term);
//     ScalorType sint = sqrtf (1.f - cost*cost);
//     ScalorType scalor = - param[epsilon] * M_PI * param[wci] * cost * sint * ri;
//     *fx = diffx * scalor;
//     *fy = diffy * scalor;
//     *fz = diffz * scalor;
//     return - param[epsilon] * cost * cost;
//   }
//   else {
//     *fx = *fy = *fz = 0.f;
//     return 0.f;
//   }
// }

// void printpoten (ScalorType * param)
// {
//   ScalorType low = 0.8;
//   ScalorType high = 3.5;
//   ScalorType h = 0.002;
//   IndexType N = (high - low) / h;
//   for (IndexType i = 0; i <= N; ++i){
//     ScalorType x, y, z;
//     ScalorType fx, fy, fz;
//     x = y = z = 0.f;
//     x = low + i * h;
//     ScalorType p = tmpcos::tmpfp (param, x, y, z, &fx, &fy, &fz);
//     printf ("%f  %f %f\n", x, p, fx);
//   }
//   exit(0);
// }    

// int main(int argc, char * argv[])
// {
//   IndexType nstep = 20;
//   char * filename;
  
//   if (argc != 3){
//     printf ("Usage:\n%s conf.gro nstep\n", argv[0]);
//     return 1;
//   }
//   if (argc != 1){
//     nstep = atoi(argv[2]);
//     filename = argv[1];
//   }
  
//   MDSystem sys;
//   sys.initConfig(filename, "lipid.headDiff.map");
//   sys.initNBForce(8);

//   ScalorType headDiff = 0.2f;
//   ScalorType head1d = 0.95f - 0.5 * headDiff;
//   ScalorType head2d = 0.95f + 0.5 * headDiff;
//   ScalorType taild = 1.0f;
//   ScalorType cosparam[mdForceNParamCosTail];  
//   ScalorType strongw = 1.6f;
//   ScalorType weekw = 1.5f;
  
//   CosTail::initParameter (cosparam, 1.f, head1d, 0.f);
//   sys.addNBForce (0, 0, mdForceCosTail, cosparam);
//   CosTail::initParameter (cosparam, 1.f, head2d, 0.f);
//   sys.addNBForce (3, 3, mdForceCosTail, cosparam);
//   CosTail::initParameter (cosparam, 1.f, 0.5*(head1d+head2d), 0.f);
//   sys.addNBForce (0, 3, mdForceCosTail, cosparam);

//   CosTail::initParameter (cosparam, 1.f, 0.5 *(head1d+taild), 0.f);
//   sys.addNBForce (0, 1, mdForceCosTail, cosparam);
//   sys.addNBForce (0, 2, mdForceCosTail, cosparam);
//   sys.addNBForce (0, 4, mdForceCosTail, cosparam);
//   sys.addNBForce (0, 5, mdForceCosTail, cosparam);
//   sys.addNBForce (0, 6, mdForceCosTail, cosparam);
//   sys.addNBForce (0, 7, mdForceCosTail, cosparam);
//   CosTail::initParameter (cosparam, 1.f, 0.5 *(head2d+taild), 0.f);
//   sys.addNBForce (3, 1, mdForceCosTail, cosparam);
//   sys.addNBForce (3, 2, mdForceCosTail, cosparam);
//   sys.addNBForce (3, 4, mdForceCosTail, cosparam);
//   sys.addNBForce (3, 5, mdForceCosTail, cosparam);
//   sys.addNBForce (3, 6, mdForceCosTail, cosparam);
//   sys.addNBForce (3, 7, mdForceCosTail, cosparam);
  
//   CosTail::initParameter (cosparam, 1.f, taild, strongw);
//   sys.addNBForce (1, 1, mdForceCosTail, cosparam);
//   sys.addNBForce (2, 2, mdForceCosTail, cosparam);
//   sys.addNBForce (4, 4, mdForceCosTail, cosparam);
//   sys.addNBForce (5, 5, mdForceCosTail, cosparam);
//   sys.addNBForce (6, 6, mdForceCosTail, cosparam);
//   sys.addNBForce (7, 7, mdForceCosTail, cosparam);

//   sys.addNBForce (1, 2, mdForceCosTail, cosparam);
//   sys.addNBForce (2, 4, mdForceCosTail, cosparam);
//   sys.addNBForce (5, 6, mdForceCosTail, cosparam);
//   sys.addNBForce (6, 7, mdForceCosTail, cosparam);

//   CosTail::initParameter (cosparam, 1.f, taild, weekw);
//   sys.addNBForce (1, 5, mdForceCosTail, cosparam);
//   sys.addNBForce (2, 6, mdForceCosTail, cosparam);
//   sys.addNBForce (2, 5, mdForceCosTail, cosparam);
//   sys.addNBForce (1, 6, mdForceCosTail, cosparam);

//   sys.addNBForce (4, 7, mdForceCosTail, cosparam);
//   sys.addNBForce (4, 6, mdForceCosTail, cosparam);
//   sys.addNBForce (2, 7, mdForceCosTail, cosparam);

//   CosTail::initParameter (cosparam, 1.f, taild, 0.f);
//   sys.addNBForce (1, 4, mdForceCosTail, cosparam);
// //  printpoten (cosparam);
//   sys.addNBForce (1, 7, mdForceCosTail, cosparam);
//   sys.addNBForce (5, 4, mdForceCosTail, cosparam);
//   sys.addNBForce (5, 7, mdForceCosTail, cosparam);  

//   sys.initBond (2);
//   ScalorType hsparam[mdForceNParamHarmonicSpring] ;
//   HarmonicSpring::initParameter (hsparam, 10.f, 4.f);
//   ScalorType feneparam[mdForceNParamFENE];
//   FENE::initParameter (feneparam, 30.f, 1.5f);
  
//   for (unsigned i = 0; i < sys.hdata.numAtom; i+=3){
//     sys.addBond (i, i+1, mdForceFENE, feneparam);
//     sys.addBond (i+2, i+1, mdForceFENE, feneparam);
//     sys.addBond (i, i+2, mdForceHarmonicSpring, hsparam);
//   }   
//   sys.buildBond();

//   ScalorType tmpsum = 0.;
//   for (IndexType i = 0; i < sys.hdata.numAtom; ++i){
//     tmpsum += 0.5 * (sys.hdata.velox[i] * sys.hdata.velox[i] +
// 		     sys.hdata.veloy[i] * sys.hdata.veloy[i] +
// 		     sys.hdata.veloz[i] * sys.hdata.veloz[i] );
//   }
//   printf ("# tmpsum is %f\n", tmpsum);

//   ScalorType maxrcut = sys.calMaxNBRcut ();
//   printf ("# max rcut is %f\n", maxrcut);
//   ScalorType nlistExten = 0.5f;
//   ScalorType rlist = maxrcut + nlistExten;
//   NeighborList nlist(sys, rlist, NThreadsPerBlock, 100);;
//   nlist.build (sys);
  
//   // Reshuffle resh (sys, nlist, NThreadsPerBlock);
//   // resh.shuffleSystem ( sys, nlist);

//   MDStatistic st(sys);

//   VelocityVerlet inte (sys, NThreadsPerBlock);;

//   VelocityRescale inte_vr (sys, NThreadsPerBlock, 1, 0.1);
// // // printf ("%f %f\n", ddata.velox[0], ddata.velox[1]);
//   // inte.removeTranslationalFreedom (ddata);
//   // // printf ("%f %f\n", ddata.velox[0], ddata.velox[1]);

//   TranslationalFreedomRemover tfremover (sys, NThreadsPerBlock);

//   InteractionEngine_interface interaction(sys, NThreadsPerBlock);;

//   MDTimer timer;
//   unsigned i;
//   ScalorType dt = 0.005;

//   try{
//     timer.tic(mdTimeTotal);
//     sys.initWriteXtc ("traj.xtc");
//     sys.writeHostDataXtc (0, 0*dt, &timer);
//     for (i = 0; i < nstep; ++i){ 
//       if (i%10 == 0){
// 	tfremover.remove (sys, &timer);
//       }
//       // if (i%1 == 0){
//       //   if (i == 0) nlist.build (sys);
//       //   else nlist.reBuild (sys);
//       //   // resh.shuffleSystem (sys, nlist);
//       // }
//       if (nlist.judgeRebuild(sys, 0.5 * nlistExten)){
// 	printf ("# Rebuild at step %09i\n", i);
// 	nlist.reBuild(sys, &timer);
//       }
   
	  
//       if ((i+1) % 1 == 0){
// 	st.clearDevice();
// 	inte.step1 (sys, dt), &timer;
// 	interaction.applyInteraction (sys, nlist,  st, &timer);
// 	inte.step2 (sys, dt, st, &timer);
// 	st.updateHost();
// 	printf ("%09d %07e %.7e %.7e %.7e %.7e %.7e\n",
// 		(i+1),  
// 		(i+1) * dt, 
// 		st.getStatistic(mdStatisticNonBondedPotential),
// 		st.getStatistic(mdStatisticBondedPotential),
// 		st.kineticEnergy() +
// 		st.getStatistic(mdStatisticNonBondedPotential) +
// 		st.getStatistic(mdStatisticBondedPotential) +
// 		st.kineticEnergy(),
// 		st.getStatistic(mdStatisticVirialXX) + 
// 		st.getStatistic(mdStatisticVirialYY) +
// 		st.getStatistic(mdStatisticVirialZZ));
// 	fflush(stdout);
//       }
//       else {
// 	inte.step1 (sys, dt, &timer);
// 	interaction.applyInteraction (sys, nlist, &timer);
// 	inte.step2 (sys, dt, &timer);
//       }
//       if ((i+1) % 200 == 0){
// 	sys.updateHost(&timer);
// 	sys.writeHostDataXtc (i+1, (i+1)*dt, &timer);
//       }
//       if ((i+1) % 200 == 0){
// 	// resh.shuffleSystem (sys, nlist, &timer);
//       }
//     }
//     sys.endWriteXtc();
//     timer.toc(mdTimeTotal);
//     timer.printRecord (stderr);
//   }
//   catch (MDExcptCuda & e){
//     sys.updateHost(&timer);
//     sys.writeHostDataXtc (i+1, (i+1)*dt, &timer);
//     timer.toc(mdTimeTotal);
//     timer.printRecord (stderr);
//     return 1;
//   }
//   catch (MDException &e){
//     fprintf (stderr, "%s\n", e.what());
//     return 1;
//   }
  
//     // resh.recoverMDData (ddata, ddata2);

//   return 0;
// }
