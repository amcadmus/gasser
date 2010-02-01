#include <stdio.h>
#include "common.h"
#include "BoxGeometry.h"
#include "MDSystem.h"
#include "RandomGenerator.h"
#include "Auxiliary.h"
#include "NeighborList_interface.h"
#include "Statistic.h"
#include "Integrator_interface.h"
#include "InteractionEngine_interface.h"
#include "tmp.h"
#include "Reshuffle_interface.h"
#include "MDTimer_interface.h"

#define NThreadsPerBlockCell	15
#define NThreadsPerBlockAtom	11

int main(int argc, char * argv[])
{   

  IndexType nstep = 20;
  char * filename;
  IndexType reshFeq = 100;
  if (argc < 5){
    printf ("Usage:\n%s nstep conf.gro reshFeq deviceNum\n", argv[0]);
    return 1;
  }
  if (argc != 1){
    nstep = atoi(argv[1]);
    filename = argv[2];
    reshFeq = atoi (argv[3]);
  }

  // CUT_DEVICE_INIT (argc, argv);
  printf ("# setting device to %d\n", atoi(argv[4]));
  cudaSetDevice (atoi(argv[4]));
  checkCUDAError ("set device");
  
  MDSystem sys;
  sys.initConfig(filename, "atom.map");  

  LennardJones6_12Parameter ljparam;
  ljparam.init (1.f, 1.f, 0.f, 3.2f);
  sys.addNonBondedInteraction (1, 1, ljparam);
  sys.buildNonBondedInteraction ();  
  
  if (argc != 1){
    nstep = atoi(argv[1]);
  }

  ScalorType rlist =3.3;
  NeighborList nlist(sys, rlist, NThreadsPerBlockCell, 10);;
  nlist.build (sys);
  
  Reshuffle resh (sys, nlist, NThreadsPerBlockCell);
  resh.shuffleSystem ( sys, nlist);

  MDStatistic st(sys);

  InteractionEngine_interface interaction(sys, NThreadsPerBlockAtom);;

  ScalorType dt = 0.001;
  ScalorType rebuildThreshold = 0.5 *(rlist - 3.2);

  // LeapFrog lpfrog (sys, NThreadsPerBlockAtom);
  // BerendsenLeapFrog blpf (sys, NThreadsPerBlockAtom, dt,
  // 			  interaction,
  // 			  nlist, rebuildThreshold);
  // blpf.TCouple (1, 0.1);
  // blpf.addPcoupleGroup (PCoupleX | PCoupleY ,
  // 			4., 1, 1);
  // blpf.addPcoupleGroup (PCoupleZ,
  // 			4., 10, 1);
  VelocityVerlet inte (sys, NThreadsPerBlockAtom);;
  VelocityRescale inte_vr (sys, NThreadsPerBlockAtom, 1, 0.1);
// // printf ("%f %f\n", ddata.velox[0], ddata.velox[1]);
  // inte.removeTranslationalFreedom (ddata);
  // // printf ("%f %f\n", ddata.velox[0], ddata.velox[1]);

  TranslationalFreedomRemover tfremover (sys, NThreadsPerBlockAtom);


  MDTimer timer;
  

  timer.tic(mdTimeTotal);
  for (unsigned i = 0; i < nstep; ++i){ 
//     if (i%10 == 0){
//       tfremover.remove (ddata, &timer);
//       if (i == 0)
// 	nlist.build (ddata, box, &timer);
//       else 
// 	nlist.reBuild (ddata, box, &timer);
//     }
    if (nlist.judgeRebuild(sys, rebuildThreshold, &timer)){
      printf("# rebuild at step %d\n", i);
      fflush(stdout);
      nlist.reBuild(sys, &timer);
      // nlist.reinit (sys, rlist, NThreadsPerBlockCell, 10);
      // nlist.build (sys);
    }
    if (reshFeq != 0 && i%reshFeq == 0){
      tfremover.remove (sys, &timer);
      resh.shuffleSystem (sys, nlist, &timer);
    }
    if ((i+1)%10 == 0){
      st.clearDevice();
      // blpf.oneStep (sys, st, &timer);
      // interaction.applyInteraction (sys, nlist, st, &timer);
      // lpfrog.step (sys, dt, st, &timer);
      inte.step1 (sys, dt, &timer);
      interaction.applyInteraction (sys, nlist, st, &timer);
      inte.step2 (sys, dt, st, &timer);
      st.updateHost();
      printf ("%07d %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e\n",
	      i+1,
	      st.getStatistic(mdStatisticNonBondedPotential),
	      st.kineticEnergy(),
	      st.getStatistic(mdStatisticNonBondedPotential) +
	      st.kineticEnergy(),
	      st.virial(),
	      st.pressure(),
	      st.pressureXX(),
	      st.pressureYY(),
	      st.pressureZZ(),
	      sys.box.size.x,
	      sys.box.size.y,
	      sys.box.size.z);
      fflush(stdout);
    }
    else {
      // blpf.oneStep (sys, &timer);      
      // interaction.applyInteraction (sys, nlist, &timer);
      // lpfrog.step (sys, dt, &timer);
      inte.step1 (sys, dt, &timer);
      interaction.applyInteraction (sys, nlist, &timer);
      inte.step2 (sys, dt, &timer);
    }    
  }
  timer.toc (mdTimeTotal);
  timer.printRecord (stderr);
  

  // IndexType * hnNei = (IndexType *) malloc (sizeof(IndexType)* nlist.dnlist.stride);
  // cudaMemcpy (hnNei, nlist.dnlist.Nneighbor, sizeof(IndexType) * nlist.dnlist.stride,
  // 	      cudaMemcpyDeviceToHost);
  // FILE * fp = fopen ("nnei.out", "w");
  // for (unsigned i = 0; i < nlist.dnlist.stride; ++i){
  //   fprintf (fp, "%d\n", hnNei[i]);
  // }
  // fclose(fp);
  
  // cpyDeviceMDDataToHost (&ddata, &hdata);
  // writeGroFile ("confout.gro", hdata.numAtom, 
  // 		resdindex, resdname,
  // 		atomname, resdindex,
  // 		hdata.coordx, hdata.coordy, hdata.coordz,
  // 		hdata.velox, hdata.veloy, hdata.veloz,
  // 		box.size.x, box.size.y, box.size.z);
  
  // destroyDeviceMDData(&ddata);
  // destroyHostMDData (&hdata);

  // free (resdindex);
  // free (atomindex);
  // free (resdname);
  // free (atomname);
}


