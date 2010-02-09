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
#include "BondList_interface.h"

#include "MDSystem_interface.h"
#include "MDException.h"
#include "Topology.h"

//#include "cutil.h"

#define NThreadsPerBlock	5

int main(int argc, char * argv[])
{

  // Topology::Atom lja (1.0, 0., 0);
  
  // return 0;
  
  
  LennardJones6_12Parameter lp;
  lp.init (1, 1, 1, 2.5);
  LennardJones6_12Parameter lp1;
  lp1.init (1, 1, 1, 2.5);
  CosTailParameter cost;
  cost.init (1, 1, 2.4);
  
  printf ("%d %d %d\n", lp.type(), lp.numParam(), lp == lp1);

  SystemNonBondedInteraction sysnb;
  printf ("is built %d\n", sysnb.isBuilt);
  sysnb.clear();
  printf ("is built %d\n", sysnb.isBuilt);

  sysnb.add (1,1,lp);
  sysnb.add (0,1,lp1);
  sysnb.add (0,0,cost);
  sysnb.add (0,2,cost);
  // sysnb.add (2,2,lp);
  printf ("is built %d\n", sysnb.isBuilt);
  sysnb.build();
  printf ("is built %d\n", sysnb.isBuilt);

  for (unsigned i = 0; i < sysnb.numInteractionItems; ++i){
    printf ("type: %d, ", sysnb.types[i]);
    for (unsigned j = sysnb.positions[i]; j < sysnb.positions[i+1]; ++j){
      printf ("%1.3e ", sysnb.parameters[j]);
    }
    printf("\n");
  }
  for (unsigned i = 0; i < sysnb.numAtomTypes; ++i){
    for (unsigned j = 0; j < sysnb.numAtomTypes; ++j){
      printf ("%d ", *sysnb.interactionTableItem(i,j));
    }
    printf("\n");    
  }

  printf ("##################################3\n");

  printf ("is built %d\n", sysnb.isBuilt);
  sysnb.clear();
  printf ("is built %d\n", sysnb.isBuilt);
  sysnb.add (1, 1, lp);
  sysnb.add (1, 0, lp);
  printf ("is built %d\n", sysnb.isBuilt);
  sysnb.build();
  printf ("is built %d\n", sysnb.isBuilt);
  
  for (unsigned i = 0; i < sysnb.numInteractionItems; ++i){
    printf ("type: %d, ", sysnb.types[i]);
    for (unsigned j = sysnb.positions[i]; j < sysnb.positions[i+1]; ++j){
      printf ("%1.3e ", sysnb.parameters[j]);
    }
    printf("\n");
  }
  for (unsigned i = 0; i < sysnb.numAtomTypes; ++i){
    for (unsigned j = 0; j < sysnb.numAtomTypes; ++j){
      printf ("%d ", *sysnb.interactionTableItem(i,j));
    }
    printf("\n");    
  }
  
  return 0;

  try{
    IndexType nstep = 20;
    char * filename;
  
    if (argc < 3){
      printf ("Usage:\n%s conf.gro nstep\n", argv[0]);
      return 1;
    }
    if (argc != 1){
      nstep = atoi(argv[2]);
      filename = argv[1];
    }

    // CUT_DEVICE_INIT (argc, argv);
    
    MDSystem sys;
    sys.initConfig(filename, "atom.map");
    // ScalorType ljparam[mdForceNParamLennardJones6_12];
    // ScalorType ljcapparam [mdForceNParamLennardJones6_12_cap];
  
    // LennardJones6_12::initParameter (ljparam, 1.f, 1.f, 0.f, 2.5f);
    // sys.addNBForce (0, 0, 
    // 		    mdForceLennardJones6_12, 
    // 		    ljparam);
    // LennardJones6_12_cap::initParameter (ljcapparam, 1.f, 1.f, 0.f, 2.5f, 100.f);
    // sys.addNBForce (1, 1,
    // 		    mdForceLennardJones6_12_cap,
    // 		    ljcapparam);
    // LennardJones6_12::initParameter (ljparam, 1.f, 1.f, 0.f, 2.5f);
    // sys.addNBForce (0, 1, 
    // 		    mdForceLennardJones6_12, 
    // 		    ljparam);

    sys.initBond ();
    ScalorType hsparam[mdForceNParamHarmonicSpring] ;
    HarmonicSpring::initParameter (hsparam, 5.f, 1);
    ScalorType feneparam[mdForceNParamFENE];
    FENE::initParameter (feneparam, 20.f, 2.3f);
  
    // for (unsigned i = 0; i < sys.hdata.numAtom; i+=2){
    //   if (i < 50){
    // 	sys.addBond (i, i+1, mdForceHarmonicSpring, hsparam);
    // 	sys.addBond (i, i+1, mdForceFENE, feneparam);
    //   }
    //   else{
    // 	sys.addBond (i, i+1, mdForceFENE, feneparam);
    // 	sys.addBond (i, i+1, mdForceHarmonicSpring, hsparam);
    //   }
    // }   
    // for (unsigned i = 0; i < sys.hdata.numAtom-1; i++){
    //   if (i < 50){
    // 	sys.addBond (i, i+1, mdForceHarmonicSpring, hsparam);
    // 	sys.addBond (i, i+1, mdForceFENE, feneparam);
    //   }
    //   else{
    // 	sys.addBond (i, i+1, mdForceFENE, feneparam);
    // 	sys.addBond (i, i+1, mdForceHarmonicSpring, hsparam);
    //   }
    // }   
    sys.buildBond();

    ScalorType ahparam [mdForceNParamAngleHarmonic];
    AngleHarmonic::initParameter (ahparam, 5.f, 2./3.*M_PI);
    sys.initAngle();
    // for (unsigned i = 0; i < sys.hdata.numAtom-2; ++i){
    //   sys.addAngle (i, i+1, i+2, mdForceAngleHarmonic, ahparam);
    // }
    // sys.addAngle(0, 1, 2, mdForceAngleHarmonic, ahparam);
    // sys.addAngle(1, 2, 3, mdForceAngleHarmonic, ahparam);
    // sys.addAngle(2, 3, 4, mdForceAngleHarmonic, ahparam);
    sys.buildAngle();

    ScalorType maxrcut = sys.calMaxNBRcut ();
    printf ("# max rcut is %f\n", maxrcut);
    ScalorType nlistExten = 0.3f;
    ScalorType rlist = maxrcut + nlistExten;
    NeighborList nlist(sys, rlist, NThreadsPerBlock, 20);;
    nlist.build (sys);
  
    Reshuffle resh (sys, nlist, NThreadsPerBlock);
    resh.shuffleSystem ( sys, nlist);

    MDStatistic st(sys);

    VelocityVerlet inte (sys, NThreadsPerBlock);;

    VelocityRescale inte_vr (sys, NThreadsPerBlock, 1, 0.1);
// // printf ("%f %f\n", ddata.velox[0], ddata.velox[1]);
    // inte.removeTranslationalFreedom (ddata);
    // // printf ("%f %f\n", ddata.velox[0], ddata.velox[1]);

    TranslationalFreedomRemover tfremover (sys, NThreadsPerBlock);

    InteractionEngine_interface interaction(sys, NThreadsPerBlock);;
  
    ScalorType dt = 0.001;
    for (unsigned i = 0; i < nstep; ++i){ 
      if (i%10 == 0){
	tfremover.remove (sys);
      }
      // if (i%10 == 0){
      //   if (i == 0)
      // 	nlist.build (ddata, box);
      //   else
      // 	nlist.reBuild (ddata, box);
      //   resh.shuffleSystem (sys, nlist);
      // }
      if (nlist.judgeRebuild(sys, 0.5 * nlistExten)){
	printf("# rebuild at step %d\n", i);
	nlist.reBuild(sys);
      }
	  
      if ((i+1) % 1 == 0){
	st.clearDevice();
	inte.step1 (sys, dt);
	interaction.applyInteraction (sys, nlist,  st);
	inte.step2 (sys, dt, st);
	st.updateHost();
	printf ("%07d %.7e %.7e %.7e %.7e\n",
		i+1, 
		st.getStatistic(mdStatisticNonBondedPotential),
		st.getStatistic(mdStatisticBondedPotential),
		st.kineticEnergy(),
		st.getStatistic(mdStatisticNonBondedPotential) +
		st.getStatistic(mdStatisticBondedPotential) +
		st.kineticEnergy());
	fflush(stdout);
        resh.shuffleSystem (sys, nlist);
      }
      else {
	inte.step1 (sys, dt);
	interaction.applyInteraction (sys, nlist);
	inte.step2 (sys, dt);
      }
    }
  
    // resh.recoverMDData (ddata, ddata2);
  }
  // catch (MDExcptCannotOpenFile *e){
  //   fprintf (stderr, "%s\n", e->what());
  //   return 1;
  // }
  catch (MDException &e){
    fprintf (stderr, "%s\n", e.what());
    return 1;
  }
    
  return 0;
}
