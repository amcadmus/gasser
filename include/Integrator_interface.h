#ifndef __Integrator_interface_h_wanghan__
#define __Integrator_interface_h_wanghan__

#include "Integrator.h"
#include "MDSystem.h"
#include "Statistic_interface.h"
#include "MDTimer_interface.h"
#include "SumVector.h"
#include "InteractionEngine_interface.h"


class TranslationalFreedomRemover 
{
  dim3 atomGridDim;
  dim3 myBlockDim;
  ScalorType * sums;
  ScalorType * buffx, * buffy, * buffz;
public:
  TranslationalFreedomRemover (const MDSystem & sys, 
			       const IndexType & NThread)
      { init (sys, NThread); }
  ~TranslationalFreedomRemover ();
  void init (const MDSystem & sys, 
	     const IndexType & NThread);
public:
  void remove (MDSystem & sys, MDTimer * timer = NULL);
};

  

class LeapFrog 
{
  dim3 atomGridDim;
  dim3 myBlockDim;
  // ScalorType * statistic_buffxx;
  // ScalorType * statistic_buffyy;
  // ScalorType * statistic_buffzz;
  SumVector<ScalorType> sum_kxx;
  SumVector<ScalorType> sum_kyy;
  SumVector<ScalorType> sum_kzz;
  IndexType sharedBuffSize;
public:
  LeapFrog () {}
  ~LeapFrog ();
  LeapFrog (const MDSystem & sys, 
	    const IndexType & NThread)
      { init (sys, NThread);}
  void init (const MDSystem & sys, 
	     const IndexType & NThread);
public:
  void step (MDSystem & sys, const ScalorType & dt, MDTimer * timer = NULL);
  void step (MDSystem & sys, const ScalorType & dt,
	     MDStatistic & st, MDTimer * timer = NULL);
public:
  void stepX (MDSystem & sys, const ScalorType & dt, MDTimer * timer = NULL);
  void stepV (MDSystem & sys, const ScalorType & dt, MDTimer * timer = NULL);
  void stepV (MDSystem & sys, const ScalorType & dt,
	      MDStatistic & st, MDTimer * timer = NULL);
public:
  // void oneStep (MDSystem & sys, const RectangularBox & box,
  // 		const NeighborList & nlist, const ScalorType & dt);
  // void oneStep (MDSystem & sys, const RectangularBox & box,
  // 		const NeighborList & nlist, const ScalorType & dt);
};


class VelocityVerlet 
{
  dim3 atomGridDim;
  dim3 myBlockDim;
  // ScalorType * statistic_buffxx;
  // ScalorType * statistic_buffyy;
  // ScalorType * statistic_buffzz;
  SumVector<ScalorType> sum_kxx;
  SumVector<ScalorType> sum_kyy;
  SumVector<ScalorType> sum_kzz;
  IndexType sharedBuffSize;
public:
  ~VelocityVerlet ();
  VelocityVerlet (const MDSystem & sys, 
		  const IndexType & NThread)
      { init (sys, NThread);}
  void init (const MDSystem & sys, 
	     const IndexType & NThread);
public:
  void step1 (MDSystem & sys, const ScalorType & dt, MDTimer * timer = NULL);
  void step2 (MDSystem & sys, const ScalorType & dt, MDTimer * timer = NULL);
  void step2 (MDSystem & sys, const ScalorType & dt,
	      MDStatistic & st, MDTimer * timer = NULL);
public:
  // void oneStep (MDSystem & sys, const RectangularBox & box,
  // 		const NeighborList & nlist, const ScalorType & dt);
  // void oneStep (MDSystem & sys, const RectangularBox & box,
  // 		const NeighborList & nlist, const ScalorType & dt);
};



class VelocityRescale
{
  dim3 atomGridDim;
  dim3 myBlockDim;
  ScalorType * kineticE;
  ScalorType hkineticE;
  ScalorType * buff;
  ScalorType Nf;
  ScalorType tau;
  ScalorType scalor1;
  ScalorType scalor2;
  ScalorType tmp2;
  ScalorType refK;
  ScalorType newK;
  ScalorType alpha;
  SumVector<ScalorType> sum_kxx;
  SumVector<ScalorType> sum_kyy;
  SumVector<ScalorType> sum_kzz;
  SumVector<ScalorType> sum_k;
  IndexType sharedBuffSize;
public:
  ~VelocityRescale();
  VelocityRescale (const MDSystem & sys, 
		   const IndexType & NThread,
		   const ScalorType & refT,
		   const ScalorType & tau)
      { init (sys, NThread, refT, tau);}
  void init (const MDSystem & sys, 
	     const IndexType & NThread,
	     const ScalorType & refT,
	     const ScalorType & tau);
public:
  void step1 (MDSystem & sys, const ScalorType & dt, MDTimer * timer = NULL);
  void step2 (MDSystem & sys, const ScalorType & dt, MDTimer * timer = NULL);
  void step2 (MDSystem & sys, const ScalorType & dt,
	      MDStatistic & st , MDTimer * timer = NULL);
};


enum PcoupleDirection {
  PCoupleX		= 1,
  PCoupleY		= 2,
  PCoupleZ		= 4,
  PCoupleXY		= 3,
  PCoupleXYZ		= 7
};
typedef int PCoupleDirection_t;

class BerendsenLeapFrog 
{
  dim3 atomGridDim;
  dim3 myBlockDim;
  bool TCoupleOn;
  bool PCoupleOn;
  ScalorType refT;
  ScalorType tauT;
  ScalorType refP[3];
  ScalorType tauP[3];
  ScalorType betaP[3];
  IndexType NPCoupleGroup;
  PCoupleDirection_t PCoupleDirections[3];
private:
  ScalorType dt;
  IndexType nstep;
  MDStatistic myst;
  ScalorType lastKineticE;
  ScalorType lastPressure[3];
  LeapFrog lpfrog;
  InteractionEngine_interface * ptr_inter;
  NeighborList * ptr_nlist;
  ScalorType rebuildThreshold;
private:
  void firstStep (MDSystem & sys, MDTimer * timer);
  void firstStep (MDSystem & sys, MDStatistic &st, MDTimer * timer);
public:
  void init (const MDSystem &sys,
	     const IndexType & NThread,
	     const ScalorType & dt,
	     InteractionEngine_interface &inter,
	     NeighborList & nlist,
	     const ScalorType & rebuildThreshold) ;
  BerendsenLeapFrog ();
  BerendsenLeapFrog (const MDSystem &sys,
		     const IndexType & NThread,
		     const ScalorType & dt,
		     InteractionEngine_interface &inter,
		     NeighborList & nlist,
		     const ScalorType & rebuildThreshold)
      : myst(sys), lpfrog(sys, NThread)
      { init (sys, NThread, dt, inter, nlist, rebuildThreshold); }
  ~BerendsenLeapFrog () {};
public:
  void TCouple (const ScalorType & refT,
		const ScalorType & tauT);
  void addPcoupleGroup (const PCoupleDirection_t & direction,
			const ScalorType & refP,
			const ScalorType & tauP,
			const ScalorType & betaP);
  void oneStep (MDSystem & sys, MDTimer * timer=NULL);
  void oneStep (MDSystem & sys, MDStatistic &st, MDTimer * timer=NULL);
};




  


#endif
