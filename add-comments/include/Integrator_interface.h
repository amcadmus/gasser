#ifndef __Integrator_interface_h_wanghan__
#define __Integrator_interface_h_wanghan__

#include "MDSystem_interface.h"
#include "Integrator.h"
#include "MDSystem.h"
#include "Statistic_interface.h"
#include "MDTimer_interface.h"
#include "SumVector.h"
#include "InteractionEngine_interface.h"
#include "Thermostat.h"
#include "Barostat.h"

/// Remove the non-zero moment of the system.

class TranslationalFreedomRemover 
{
  dim3 atomGridDim;
  dim3 myBlockDim;
  ScalorType * sums;
  IndexType sharedBuffSize;
  SumVector<ScalorType > sum_x;
  SumVector<ScalorType > sum_y;
  SumVector<ScalorType > sum_z;
public:
  /** 
   * Constructor.
   * 
   * @param sys Remove the non-zero moment of this system, which
   * should be the same as that used in the TranslationalFreedomRemover::remove.
   * @param NThread Number of threads in a block.
   */
  TranslationalFreedomRemover (const MDSystem & sys, 
			       const IndexType & NThread)
      { init (sys, NThread); }
  ~TranslationalFreedomRemover ();
  /** 
   * Initializer.
   * 
   * @param sys Remove the non-zero moment of this system, which
   * should be the same as that used in the TranslationalFreedomRemover::remove.
   * @param NThread Number of threads in a block.
   */
  void init (const MDSystem & sys, 
	     const IndexType & NThread);
public:
  /** 
   * Remove the non-zero moment of the system.
   * 
   * @param sys Remove the non-zero moment of this system.
   * @param timer Timer measuring the performance.
   */
  void remove (MDSystem & sys, MDTimer * timer = NULL);
};

/// The leap frog integrator.  
/**
 * The formula of the leap frog integrator contain two steps, the
 * velocity step:
 * \f[
 * v_i (t+\frac{\Delta t}2) = v_i (t-\frac{\Delta t}2) + \Delta t\frac1{m_i}F_i(t)
 * \f]
 * and the position step:
 * \f[
 * r_i (t+{\Delta t}) = r_i (t) + \Delta t\, v_i (t+\frac{\Delta t}2)
 * \f]
 */

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
  /** 
   * Constructor.
   * 
   * @param sys The MDSystem to integrate, which should be the same as
   * that used in the following calls.
   * @param NThread Number of threads in a block.
   */
  LeapFrog (const MDSystem & sys, 
	    const IndexType & NThread)
      { reinit (sys, NThread);}
  /** 
   * Reinitializer.
   * 
   * @param sys The MDSystem to integrate, which should be the same as
   * that used in the following calls.
   * @param NThread Number of threads in a block.
   */
  void reinit (const MDSystem & sys, 
	       const IndexType & NThread);
public:
  /** 
   * Integrate one step.
   * 
   * @param sys The MDSystem to integrate.
   * @param dt The time step.
   * @param timer Timer measuring the performance.
   */
  void step (MDSystem & sys,
	     const ScalorType & dt,
	     MDTimer * timer = NULL);
  /** 
   * Integrate one step with statistic.
   * 
   * @param sys The MDSystem to integrate.
   * @param dt The time step.
   * @param st The statistic.
   * @param timer Timer measuring the performance.
   */
  void step (MDSystem & sys,
	     const ScalorType & dt,
	     MDStatistic & st,
	     MDTimer * timer = NULL);
public:
  /** 
   * Integrate the position forward.
   * 
   * @param sys The MDSystem to integrate.
   * @param dt The time step.
   * @param timer Timer measuring the performance.
   */
  void stepX (MDSystem & sys,
	      const ScalorType & dt,
	      MDTimer * timer = NULL);
  /** 
   * Integrate the velocity forward.
   * 
   * @param sys The MDSystem to integrate.
   * @param dt The time step.
   * @param timer Timer measuring the performance.
   */
  void stepV (MDSystem & sys,
	      const ScalorType & dt,
	      MDTimer * timer = NULL);
  /** 
   * Integrate the velocity forward, with statistic.
   * 
   * @param sys The MDSystem to integrate.
   * @param dt The time step.
   * @param st The statistic.
   * @param timer Timer measuring the performance.
   */
  void stepV (MDSystem & sys,
	      const ScalorType & dt,
	      MDStatistic & st,
	      MDTimer * timer = NULL);
  void stepV_VCouple (MDSystem & sys,
		      const ScalorType & dt,
		      const ScalorType * lambda,
		      MDTimer * timer = NULL);
  void stepV_VCouple (MDSystem & sys,
		      const ScalorType & dt,
		      const ScalorType * lambda,
		      MDStatistic & st,
		      MDTimer * timer = NULL);
public:
  // void oneStep (MDSystem & sys, const RectangularBox & box,
  // 		const NeighborList & nlist, const ScalorType & dt);
  // void oneStep (MDSystem & sys, const RectangularBox & box,
  // 		const NeighborList & nlist, const ScalorType & dt);
};

class LeapFrog_TPCouple_Rescale
{
  dim3 atomGridDim;
  dim3 myBlockDim;
// private:
//   IndexType NPcoupleGroup;
//   PcoupleDirection_t PcoupleDirections[3];
//   IndexType NDirInGroup[3];
//   IndexType DirIndex[3][3];
private:
  ScalorType			dt;
  IndexType			nstep;
  MDStatistic			myst;
  LeapFrog			lpfrog;
  InteractionEngine *		ptr_inter;
  CellList *			ptr_clist;
  NeighborList *		ptr_nlist;
  BondedInteractionList *	ptr_bdInterList;
  ScalorType			rebuildThreshold;
  void nullPointers ();
private:
  const Thermostat_VRescale * ptr_thermostat;
  const Barostat_XRescale * ptr_barostat;
private:
  void firstStep (MDSystem & sys, MDTimer * timer);
  void firstStep (MDSystem & sys, MDStatistic &st, MDTimer * timer);
public:
  void init (const MDSystem &sys,
	     const IndexType & NThread,
	     const ScalorType & dt,
	     InteractionEngine &inter,
	     NeighborList & nlist,
	     const ScalorType & rebuildThreshold,
	     BondedInteractionList * ptr_bdInterList = NULL) ;
  LeapFrog_TPCouple_Rescale ();
  LeapFrog_TPCouple_Rescale (const MDSystem &sys,
		     const IndexType & NThread,
		     const ScalorType & dt,
		     InteractionEngine &inter,
		     NeighborList & nlist,
		     const ScalorType & rebuildThreshold,
		     BondedInteractionList * ptr_bdInterList = NULL)
      : myst(sys), lpfrog(sys, NThread)
      { init (sys, NThread, dt, inter, nlist, rebuildThreshold, ptr_bdInterList); }
  ~LeapFrog_TPCouple_Rescale () {};
public:
  void addThermostat (const Thermostat_VRescale & thermostat);
  void disableThermostat ();
  void addBarostat (const Barostat_XRescale & barostat);
  void disableBarostat ();
  void oneStep (MDSystem & sys, MDTimer * timer=NULL);
  void oneStep (MDSystem & sys, MDStatistic &st, MDTimer * timer=NULL);
};

class LeapFrog_TPCouple_VCouple
{
  dim3 atomGridDim;
  dim3 myBlockDim;
// private:
//   IndexType NPcoupleGroup;
//   PcoupleDirection_t PcoupleDirections[3];
//   IndexType NDirInGroup[3];
//   IndexType DirIndex[3][3];
private:
  bool				isFirstStep;
  LeapFrog			lpfrog;
  void nullPointers ();
private:
  const Thermostat_VCouple *	ptr_thermostat;
  const Barostat_VCouple *	ptr_barostat;
private:
  // void firstStep (MDSystem & sys, MDTimer * timer);
  void firstStep (MDSystem & sys,
		  const ScalorType & dt,
		  MDStatistic &st,
		  MDTimer * timer);
public:
  LeapFrog_TPCouple_VCouple ();
  LeapFrog_TPCouple_VCouple (const MDSystem &sys,
			     const IndexType & NThread);
  ~LeapFrog_TPCouple_VCouple () {};
  void reinit (const MDSystem &sys,
	       const IndexType & NThread) ;
public:
  void addThermostat (const Thermostat_VCouple & thermostat);
  void disableThermostat ();
  void addBarostat (const Barostat_VCouple & barostat);
  void disableBarostat ();
public:
  void reset ();
  // void oneStep (MDSystem & sys,
  // 		MDTimer * timer=NULL);
  void oneStep (MDSystem & sys,
		const ScalorType & dt,
		MDStatistic &lastStepSt,
		MDStatistic &st,
		MDTimer * timer=NULL);
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
  InteractionEngine * ptr_inter;
  NeighborList * ptr_nlist;
  BondedInteractionList * ptr_bdInterList;
  ScalorType rebuildThreshold;
private:
  void firstStep (MDSystem & sys, MDTimer * timer);
  void firstStep (MDSystem & sys, MDStatistic &st, MDTimer * timer);
public:
  void init (const MDSystem &sys,
	     const IndexType & NThread,
	     const ScalorType & dt,
	     InteractionEngine &inter,
	     NeighborList & nlist,
	     const ScalorType & rebuildThreshold,
	     BondedInteractionList * ptr_bdInterList = NULL) ;
  BerendsenLeapFrog ();
  BerendsenLeapFrog (const MDSystem &sys,
		     const IndexType & NThread,
		     const ScalorType & dt,
		     InteractionEngine &inter,
		     NeighborList & nlist,
		     const ScalorType & rebuildThreshold,
		     BondedInteractionList * ptr_bdInterList = NULL)
      : myst(sys), lpfrog(sys, NThread)
      { init (sys, NThread, dt, inter, nlist, rebuildThreshold, ptr_bdInterList); }
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
