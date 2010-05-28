#ifndef __Parallel_Integrator_h_wanghan__
#define __Parallel_Integrator_h_wanghan__

#include "common.h"

#ifdef DEVICE_CODE
#include "Parallel_CellList.h"
#include "Parallel_Statistic.h"
#include "SumVector.h"
#endif

namespace Parallel {
  class HostSystemMomentum
  {
    double * p;
    double * sump;
    ScalorType * p_f;
public:
    HostSystemMomentum ();
    ~HostSystemMomentum ();
public:
    void sumAll ();
    void setMomentunX (const ScalorType & x) {p_f[0] = x;}
    void setMomentunY (const ScalorType & y) {p_f[1] = y;}
    void setMomentunZ (const ScalorType & z) {p_f[2] = z;}
    const ScalorType & getMomentumX () const {return p_f[0];}
    const ScalorType & getMomentumY () const {return p_f[1];}
    const ScalorType & getMomentumZ () const {return p_f[2];}
  };
  class HostSystemMass
  {
    double mass;
    double sumMass;
    ScalorType p_mass;
public:
    void sumAll();
    void setMass (const ScalorType & x) {p_mass = x;}
    const ScalorType & getMass () const {return p_mass;}
  };
}

#ifdef DEVICE_CODE

namespace Parallel {
  class TranslationalFreedomRemover 
  {
    dim3 gridDim;
    IndexType numThreadsInCell;
    ScalorType * sums;
    ScalorType * hsums;
    HostSystemMomentum hmomentum;
    ScalorType * sumM;
    ScalorType totalMassi;
    bool malloced;
    IndexType sharedBuffSize;
    SumVector<ScalorType > sum_x;
    SumVector<ScalorType > sum_y;
    SumVector<ScalorType > sum_z;
    void clear ();
public:
    TranslationalFreedomRemover ()
	: malloced (false) {}
    TranslationalFreedomRemover (const DeviceCellListedMDData & data)
	: malloced (false) { reinit (data);}
    ~TranslationalFreedomRemover ();
    void reinit (const DeviceCellListedMDData & data);
public:
    void remove (DeviceCellListedMDData & data);
  };
    
  namespace Integrator {
    class VelocityVerlet 
    {
      dim3 gridDim;
      SumVector<ScalorType> sum_kxx;
      SumVector<ScalorType> sum_kyy;
      SumVector<ScalorType> sum_kzz;
      IndexType sharedBuffSize;
  public:
      VelocityVerlet () {};
      VelocityVerlet (const DeviceCellListedMDData & sys) {reinit (sys);}
      void reinit (const DeviceCellListedMDData & sys);
  public:
      void step1 (DeviceCellListedMDData & sys,
		  const ScalorType & dt);
      void step2 (DeviceCellListedMDData & sys,
		  const ScalorType & dt);
      void step2 (DeviceCellListedMDData & sys,
		  const ScalorType & dt,
		  DeviceStatistic & st);
    };
  }


  class BerendsenLeapFrog 
  {
    dim3 gridDim;
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
    BondedInteractionList * ptr_bdInterList;
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
	       const ScalorType & rebuildThreshold,
	       BondedInteractionList * ptr_bdInterList = NULL) ;
    BerendsenLeapFrog ();
    BerendsenLeapFrog (const MDSystem &sys,
		       const IndexType & NThread,
		       const ScalorType & dt,
		       InteractionEngine_interface &inter,
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


  namespace CudaGlobal{
    __global__ void
    prepareCalTotalMass (const IndexType * numAtomInCell,
			 const ScalorType * mass,
			 ScalorType * buff);
    __global__ void
    prepareRemoveTranslationalFreedom (const IndexType * numAtomInCell,
				       const ScalorType * mass,
				       const ScalorType * velox,
				       const ScalorType * veloy,
				       const ScalorType * veloz,
				       ScalorType * st_buff_x,
				       ScalorType * st_buff_y,
				       ScalorType * st_buff_z);
    __global__ void
    removeTranslationalFreedom (const IndexType * numAtomInCell,
				const ScalorType totalMassi,
				const ScalorType * sums,
				ScalorType * velox,
				ScalorType * veloy,
				ScalorType * veloz);
    __global__ void
    velocityVerlet_step1 (const IndexType * numAtomInCell,
			  const ScalorType * forcx,
			  const ScalorType * forcy,
			  const ScalorType * forcz,
			  const ScalorType * mass,
			  const ScalorType dt,
			  CoordType * coord,
			  ScalorType * velox,
			  ScalorType * veloy,
			  ScalorType * veloz);
    __global__ void
    velocityVerlet_step2 (const IndexType * numAtomInCell,
			  const ScalorType * forcx,
			  const ScalorType * forcy,
			  const ScalorType * forcz,
			  const ScalorType * mass,
			  const ScalorType   dt,
			  ScalorType * velox,
			  ScalorType * veloy,
			  ScalorType * veloz);
    __global__ void 
    velocityVerlet_step2 (const IndexType * numAtomInCell,
			  const ScalorType * forcx,
			  const ScalorType * forcy,
			  const ScalorType * forcz,
			  const ScalorType * mass,
			  const ScalorType   dt,
			  ScalorType * velox,
			  ScalorType * veloy,
			  ScalorType * veloz,
			  ScalorType * statistic_buffxx,
			  ScalorType * statistic_buffyy,
			  ScalorType * statistic_buffzz);
  }
}

#endif

#endif
