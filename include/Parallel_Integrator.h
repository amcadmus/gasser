#ifndef __Parallel_Integrator_h_wanghan__
#define __Parallel_Integrator_h_wanghan__

#define DEVICE_CODE
#include "SumVector.h"
#include "Parallel_CellList.h"
#include "Parallel_Statistic.h"

namespace Parallel {
  class TranslationalFreedomRemover 
  {
    dim3 gridDim;
    IndexType numThreadsInCell;
    ScalorType * sums;
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
      VelocityVerlet (const DeviceCellListedMDData & sys) { init (sys);}
      void init (const DeviceCellListedMDData & sys);
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
