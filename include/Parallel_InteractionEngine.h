#ifndef __Parallel_InteractionEngine_h__
#define __Parallel_InteractionEngine_h__

#define DEVICE_CODE

#include "Parallel_CellList.h"
#include "SumVector.h"
#include "Parallel_Statistic.h"
#include "Parallel_BondList.h"
#include "SystemNonBondedInteraction.h"
#include "SystemBondedInteraction.h"

using namespace RectangularBoxGeometry;

namespace Parallel{

#ifdef DEVICE_CODE
  class InteractionEngine
  {
    bool inited ;
    dim3 gridDim;
    IndexType totalNumCell;
    IndexType devideLevel;
    IndexType numThreadsInCell;
    bool hasBond;
    bool hasAngle;
    IndexType calBondInteraction_sbuffSize;
    IndexType calAngleInteraction_sbuffSize;
    IndexType applyNonBondedInteraction_CellList_sbuffSize;
    SumVector<ScalorType> sum_nb_p;
    SumVector<ScalorType> sum_nb_vxx;
    SumVector<ScalorType> sum_nb_vyy;
    SumVector<ScalorType> sum_nb_vzz;
    SumVector<ScalorType> sum_b_p;
    SumVector<ScalorType> sum_b_vxx;
    SumVector<ScalorType> sum_b_vyy;
    SumVector<ScalorType> sum_b_vzz;
    SumVector<ScalorType> sum_angle_p;
    MDError err;
    IndexType * ghostCellIndex;
    IndexType   numGhostCell;
    ScalorType rcut;
private:
    void clear ();
    // void initNonBondedInteraction (const MDSystem & sys);
public:
    InteractionEngine ();
    InteractionEngine (const DeviceCellListedMDData & data) ;
    ~InteractionEngine ();
    void reinit (const DeviceCellListedMDData & data);
    void registNonBondedInteraction (const SystemNonBondedInteraction & sysNbInter);
    void registBondedInteraction    (const SystemBondedInteraction    & sysBdInter);
public:
    ScalorType getMaxRcut () const {return rcut;}
    void clearInteraction (DeviceCellListedMDData & data);
    void applyNonBondedInteraction (DeviceCellListedMDData & data,
				    const DeviceCellRelation & relation);
    void applyNonBondedInteraction (DeviceCellListedMDData & data,
				    const DeviceCellRelation & relation,
				    DeviceStatistic & st);
    void applyBondedInteraction (DeviceCellListedMDData & data,
    				 DeviceBondList & bdlist);
    void applyBondedInteraction (DeviceCellListedMDData & data,
    				 DeviceBondList & bdlist,
    				 DeviceStatistic & st);
  };

  namespace CudaGlobal {
    __global__ void 
    calNonBondedInteraction (const CoordType * coord,
			     const TypeType  * type,
			     const HostVectorType boxSize,
			     const HostVectorType boxSizei,
			     const ScalorType  rlist,
			     const IndexType * numAtomInCell,
			     const IndexType * numNeighborCell,
			     const IndexType * neighborCellIndex,
			     const CoordNoiType * neighborShiftNoi,
			     const IndexType   stride,
			     ScalorType * forcx,
			     ScalorType * forcy,
			     ScalorType * forcz,
			     mdError_t * ptr_de);
    __global__ void 
    calNonBondedInteraction (const CoordType * coord,
			     const TypeType  * type,
			     const HostVectorType boxSize,
			     const HostVectorType boxSizei,
			     const ScalorType  rlist,
			     const IndexType * numAtomInCell,
			     const IndexType * numNeighborCell,
			     const IndexType * neighborCellIndex,
			     const CoordNoiType * neighborShiftNoi,
			     const IndexType   stride,
			     ScalorType * forcx,
			     ScalorType * forcy,
			     ScalorType * forcz,
			     ScalorType * statistic_nb_buff0,
			     ScalorType * statistic_nb_buff1,
			     ScalorType * statistic_nb_buff2,
			     ScalorType * statistic_nb_buff3,
			     mdError_t * ptr_de);
    __global__ void 
    clearGhostBond (const IndexType * ghostCellIndex,
		    IndexType * myNumBond);
    __global__ void 
    calBondInteraction (const CoordType * coord,
			const HostVectorType boxSize,
			const HostVectorType boxSizei,
			const IndexType * numAtomInCell,
			const IndexType * numBond,
			const IndexType * bondIndex,
			const IndexType   bondTopStride,
			const IndexType * bondNeighbor_localIndex,
			const IndexType   bondListStride,
			ScalorType * forcx,
			ScalorType * forcy,
			ScalorType * forcz,
			mdError_t * ptr_de,
			ScalorType * errsrc);
    __global__ void 
    calBondInteraction (const CoordType * coord,
			const HostVectorType boxSize,
			const HostVectorType boxSizei,
			const IndexType * numAtomInCell,
			const IndexType * numBond,
			const IndexType * bondIndex,
			const IndexType   bondTopStride,
			const IndexType * bondNeighbor_localIndex,
			const IndexType   bondListStride,
			ScalorType * forcx,
			ScalorType * forcy,
			ScalorType * forcz,
			ScalorType * statistic_b_buff0,
			ScalorType * statistic_b_buff1,
			ScalorType * statistic_b_buff2,
			ScalorType * statistic_b_buff3,
			mdError_t * ptr_de,
			ScalorType * errsrc);
  }
  namespace CudaDevice {
    __device__ IndexType
    calNonBondedForceIndex (const IndexType * table, 
			    const IndexType numAtomType,
			    const TypeType type0, 
			    const TypeType type1);
  }
  

#endif // DEVICE_CODE
}



#ifdef DEVICE_CODE
static __device__ IndexType Parallel::CudaDevice::
calNonBondedForceIndex (const IndexType * table, 
			const IndexType numType,
			const TypeType atom0, 
			const TypeType atom1);
#endif


#endif
