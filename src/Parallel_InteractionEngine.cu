#define DEVICE_CODE

#include "common.h"
#include "Parallel_InteractionEngine.h"
#include "Parallel_Interface.h"
#include "NonBondedInteraction.h"
#include "compile_error_mixcode.h"

__constant__
InteractionType nonBondedInteractionType [MaxNumberNonBondedInteraction];
__constant__
ScalorType nonBondedInteractionParameter [MaxNumberNonBondedInteractionParameter];
__constant__
IndexType nonBondedInteractionParameterPosition [MaxNumberNonBondedInteraction];
__constant__
IndexType const_nonBondedInteractionTableLength[1];
__constant__
IndexType const_numAtomType[1];
__constant__
IndexType const_nonBondedInteractionTable [MaxLengthNonBondedInteractionTable];


Parallel::InteractionEngine::
InteractionEngine (const DeviceCellListedMDData & ddata)
    : hasBond (false), hasAngle(false)
{
  init (ddata);
}

void Parallel::InteractionEngine::
init (const DeviceCellListedMDData & ddata)
{
  totalNumCell = ddata.getNumCell().x *
      ddata.getNumCell().y * ddata.getNumCell().z;
  gridDim = toGridDim (totalNumCell);
  IndexType numThreadsInCell = Parallel::Interface::numThreadsInCell();
  
  sum_nb_p.init (totalNumCell*numThreadsInCell, NThreadForSum);
  sum_nb_vxx.init (totalNumCell*numThreadsInCell, NThreadForSum);
  sum_nb_vyy.init (totalNumCell*numThreadsInCell, NThreadForSum);
  sum_nb_vzz.init (totalNumCell*numThreadsInCell, NThreadForSum);
  sum_b_p.init (totalNumCell, NThreadForSum);
  sum_b_vxx.init (totalNumCell, NThreadForSum);
  sum_b_vyy.init (totalNumCell, NThreadForSum);
  sum_b_vzz.init (totalNumCell, NThreadForSum);
  sum_angle_p.init (totalNumCell, NThreadForSum);
}


void Parallel::InteractionEngine::
registNonBondedInteraction (const SystemNonBondedInteraction & sysNbInter)
{
  if (! sysNbInter.beBuilt()) {
    throw MDExcptUnbuiltNonBondedInteraction ("InteractionEngine_interface");
  }
  if (sysNbInter.numberOfInteraction() > MaxNumberBondedInteraction ){
    throw MDExcptExceedConstantMemLimit (
	"InteractionEngine::registNonBondedInteraction",
	"nonBonedInteractionType",
	MaxNumberNonBondedInteraction * sizeof(InteractionType));
  }
  if (sysNbInter.numberOfParameter() > MaxNumberNonBondedInteractionParameter ){
    throw MDExcptExceedConstantMemLimit (
	"InteractionEngine::registNonBondedInteraction",
	"nonBondedInteractionParameter",
	MaxNumberNonBondedInteractionParameter * sizeof(ScalorType));
  }

  cudaMemcpyToSymbol (nonBondedInteractionType,
		      sysNbInter.interactionType(), 
  		      sizeof(InteractionType) * sysNbInter.numberOfInteraction());
  cudaMemcpyToSymbol (nonBondedInteractionParameterPosition,
		      sysNbInter.interactionParameterPosition(),
  		      sizeof(ScalorType) * sysNbInter.numberOfInteraction());
  cudaMemcpyToSymbol (nonBondedInteractionParameter,
		      sysNbInter.interactionParameter(),
		      sizeof(IndexType) * sysNbInter.numberOfParameter());
  checkCUDAError ("InteractionEngine::init, init NB force setting");

  IndexType tableSize = sysNbInter.interactionTableSize();
  IndexType tmpNumAtomType = sysNbInter.numberOfAtomTypes();
  if (tableSize > MaxLengthNonBondedInteractionTable){
    throw MDExcptExceedConstantMemLimit(
	"InteractionEngine::registNonBondedInteraction",
	"nonBondedInteractionTable",
	MaxLengthNonBondedInteractionTable * sizeof (ScalorType));
  }
  cudaMemcpyToSymbol (const_nonBondedInteractionTableLength,
  		      &tableSize,
  		      sizeof (IndexType));
  checkCUDAError ("InteractionEngine::init, const_nonBondedInteractionTableLength");
  cudaMemcpyToSymbol (const_numAtomType,
		      &tmpNumAtomType,
		      sizeof (IndexType));
  checkCUDAError ("InteractionEngine::init, const_numAtomType");
  cudaMemcpyToSymbol (const_nonBondedInteractionTable,
  		      sysNbInter.interactionTable(),
  		      sizeof (IndexType) * tableSize);
  checkCUDAError ("InteractionEngine::init, const_nonBondedInteractionTable");

  IndexType numThreadsInCell = Parallel::Interface::numThreadsInCell();
  applyNonBondedInteraction_CellList_sbuffSize =
      sizeof(CoordType) * numThreadsInCell +
      sizeof(TypeType)  * numThreadsInCell;
  checkCUDAError ("InteractionEngine::init, init nonBondedInteractionTable");
}

void Parallel::InteractionEngine::
applyNonBondedInteraction (DeviceCellListedMDData & ddata,
			   const DeviceCellRelation & relation)
{
  Parallel::CudaGlobal::calNonBondedInteraction
      <<<gridDim, Parallel::Interface::numThreadsInCell()>>> (
	  ddata.dptr_coordinate(),
	  ddata.dptr_type(),
	  ddata.getGlobalBox().size,
	  ddata.getGlobalBox().sizei,
	  ddata.getRlist(),
	  ddata.dptr_numAtomInCell(),
	  relation.dptr_numNeighborCell(),
	  relation.dptr_neighborCellIndex(),
	  totalNumCell,
	  ddata.dptr_forceX(),
	  ddata.dptr_forceY(),
	  ddata.dptr_forceZ(),
	  sum_nb_p.getBuff(),
	  sum_nb_vxx.getBuff(),
	  sum_nb_vyy.getBuff(),
	  sum_nb_vzz.getBuff(),
	  err.ptr_de);
  checkCUDAError ("InteractionEngine::applyNonBondedInteraction");
}

void Parallel::InteractionEngine::
applyNonBondedInteraction (DeviceCellListedMDData & ddata,
			   const DeviceCellRelation & relation,
			   DeviceStatistic & st)
{
  Parallel::CudaGlobal::calNonBondedInteraction
      <<<gridDim, Parallel::Interface::numThreadsInCell()>>> (
	  ddata.dptr_coordinate(),
	  ddata.dptr_type(),
	  ddata.getGlobalBox().size,
	  ddata.getGlobalBox().sizei,
	  ddata.getRlist(),
	  ddata.dptr_numAtomInCell(),
	  relation.dptr_numNeighborCell(),
	  relation.dptr_neighborCellIndex(),
	  totalNumCell,
	  ddata.dptr_forceX(),
	  ddata.dptr_forceY(),
	  ddata.dptr_forceZ(),
	  sum_nb_p.getBuff(),
	  sum_nb_vxx.getBuff(),
	  sum_nb_vyy.getBuff(),
	  sum_nb_vzz.getBuff(),
	  err.ptr_de);
  checkCUDAError ("InteractionEngine::applyNonBondedInteraction");
  sum_nb_p.sumBuffAdd   (st.dptr_statisticData(), mdStatisticNonBondedPotential, 0);
  sum_nb_vxx.sumBuffAdd (st.dptr_statisticData(), mdStatisticVirialXX, 0);
  sum_nb_vyy.sumBuffAdd (st.dptr_statisticData(), mdStatisticVirialYY, 0);
  sum_nb_vzz.sumBuffAdd (st.dptr_statisticData(), mdStatisticVirialZZ, 0);
}
	  
__global__ void Parallel::CudaGlobal::
calNonBondedInteraction (const CoordType * coord,
			 const TypeType  * type,
			 const HostVectorType boxSize,
			 const HostVectorType boxSizei,
			 const ScalorType  rlist,
			 const IndexType * numAtomInCell,
			 const IndexType * numNeighborCell,
			 const IndexType * neighborCellIndex,
			 const IndexType   stride,
			 ScalorType * forcx,
			 ScalorType * forcy,
			 ScalorType * forcz,
			 ScalorType * statistic_nb_buff0,
			 ScalorType * statistic_nb_buff1,
			 ScalorType * statistic_nb_buff2,
			 ScalorType * statistic_nb_buff3,
			 mdError_t * ptr_de)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;

  IndexType this_numAtomInCell;
  IndexType this_numNeighborCell;  
  IndexType target_cellIndex;
  IndexType target_numAtomInCell;

  this_numAtomInCell = numAtomInCell[bid];
  if (this_numAtomInCell == 0) return;
  this_numNeighborCell = numNeighborCell[bid];

  IndexType ii = bid * blockDim.x + tid;
  CoordType refCoord;
  TypeType refType;
  ScalorType fsumx (0.f), fsumy(0.f), fsumz(0.f);
  ScalorType myPoten (0.0f), myVxx (0.0f), myVyy (0.0f), myVzz (0.0f);

  if (tid < this_numAtomInCell){
    refCoord = coord[ii];
    refType = type[ii];
  }  
  ScalorType rlist2 = rlist * rlist;
  
  extern __shared__ volatile char pub_sbuff[];
  CoordType * targetCoord =
      (CoordType *) & pub_sbuff;
  TypeType * targetType =
      (TypeType *) & targetCoord[blockDim.x];
  
  for (IndexType kk = 0; kk < this_numNeighborCell; ++kk){
    target_cellIndex = neighborCellIndex[bid * stride + kk];
    target_numAtomInCell = numAtomInCell[target_cellIndex];
    if (target_numAtomInCell == 0) continue;
    IndexType tmpLower = target_cellIndex * blockDim.x;
    IndexType jj = tmpLower + tid;
    if (tid < target_numAtomInCell) {
      targetCoord[tid] = coord[jj];
      targetType[tid] = type[jj];
    }
    if (tid < this_numAtomInCell){
      IndexType tmpUpper = tmpLower + target_numAtomInCell;
      for (IndexType ll = tmpLower; ll < tmpUpper; ++ll){
	if (ll != ii) {
	  ScalorType diffx = targetCoord[ll].x - refCoord.x;
	  ScalorType diffy = targetCoord[ll].y - refCoord.y;
	  ScalorType diffz = targetCoord[ll].z - refCoord.z;
	  shortestImage (boxSize.x, boxSizei.x, &diffx);
	  shortestImage (boxSize.y, boxSizei.y, &diffy);
	  shortestImage (boxSize.z, boxSizei.z, &diffz);
	  if (diffx*diffx+diffy*diffy+diffz*diffz < rlist2) {
	    IndexType fidx(0);
	    fidx = Parallel::CudaDevice::calNonBondedForceIndex (
		const_nonBondedInteractionTable,
		const_numAtomType[0],
		refType,
		targetType[ll]);
	    ScalorType fx, fy, fz, dp;
	    nbForcePoten (nonBondedInteractionType[fidx],
			  &nonBondedInteractionParameter
			  [nonBondedInteractionParameterPosition[fidx]],
			  diffx, diffy, diffz,
			  &fx, &fy, &fz, &dp);
	    myPoten += dp;
	    myVxx += fx * diffx;
	    myVyy += fy * diffy;
	    myVzz += fz * diffz;
	    fsumx += fx;
	    fsumy += fy;
	    fsumz += fz;
	  }
	  // __syncthreads();
	}
      }
    }
  }

  if (tid < this_numAtomInCell){
    forcx[ii] += fsumx;
    forcy[ii] += fsumy;
    forcz[ii] += fsumz;
  }
  
  statistic_nb_buff0[ii] = myPoten * 0.5f;
  statistic_nb_buff1[ii] = myVxx * 0.5f;
  statistic_nb_buff2[ii] = myVyy * 0.5f;
  statistic_nb_buff3[ii] = myVzz * 0.5f;
}

	  
    



    
