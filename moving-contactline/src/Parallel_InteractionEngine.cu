#define DEVICE_CODE

#include "common.h"
#include "Parallel_InteractionEngine.h"
#include "Parallel_Interface.h"
#include "NonBondedInteraction.h"
#include "Parallel_Auxiliary.h"
#include "BondInteraction.h"

#include "compile_error_mixcode.h"

texture<CoordType,  1, cudaReadModeElementType> global_texRef_interaction_coord;
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
__constant__
InteractionType bondedInteractionType [MaxNumberBondedInteraction];
__constant__
IndexType bondedInteractionParameterPosition [MaxNumberBondedInteraction];
__constant__
ScalorType bondedInteractionParameter [MaxNumberBondedInteractionParamemter];


Parallel::InteractionEngine::
InteractionEngine ()
    : inited(false), hasBond (false), hasAngle(false), rcut (0.f)
{
}

Parallel::InteractionEngine::
InteractionEngine (const DeviceCellListedMDData & ddata)
    : inited(false), hasBond (false), hasAngle(false), rcut (0.f)
{
  reinit (ddata);
}

void Parallel::InteractionEngine::
reinit (const DeviceCellListedMDData & ddata)
{
  clear ();
  
  totalNumCell = ddata.getNumCell().x * ddata.getNumCell().y * ddata.getNumCell().z;
  devideLevel = ddata.getDevideLevel();
  gridDim = toGridDim (totalNumCell);
  numThreadsInCell = Parallel::Interface::numThreadsInCell();
  
  // size_t sizetype = sizeof(TypeType)*sys.ddata.numMem;
  cudaBindTexture(0,
		  global_texRef_interaction_coord,
		  ddata.dptr_coordinate(),
		  sizeof(CoordType) * totalNumCell * numThreadsInCell);
		  
  checkCUDAError ("InteractionEngine::init, bind texture");
  
  sum_nb_p.reinit (totalNumCell, NThreadForSum);
  sum_nb_vxx.reinit (totalNumCell, NThreadForSum);
  sum_nb_vyy.reinit (totalNumCell, NThreadForSum);
  sum_nb_vzz.reinit (totalNumCell, NThreadForSum);
  sum_b_p.reinit (totalNumCell, NThreadForSum);
  sum_b_vxx.reinit (totalNumCell, NThreadForSum);
  sum_b_vyy.reinit (totalNumCell, NThreadForSum);
  sum_b_vzz.reinit (totalNumCell, NThreadForSum);
  sum_angle_p.reinit (totalNumCell, NThreadForSum);

  // prepare for clear ghost cell bond
  // tell the engine who are ghost cells
  SubCellList ghostCells;
  ddata.buildSubListGhostCell (ghostCells);
  numGhostCell = ghostCells.size();
  size_t size = numGhostCell * sizeof(IndexType);
  IndexType * host_ghostCellIndex = (IndexType *) malloc (size);
  if (host_ghostCellIndex == NULL){
    throw MDExcptFailedMallocOnHost ("InteractionEngine::reinit",
				     "host_ghostCellIndex", size);
  }
  for (IndexType i = 0; i < numGhostCell; ++i){
    host_ghostCellIndex[i] = ghostCells[i];
  }
  cudaMalloc ((void**)&ghostCellIndex, size);
  checkCUDAError ("InteractionEngine::reinit, malloc ghostCellIndex");
  cudaMemcpy (ghostCellIndex, host_ghostCellIndex, size, cudaMemcpyHostToDevice);
  checkCUDAError ("InteractionEngine::reinit, copy ghostCellIndex");
  free (host_ghostCellIndex);

  inited = true;
}

void Parallel::InteractionEngine::
clear ()
{
  if (inited){
    cudaFree (ghostCellIndex);
    cudaUnbindTexture(global_texRef_interaction_coord);
    checkCUDAError ("InteractionEngine::clear()");
    inited = false;
  }
}


Parallel::InteractionEngine::
~InteractionEngine ()
{
  clear ();
}

void Parallel::InteractionEngine::
registNonBondedInteraction (const SystemNonBondedInteraction & sysNbInter)
{
  rcut = sysNbInter.maxRcut();
  
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

  size_t tmpSize;
  sizeof(TypeType) > sizeof(ScalorType) ?
      tmpSize = sizeof(TypeType)   * numThreadsInCell :
      tmpSize = sizeof(ScalorType) * numThreadsInCell ;
  applyNonBondedInteraction_CellList_sbuffSize =
      sizeof(CoordType) * numThreadsInCell + tmpSize;
  checkCUDAError ("InteractionEngine::init, init nonBondedInteractionTable");
}


void Parallel::InteractionEngine::
registBondedInteraction (const SystemBondedInteraction & sysBdInter)
{
  if (sysBdInter.hasBond() ){
    hasBond = true;
  }
  if (sysBdInter.hasAngle()){
    hasAngle = true;
  }

  if (sysBdInter.numberOfInteraction() > MaxNumberBondedInteraction ){
    throw MDExcptExceedConstantMemLimit (
	"InteractionEngine::registBondedInteraction",
	"bondedInteractionType",
	MaxNumberBondedInteraction * sizeof(InteractionType));
  }
  if (sysBdInter.numberOfParameter() > MaxNumberBondedInteractionParamemter ){
    throw MDExcptExceedConstantMemLimit (
	"InteractionEngine::registBondedInteraction",
	"bondedInteractionParameter",
	MaxNumberBondedInteractionParamemter * sizeof(ScalorType));
  }

  if (hasBond || hasAngle){
    cudaMemcpyToSymbol (bondedInteractionType,
			sysBdInter.interactionType(),
			sizeof(InteractionType) * sysBdInter.numberOfInteraction());
    cudaMemcpyToSymbol (bondedInteractionParameterPosition,
			sysBdInter.interactionParameterPosition(),
			sizeof(ScalorType) * sysBdInter.numberOfInteraction());
    cudaMemcpyToSymbol (bondedInteractionParameter,
			sysBdInter.interactionParameter(),
			sizeof(IndexType) * sysBdInter.numberOfParameter());
    checkCUDAError ("InteractionEngine::init, init bond force setting");
    // cal shared buff size
    calBondInteraction_sbuffSize  = numThreadsInCell * sizeof(ScalorType);
    calAngleInteraction_sbuffSize = numThreadsInCell * sizeof(ScalorType);
  }
}  


void Parallel::InteractionEngine::
applyNonBondedInteraction (DeviceCellListedMDData & ddata,
			   const DeviceCellRelation & relation,
			   DeviceStatistic & st)
{
  Parallel::CudaGlobal::calNonBondedInteraction
      <<<gridDim, Parallel::Interface::numThreadsInCell(),
      applyNonBondedInteraction_CellList_sbuffSize>>> (
	  ddata.dptr_coordinate(),
	  ddata.dptr_type(),
	  ddata.getGlobalBox().size,
	  ddata.getGlobalBox().sizei,
	  ddata.getRlist(),
	  ddata.dptr_numAtomInCell(),
	  relation.dptr_numNeighborCell(),
	  relation.dptr_neighborCellIndex(),
	  relation.dptr_neighborShiftNoi(),
	  relation.stride_neighborCellIndex(),
	  ddata.dptr_forceX(),
	  ddata.dptr_forceY(),
	  ddata.dptr_forceZ(),
	  sum_nb_p.buff,
	  sum_nb_vxx.buff,
	  sum_nb_vyy.buff,
	  sum_nb_vzz.buff,
	  err.ptr_de);
  checkCUDAError ("InteractionEngine::applyNonBondedInteraction");
  sum_nb_p.sumBuffAdd   (st.dptr_statisticData(), mdStatistic_NonBondedPotential, 0);
  sum_nb_vxx.sumBuffAdd (st.dptr_statisticData(), mdStatistic_VirialXX, 0);
  sum_nb_vyy.sumBuffAdd (st.dptr_statisticData(), mdStatistic_VirialYY, 0);
  sum_nb_vzz.sumBuffAdd (st.dptr_statisticData(), mdStatistic_VirialZZ, 0);
}

void Parallel::InteractionEngine::
applyNonBondedInteraction (DeviceCellListedMDData & ddata,
			   const DeviceCellRelation & relation)
{
  Parallel::CudaGlobal::calNonBondedInteraction
      <<<gridDim, Parallel::Interface::numThreadsInCell(),
      applyNonBondedInteraction_CellList_sbuffSize>>> (
	  ddata.dptr_coordinate(),
	  ddata.dptr_type(),
	  ddata.getGlobalBox().size,
	  ddata.getGlobalBox().sizei,
	  ddata.getRlist(),
	  ddata.dptr_numAtomInCell(),
	  relation.dptr_numNeighborCell(),
	  relation.dptr_neighborCellIndex(),
	  relation.dptr_neighborShiftNoi(),
	  relation.stride_neighborCellIndex(),
	  ddata.dptr_forceX(),
	  ddata.dptr_forceY(),
	  ddata.dptr_forceZ(),
	  err.ptr_de);
  checkCUDAError ("InteractionEngine::applyNonBondedInteraction");
}



__device__ IndexType Parallel::CudaDevice::
calNonBondedForceIndex (const IndexType * table, 
			const IndexType numType,
			const TypeType atom0, 
			const TypeType atom1)
{
  IndexType i, j;
  atom0 <= atom1 ?
      (i = atom0, j = atom1) :
      (i = atom1, j = atom0) ;
  return table[i * numType + j - ((i*(i+1)) >> 1)];
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
			 const CoordNoiType * neighborShiftNoi,
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
  IndexType ii = bid * blockDim.x + tid;

  this_numNeighborCell = numNeighborCell[bid];
  if (this_numNeighborCell == 0) {
    if (threadIdx.x == 0){
      statistic_nb_buff0[bid] = 0.f;
      statistic_nb_buff1[bid] = 0.f;
      statistic_nb_buff2[bid] = 0.f;
      statistic_nb_buff3[bid] = 0.f;
    }
    return;
  }  
  this_numAtomInCell = numAtomInCell[bid];
  if (this_numAtomInCell == 0) {
    if (threadIdx.x == 0){
      statistic_nb_buff0[bid] = 0.f;
      statistic_nb_buff1[bid] = 0.f;
      statistic_nb_buff2[bid] = 0.f;
      statistic_nb_buff3[bid] = 0.f;
    }
    return;
  }  
    
  // if (tid == 0){
  //   printf ("bid: %d, numNei: %d\n", bid, this_numNeighborCell);
  // }    

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
      (CoordType *) & pub_sbuff[0];
  TypeType * targetType =
      (TypeType *) & targetCoord[blockDim.x];
  ScalorType * st_buff =
      (ScalorType *) & targetCoord[blockDim.x];
  
  // IndexType count = 0;
  for (IndexType kk = 0; kk < this_numNeighborCell; ++kk){
    __syncthreads();
    IndexType target_cellIndex;
    CoordType target_shift;
    IndexType target_numAtomInCell;
    target_cellIndex = neighborCellIndex[bid * stride + kk];
    target_shift.x   = neighborShiftNoi [bid * stride + kk].x * boxSize.x;
    target_shift.y   = neighborShiftNoi [bid * stride + kk].y * boxSize.y;
    target_shift.z   = neighborShiftNoi [bid * stride + kk].z * boxSize.z;
    target_numAtomInCell = numAtomInCell[target_cellIndex];
    if (target_numAtomInCell == 0) continue;
    IndexType indexShift = target_cellIndex * blockDim.x;
    IndexType jj = indexShift + tid;
    if (tid < target_numAtomInCell) {
      targetCoord[tid] = coord[jj];
      targetType[tid] = type[jj];
    }
    __syncthreads();
    // if (ii == 640*4){
    //   printf ("%f %f %f\n", targetCoord[0].x, targetCoord[0].y, targetCoord[0].z);
    // }    
    target_shift.x -= refCoord.x;
    target_shift.y -= refCoord.y;
    target_shift.z -= refCoord.z;
    if (tid < this_numAtomInCell){
      for (IndexType ll = 0; ll < target_numAtomInCell; ++ll){
	ScalorType diffx = targetCoord[ll].x + target_shift.x;
	ScalorType diffy = targetCoord[ll].y + target_shift.y;
	ScalorType diffz = targetCoord[ll].z + target_shift.z;
	if (diffx*diffx+diffy*diffy+diffz*diffz < rlist2 && ll + indexShift != ii) {
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
	  // printf ("%f, %f %f %f,  %f %f %f,  %f %f %f,  %f\n",
	  // 	  sqrtf(diffx*diffx+diffy*diffy+diffz*diffz),
	  // 	  refCoord.x, refCoord.y, refCoord.z,
	  // 	  targetCoord[ll].x, targetCoord[ll].y, targetCoord[ll].z,
	  // 	  diffx, diffy, diffz,
	  // 	  dp);
	  myPoten += dp;
	  myVxx += fx * diffx;
	  myVyy += fy * diffy;
	  myVzz += fz * diffz;
	  fsumx += fx;
	  fsumy += fy;
	  fsumz += fz;
	}
      }
    }
  }
  
  if (tid < this_numAtomInCell){
    forcx[ii] += fsumx;
    forcy[ii] += fsumy;
    forcz[ii] += fsumz;
  }

  __syncthreads();
  if (threadIdx.x < this_numAtomInCell){
    st_buff[threadIdx.x] = 0.5f * myPoten;
  }
  else {
    st_buff[threadIdx.x] = 0.f;
  }
  __syncthreads();
  sumVectorBlockBuffer_2 (st_buff);
  if (threadIdx.x == 0) statistic_nb_buff0[bid] = st_buff[0];
  __syncthreads();

  if (threadIdx.x < this_numAtomInCell){
    st_buff[threadIdx.x] = 0.5f * myVxx;
  }
  else {
    st_buff[threadIdx.x] = 0.f;
  }
  __syncthreads();
  sumVectorBlockBuffer_2 (st_buff);
  if (threadIdx.x == 0) statistic_nb_buff1[bid] = st_buff[0];
  __syncthreads();

  if (threadIdx.x < this_numAtomInCell){
    st_buff[threadIdx.x] = 0.5f * myVyy;
  }
  else {
    st_buff[threadIdx.x] = 0.f;
  }
  __syncthreads();
  sumVectorBlockBuffer_2 (st_buff);
  if (threadIdx.x == 0) statistic_nb_buff2[bid] = st_buff[0];
  __syncthreads();

  if (threadIdx.x < this_numAtomInCell){
    st_buff[threadIdx.x] = 0.5f * myVzz;
  }
  else {
    st_buff[threadIdx.x] = 0.f;
  }
  __syncthreads();
  sumVectorBlockBuffer_2 (st_buff);
  if (threadIdx.x == 0) statistic_nb_buff3[bid] = st_buff[0];
  
  // statistic_nb_buff0[ii] = myPoten * 0.5f;
  // statistic_nb_buff1[ii] = myVxx * 0.5f;
  // statistic_nb_buff2[ii] = myVyy * 0.5f;
  // statistic_nb_buff3[ii] = myVzz * 0.5f;
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
			 const CoordNoiType * neighborShiftNoi,
			 const IndexType   stride,
			 ScalorType * forcx,
			 ScalorType * forcy,
			 ScalorType * forcz,
			 mdError_t * ptr_de)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;

  IndexType this_numAtomInCell;
  IndexType this_numNeighborCell;  
  IndexType ii = bid * blockDim.x + tid;

  this_numNeighborCell = numNeighborCell[bid];
  if (this_numNeighborCell == 0) {
    return;
  }  
  this_numAtomInCell = numAtomInCell[bid];
  if (this_numAtomInCell == 0) {
    return;
  }  
    
  // if (tid == 0){
  //   printf ("bid: %d, numNei: %d\n", bid, this_numNeighborCell);
  // }    

  CoordType refCoord;
  TypeType refType;
  ScalorType fsumx (0.f), fsumy(0.f), fsumz(0.f);

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
  
  // IndexType count = 0;
  for (IndexType kk = 0; kk < this_numNeighborCell; ++kk){
    __syncthreads();
    IndexType target_cellIndex;
    CoordType target_shift;
    IndexType target_numAtomInCell;
    target_cellIndex = neighborCellIndex[bid * stride + kk];
    target_shift.x   = neighborShiftNoi [bid * stride + kk].x * boxSize.x;
    target_shift.y   = neighborShiftNoi [bid * stride + kk].y * boxSize.y;
    target_shift.z   = neighborShiftNoi [bid * stride + kk].z * boxSize.z;
    target_numAtomInCell = numAtomInCell[target_cellIndex];
    if (target_numAtomInCell == 0) continue;
    IndexType indexShift = target_cellIndex * blockDim.x;
    IndexType jj = indexShift + tid;
    if (tid < target_numAtomInCell) {
      targetCoord[tid] = coord[jj];
      targetType[tid] = type[jj];
    }
    __syncthreads();
    target_shift.x -= refCoord.x;
    target_shift.y -= refCoord.y;
    target_shift.z -= refCoord.z;
    if (tid < this_numAtomInCell){
      for (IndexType ll = 0; ll < target_numAtomInCell; ++ll){
	ScalorType diffx = targetCoord[ll].x + target_shift.x;
	ScalorType diffy = targetCoord[ll].y + target_shift.y;
	ScalorType diffz = targetCoord[ll].z + target_shift.z;
	if ((diffx*diffx+diffy*diffy+diffz*diffz) < rlist2 && ll + indexShift != ii) {
	  IndexType fidx;
	  fidx = Parallel::CudaDevice::calNonBondedForceIndex (
	      const_nonBondedInteractionTable,
	      const_numAtomType[0],
	      refType,
	      targetType[ll]);
	  ScalorType fx, fy, fz;
	  nbForce (nonBondedInteractionType[fidx],
		   &nonBondedInteractionParameter
		   [nonBondedInteractionParameterPosition[fidx]],
		   diffx,
		   diffy,
		   diffz,
		   &fx, &fy, &fz);
	  // printf ("%f, %f %f %f,  %f %f %f,  %f %f %f\n",
	  // 	  sqrtf(diffx*diffx+diffy*diffy+diffz*diffz),
	  // 	  refCoord.x, refCoord.y, refCoord.z,
	  // 	  targetCoord[ll].x, targetCoord[ll].y, targetCoord[ll].z,
	  // 	  diffx, diffy, diffz);
	  fsumx += fx;
	  fsumy += fy;
	  fsumz += fz;
	}
      }
    }
  }
  // printf ("bid: %d, tid: %d, num eff: %d. fsum %f\n", bid, tid, count, fsumx);

  if (tid < this_numAtomInCell){
    forcx[ii] += fsumx;
    forcy[ii] += fsumy;
    forcz[ii] += fsumz;
  }  
}


void Parallel::InteractionEngine::
clearInteraction (DeviceCellListedMDData & data)
{
  Parallel::Auxiliary::setValue
      <<<gridDim, Parallel::Interface::numThreadsInCell()>>> (
	  data.dptr_forceX(), 0.f);
  Parallel::Auxiliary::setValue
      <<<gridDim, Parallel::Interface::numThreadsInCell()>>> (
	  data.dptr_forceY(), 0.f);
  Parallel::Auxiliary::setValue
      <<<gridDim, Parallel::Interface::numThreadsInCell()>>> (
	  data.dptr_forceZ(), 0.f);
  checkCUDAError ("InteractionEngine::clearInteraction");
}
    
void Parallel::InteractionEngine::
applyBondedInteraction (DeviceCellListedMDData & ddata,
			DeviceBondList & dbdlist)
{
  if (hasBond){
    Parallel::CudaGlobal::clearGhostBond
	<<<numGhostCell, numThreadsInCell>>> (
	    ghostCellIndex,
	    ddata.dptr_numBond());
    Parallel::CudaGlobal::calBondInteraction
	<<<gridDim, numThreadsInCell>>>(
	    ddata.dptr_coordinate(),
	    ddata.getGlobalBox().size,
	    ddata.getGlobalBox().sizei,
	    ddata.dptr_numAtomInCell(),
	    ddata.dptr_numBond(),
	    ddata.dptr_bondIndex(),
	    ddata.bondTopStride(),
	    dbdlist.dptr_bondNeighbor_localIndex(),
	    dbdlist.getStride(),
	    ddata.dptr_forceX(),
	    ddata.dptr_forceY(),
	    ddata.dptr_forceZ(),
	    err.ptr_de,
	    err.ptr_dscalor);
    checkCUDAError ("InteractionEngine::applyBondedInteraction");
    err.check ("InteractionEngine::applyBondedInteraction no st");
  }
}


void Parallel::InteractionEngine::
applyBondedInteraction (DeviceCellListedMDData & ddata,
			DeviceBondList & dbdlist,
			DeviceStatistic & st)
{
  if (hasBond){
    Parallel::CudaGlobal::clearGhostBond
	<<<numGhostCell, numThreadsInCell>>> (
	    ghostCellIndex,
	    ddata.dptr_numBond());
    Parallel::CudaGlobal::calBondInteraction
	<<<gridDim, numThreadsInCell, calBondInteraction_sbuffSize>>>(
	    ddata.dptr_coordinate(),
	    ddata.getGlobalBox().size,
	    ddata.getGlobalBox().sizei,
	    ddata.dptr_numAtomInCell(),
	    ddata.dptr_numBond(),
	    ddata.dptr_bondIndex(),
	    ddata.bondTopStride(),
	    dbdlist.dptr_bondNeighbor_localIndex(),
	    dbdlist.getStride(),
	    ddata.dptr_forceX(),
	    ddata.dptr_forceY(),
	    ddata.dptr_forceZ(),
	    sum_b_p.buff,
	    sum_b_vxx.buff,
	    sum_b_vyy.buff,
	    sum_b_vzz.buff,
	    err.ptr_de,
	    err.ptr_dscalor);
    checkCUDAError ("InteractionEngine::applyBondedInteraction");
    err.check ("InteractionEngine::applyBondedInteraction with st");
    sum_b_p.sumBuffAdd   (st.dptr_statisticData(), mdStatistic_BondedPotential, 0);
    sum_b_vxx.sumBuffAdd (st.dptr_statisticData(), mdStatistic_VirialXX, 0);
    sum_b_vyy.sumBuffAdd (st.dptr_statisticData(), mdStatistic_VirialYY, 0);
    sum_b_vzz.sumBuffAdd (st.dptr_statisticData(), mdStatistic_VirialZZ, 0);
  }
}

__global__ void Parallel::CudaGlobal::
clearGhostBond (const IndexType * ghostCellIndex,
		IndexType * myNumBond)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType this_cellIndex = ghostCellIndex[bid];
  IndexType ii = this_cellIndex * blockDim.x + tid;
  myNumBond[ii] = 0;
}

__global__ void Parallel::CudaGlobal::
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
		    ScalorType * errsrc)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;
  
  IndexType my_numBond = numBond[ii];
  if (my_numBond == 0) return;
  // if (__all(my_numBond == 0)) return;
  IndexType this_numAtomInCell = numAtomInCell[bid];
  if (tid >= this_numAtomInCell) return;
  
  ScalorType fsumx = 0.0f;
  ScalorType fsumy = 0.0f;
  ScalorType fsumz = 0.0f;
  CoordType ref = coord[ii];

  for (IndexType jj = 0; jj < my_numBond; ++jj){
    IndexType list_index = bondListStride * jj + ii;
    IndexType top_index = bondTopStride * jj + ii;
    IndexType target_index = bondNeighbor_localIndex[list_index];
    IndexType my_bondIndex   = bondIndex[top_index];
    CoordType target_coord = tex1Dfetch (global_texRef_interaction_coord, target_index);
    ScalorType diffx, diffy, diffz;
    ScalorType fx, fy, fz;
    diffx = target_coord.x - ref.x;
    diffy = target_coord.y - ref.y;
    diffz = target_coord.z - ref.z;
    shortestImage (boxSize.x, boxSizei.x, &diffx);
    shortestImage (boxSize.y, boxSizei.y, &diffy);
    shortestImage (boxSize.z, boxSizei.z, &diffz);
    if (bondedInteractionType[my_bondIndex] == mdForceFENE){
      ScalorType rinf2 = bondedInteractionParameter
    	  [bondedInteractionParameterPosition[my_bondIndex] + 1];
      ScalorType diff2 = diffx*diffx + diffy*diffy + diffz*diffz;
      if (diff2 > rinf2){
    	*ptr_de = mdErrorBreakFENEBond;
    	errsrc[0] = target_coord.x;
    	errsrc[1] = target_coord.y;
    	errsrc[2] = target_coord.z;
    	errsrc[3] = ref.x;
    	errsrc[4] = ref.y;
    	errsrc[5] = ref.z;
    	errsrc[6] = diffx;
    	errsrc[7] = diffy;
    	errsrc[8] = diffz;
    	continue;
      }
    }	  
    bondForce (bondedInteractionType[my_bondIndex],
	       &bondedInteractionParameter
	       [bondedInteractionParameterPosition[my_bondIndex]],
	       diffx, diffy, diffz, &fx, &fy, &fz);
    fsumx += fx;
    fsumy += fy;
    fsumz += fz;
  }

  forcx[ii] += fsumx;
  forcy[ii] += fsumy;
  forcz[ii] += fsumz;  
}


__global__ void Parallel::CudaGlobal::
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
		    ScalorType * errsrc)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;
  
  IndexType my_numBond = numBond[ii];
  // if (__all(my_numBond == 0)) return;
  IndexType this_numAtomInCell = numAtomInCell[bid];

  extern __shared__ volatile ScalorType buff[];
  buff[tid] = 0.f;
  __syncthreads();
  
  ScalorType fsumx = 0.0f;
  ScalorType fsumy = 0.0f;
  ScalorType fsumz = 0.0f;
  ScalorType myPoten = 0.0f, myVxx = 0.0f, myVyy = 0.0f, myVzz = 0.0f;
  
  if (tid < this_numAtomInCell){
    CoordType ref = coord[ii];
    for (IndexType jj = 0; jj < my_numBond; ++jj){
      IndexType list_index = bondListStride * jj + ii;
      IndexType top_index  = bondTopStride  * jj + ii;
      IndexType target_index = bondNeighbor_localIndex[list_index];
      IndexType my_bondIndex = bondIndex[top_index];
      CoordType target_coord = tex1Dfetch (global_texRef_interaction_coord, target_index);
      ScalorType diffx, diffy, diffz;
      ScalorType fx, fy, fz, dp;
      diffx = target_coord.x - ref.x;
      diffy = target_coord.y - ref.y;
      diffz = target_coord.z - ref.z;
      // printf ("%f %f %f    %f %f %f\n",
      // 	      target_coord.x, target_coord.y, target_coord.z,
      // 	      ref.x, ref.y, ref.z);
      shortestImage (boxSize.x, boxSizei.x, &diffx);
      shortestImage (boxSize.y, boxSizei.y, &diffy);
      shortestImage (boxSize.z, boxSizei.z, &diffz);
      if (bondedInteractionType[my_bondIndex] == mdForceFENE){
	ScalorType rinf2 = bondedInteractionParameter
	    [bondedInteractionParameterPosition[my_bondIndex] + 1];
	ScalorType diff2 = diffx*diffx + diffy*diffy + diffz*diffz;
	if (diff2 > rinf2){
	  *ptr_de = mdErrorBreakFENEBond;
	  errsrc[0] = target_coord.x;
	  errsrc[1] = target_coord.y;
	  errsrc[2] = target_coord.z;
	  errsrc[3] = ref.x;
	  errsrc[4] = ref.y;
	  errsrc[5] = ref.z;
	  errsrc[6] = diffx;
	  errsrc[7] = diffy;
	  errsrc[8] = diffz;
	  errsrc[9] = my_numBond;
	  errsrc[10] = this_numAtomInCell;
	  errsrc[11] = bid;
	  errsrc[12] = tid;
	  continue;
	}
      }	  
      bondForcePoten (bondedInteractionType[my_bondIndex],
		      &bondedInteractionParameter
		      [bondedInteractionParameterPosition[my_bondIndex]],
		      diffx, diffy, diffz, &fx, &fy, &fz, &dp);
      // printf ("%f\t%f\t%f\n", diffx, sqrtf(diffx*diffx+diffy*diffy+diffz*diffz), dp);
      myPoten += dp;
      myVxx += fx * diffx;
      myVyy += fy * diffy;
      myVzz += fz * diffz;
      fsumx += fx;
      fsumy += fy;
      fsumz += fz;
    }

    forcx[ii] += fsumx;
    forcy[ii] += fsumy;
    forcz[ii] += fsumz;
  }

  __syncthreads();
  buff[tid] = myPoten * 0.5f;
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) statistic_b_buff0[bid] = buff[0];
  __syncthreads();
  buff[tid] = myVxx * 0.5f;
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) statistic_b_buff1[bid] = buff[0];
  __syncthreads();
  buff[tid] = myVyy * 0.5f;
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) statistic_b_buff2[bid] = buff[0];
  __syncthreads();
  buff[tid] = myVzz * 0.5f;
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) statistic_b_buff3[bid] = buff[0];
  __syncthreads();
}



    
