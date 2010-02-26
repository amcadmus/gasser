#include "InteractionEngine_interface.h"
#include "NonBondedInteraction.h"
#include "BondInteraction.h"
#include "AngleInteraction.h"
// #include "CellList_interface.h"

texture<CoordType,  1, cudaReadModeElementType> global_texRef_interaction_coord;
texture<TypeType ,  1, cudaReadModeElementType> global_texRef_interaction_type;
__constant__
InteractionType nonBondedInteractionType [MaxNumberNonBondedInteraction];
__constant__
ScalorType nonBondedInteractionParameter [MaxNumberNonBondedInteractionParameter];
__constant__
IndexType nonBondedInteractionParameterPosition [MaxNumberNonBondedInteraction];
__constant__
InteractionType bondedInteractionType [MaxNumberBondedInteraction];
__constant__
IndexType bondedInteractionParameterPosition [MaxNumberBondedInteraction];
__constant__
ScalorType bondedInteractionParameter [MaxNumberBondedInteractionParamemter];
// __constant__
// IndexType nonBondedInteractionTableLength;
__constant__
IndexType numAtomType[1];
__constant__
IndexType nonBondedInteractionTable [MaxLengthNonBondedInteractionTable];


void InteractionEngine_interface::init (const MDSystem  & sys,
					const IndexType & NTread)
{
  hasBond = false;
  hasAngle = false;
  myBlockDim.y = 1;
  myBlockDim.z = 1;
  myBlockDim.x = NTread;
  IndexType nob;
  if (sys.ddata.numAtom % myBlockDim.x == 0){
    nob = sys.ddata.numAtom / myBlockDim.x;
  } else {
    nob = sys.ddata.numAtom / myBlockDim.x + 1;
  }
  atomGridDim = toGridDim (nob);

  // size_t sizetype = sizeof(TypeType)*sys.ddata.numMem;
  cudaBindTexture(0, global_texRef_interaction_coord, sys.ddata.coord,
		  sizeof(CoordType) * sys.ddata.numMem);
  cudaBindTexture(0, global_texRef_interaction_type, sys.ddata.type,
		  sizeof(TypeType)  * sys.ddata.numMem);
  checkCUDAError ("InteractionEngine::init, bind texture");
  
  // init sum vectors
  sum_nb_p.init (sys.ddata.numAtom, NThreadForSum);
  sum_nb_vxx.init (sys.ddata.numAtom, NThreadForSum);
  sum_nb_vyy.init (sys.ddata.numAtom, NThreadForSum);
  sum_nb_vzz.init (sys.ddata.numAtom, NThreadForSum);
  sum_b_p.init (nob, NThreadForSum);
  sum_b_vxx.init (nob, NThreadForSum);
  sum_b_vyy.init (nob, NThreadForSum);
  sum_b_vzz.init (nob, NThreadForSum);
  sum_angle_p.init (nob, NThreadForSum);
  for (IndexType i = 0; i < 8; ++i){
    cudaStreamCreate(&sum_stream[i]);
  }
  checkCUDAError ("InteractionEngine::init init sum statistic");
  
  // init nb force param
}


static IndexType hroundUp4 (IndexType x)
{
  if (x & 3 == 0){
    return x;
  }
  else {
    return ((x >> 2) + 1) << 2;
  }
}

void InteractionEngine_interface::
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
  // cudaMemcpyToSymbol (&nonBondedInteractionTableLength,
  // 		      &tableSize,
  // 		      sizeof (IndexType));
  cudaMemcpyToSymbol (numAtomType,
		      &tmpNumAtomType,
		      sizeof (IndexType));
  cudaMemcpyToSymbol (nonBondedInteractionTable,
		      sysNbInter.interactionTable(),
		      sizeof (IndexType) * tableSize);
  checkCUDAError ("InteractionEngine::init, init nonBondedInteractionTable");

  applyNonBondedInteraction_CellList_sbuffSize =
	sizeof(IndexType) *      hroundUp4(myBlockDim.x) +
	sizeof(CoordType) * hroundUp4(myBlockDim.x) +
	sizeof(TypeType)  *      hroundUp4(myBlockDim.x);

}


void InteractionEngine_interface::
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
    calBondInteraction_sbuffSize  = myBlockDim.x * sizeof(ScalorType);
    calAngleInteraction_sbuffSize = myBlockDim.x * sizeof(ScalorType);
  }
}

InteractionEngine_interface::~InteractionEngine_interface()
{
  cudaUnbindTexture(global_texRef_interaction_coord);
  for (IndexType i = 0; i < 8; ++i){
    cudaStreamDestroy(sum_stream[i]);
  }
}

void InteractionEngine_interface::clearInteraction (MDSystem & sys)
{
  clearForce
      <<<atomGridDim, myBlockDim>>>(
	  sys.ddata.numAtom,
	  sys.ddata.forcx, sys.ddata.forcy, sys.ddata.forcz);
  checkCUDAError ("InteractionEngine::clearInteraction");
}


void InteractionEngine_interface::
applyNonBondedInteraction  (MDSystem & sys,
			    const NeighborList & nlist,
			    MDTimer *timer )
{
  if (timer != NULL) timer->tic(mdTimeNonBondedInteraction);
  calNonBondedInteraction
      <<<atomGridDim, myBlockDim>>> (
	  sys.ddata.numAtom,
	  sys.ddata.coord,
	  sys.ddata.forcx,  sys.ddata.forcy,  sys.ddata.forcz,
	  sys.ddata.type, 
	  sys.box, nlist.dnlist,
	  err.ptr_de,
	  err.ptr_dindex,
	  err.ptr_dscalor);
  checkCUDAError ("InteractionEngine::applyInteraction nb");
  err.check ("interaction engine nb");	
  if (timer != NULL) timer->toc(mdTimeNonBondedInteraction);
}

void InteractionEngine_interface::
applyNonBondedInteractionCell  (MDSystem & sys,
				const NeighborList & nlist,
				MDTimer *timer )
{
  if (timer != NULL) timer->tic(mdTimeNonBondedInteraction);
  calNonBondedInteraction
      <<<nlist.cellGridDim, nlist.myBlockDim,
      applyNonBondedInteraction_CellList_sbuffSize>>> (
	  sys.ddata.numAtom,
	  sys.ddata.coord,
	  sys.ddata.forcx,
	  sys.ddata.forcy,
	  sys.ddata.forcz,
	  sys.ddata.type, 
	  sys.box,
	  nlist.dclist,
	  err.ptr_de);
  checkCUDAError ("InteractionEngine::applyInteraction nb");
  err.check ("interaction engine nb");	
  if (timer != NULL) timer->toc(mdTimeNonBondedInteraction);
}


void InteractionEngine_interface::
applyBondedInteraction (MDSystem & sys,
			const BondedInteractionList & bdlist,
			MDTimer *timer )
{
  if (hasBond) {
    if (timer != NULL) timer->tic(mdTimeBondedInteraction);
    calBondInteraction
	<<<atomGridDim, myBlockDim>>> (
	    sys.ddata.numAtom,
	    sys.ddata.coord,
	    sys.ddata.forcx,  sys.ddata.forcy,  sys.ddata.forcz,
	    sys.box,
	    bdlist.deviceBondList());
    checkCUDAError ("InteractionEngine::applyInteraction bonded");
    err.check ("interaction engine b");	
    if (timer != NULL) timer->toc(mdTimeBondedInteraction);
  }
  if (hasAngle){
    if (timer != NULL) timer->tic(mdTimeAngleInteraction);
    calAngleInteraction
	<<<atomGridDim, myBlockDim>>> (
	    sys.ddata.numAtom,
	    sys.ddata.coord,
	    sys.ddata.forcx,  sys.ddata.forcy,  sys.ddata.forcz,
	    sys.box,
	    bdlist.deviceAngleList());
    checkCUDAError ("InteractionEngine::applyInteraction angle");
    err.check ("interaction engine angle");	
    if (timer != NULL) timer->toc(mdTimeAngleInteraction);
  }
}
  
void InteractionEngine_interface::
applyNonBondedInteraction (MDSystem & sys,
			   const NeighborList & nlist,
			   MDStatistic & st,
			   MDTimer *timer )
{
  if (timer != NULL) timer->tic(mdTimeNBInterStatistic);
  calNonBondedInteraction
      <<<atomGridDim, myBlockDim>>> (
	  sys.ddata.numAtom,
	  sys.ddata.coord,
	  sys.ddata.forcx,  sys.ddata.forcy,  sys.ddata.forcz,
	  sys.ddata.type, 
	  sys.box, nlist.dnlist
	  ,
	  sum_nb_p.getBuff(),
	  sum_nb_vxx.getBuff(),
	  sum_nb_vyy.getBuff(),
	  sum_nb_vzz.getBuff(),
	  err.ptr_de,
	  err.ptr_dindex,
	  err.ptr_dscalor
	  );
  checkCUDAError ("InteractionEngine::applyInteraction nb (with statistic)");
  err.check ("interaction engine nb");	
  cudaThreadSynchronize();
  sum_nb_p.sumBuffAdd(st.ddata, mdStatisticNonBondedPotential, 0);
  sum_nb_vxx.sumBuffAdd(st.ddata, mdStatisticVirialXX, 1);
  sum_nb_vyy.sumBuffAdd(st.ddata, mdStatisticVirialYY, 2);
  sum_nb_vzz.sumBuffAdd(st.ddata, mdStatisticVirialZZ, 3);
  cudaThreadSynchronize();
  if (timer != NULL) timer->toc(mdTimeNBInterStatistic);
}


void InteractionEngine_interface::
applyNonBondedInteractionCell (MDSystem & sys,
			       const NeighborList & nlist,
			       MDStatistic & st,
			       MDTimer *timer )
{
  if (timer != NULL) timer->tic(mdTimeNBInterStatistic);
  calNonBondedInteraction
      <<<nlist.cellGridDim, nlist.myBlockDim,
      applyNonBondedInteraction_CellList_sbuffSize>>> (
	  sys.ddata.numAtom,
	  sys.ddata.coord,
	  sys.ddata.forcx,
	  sys.ddata.forcy,
	  sys.ddata.forcz,
	  sys.ddata.type, 
	  sys.box,
	  nlist.dclist
	  ,
	  sum_nb_p.getBuff(),
	  sum_nb_vxx.getBuff(),
	  sum_nb_vyy.getBuff(),
	  sum_nb_vzz.getBuff(),
	  err.ptr_de
	  );
  checkCUDAError ("InteractionEngine::applyInteraction nb (with statistic)");
  err.check ("interaction engine nb");	
  cudaThreadSynchronize();
  sum_nb_p.sumBuffAdd(st.ddata, mdStatisticNonBondedPotential, 0);
  sum_nb_vxx.sumBuffAdd(st.ddata, mdStatisticVirialXX, 1);
  sum_nb_vyy.sumBuffAdd(st.ddata, mdStatisticVirialYY, 2);
  sum_nb_vzz.sumBuffAdd(st.ddata, mdStatisticVirialZZ, 3);
  cudaThreadSynchronize();
  if (timer != NULL) timer->toc(mdTimeNBInterStatistic);
}


void InteractionEngine_interface::
applyBondedInteraction (MDSystem & sys,
			const BondedInteractionList & bdlist,
			MDStatistic & st,
			MDTimer *timer)
{
  if (hasBond) {
    if (timer != NULL) timer->tic(mdTimeBInterStatistic);
    calBondInteraction
	<<<atomGridDim, myBlockDim,
	calBondInteraction_sbuffSize>>> (
	    sys.ddata.numAtom,
	    sys.ddata.coord,
	    sys.ddata.forcx,  sys.ddata.forcy,  sys.ddata.forcz,
	    sys.box,
	    bdlist.deviceBondList()
	    ,
	    sum_b_p.getBuff(),
	    sum_b_vxx.getBuff(),
	    sum_b_vyy.getBuff(),
	    sum_b_vzz.getBuff(),
	    err.ptr_de
	    );
    checkCUDAError ("InteractionEngine::applyInteraction bonded (with statistic)");
    err.check ("interaction engine");	
    if (timer != NULL) timer->toc(mdTimeBInterStatistic);
  }
  if (hasAngle){
    if (timer != NULL) timer->tic(mdTimeAngleInterStatistic);
    calAngleInteraction
	<<<atomGridDim, myBlockDim,
	calAngleInteraction_sbuffSize>>> (
	    sys.ddata.numAtom,
	    sys.ddata.coord,
	    sys.ddata.forcx,  sys.ddata.forcy,  sys.ddata.forcz,
	    sys.box,
	    bdlist.deviceAngleList(),
	    sum_angle_p.getBuff(),
	    err.ptr_de);
    checkCUDAError ("InteractionEngine::applyInteraction angle");
    err.check ("interaction engine angle");	
    if (timer != NULL) timer->toc(mdTimeAngleInterStatistic);
  }
  if (hasBond) {
    if (timer != NULL) timer->tic(mdTimeBInterStatistic);
    cudaThreadSynchronize();
    sum_b_p.sumBuffAdd(st.ddata, mdStatisticBondedPotential, 4);
    sum_b_vxx.sumBuffAdd(st.ddata, mdStatisticVirialXX, 5);
    sum_b_vyy.sumBuffAdd(st.ddata, mdStatisticVirialYY, 6);
    sum_b_vzz.sumBuffAdd(st.ddata, mdStatisticVirialZZ, 7);
    cudaThreadSynchronize();
    if (timer != NULL) timer->toc(mdTimeBInterStatistic);
  }
  if (hasAngle){
    if (timer != NULL) timer->tic(mdTimeAngleInterStatistic);
    sum_angle_p.sumBuffAdd(st.ddata, mdStatisticBondedPotential, 4);
    if (timer != NULL) timer->toc(mdTimeAngleInterStatistic);
  }
  checkCUDAError ("InteractionEngine::applyInteraction sum statistic (with statistic)");
}



__global__ void clearForce (const IndexType numAtom,
			    ScalorType * forcx,
			    ScalorType * forcy, 
			    ScalorType * forcz)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  if (ii < numAtom) {
    forcx[ii] = 0.0f;
    forcy[ii] = 0.0f;
    forcz[ii] = 0.0f;
  }
}


// __global__ void
// calNonBondedInteraction (const CoordType * coord,
// 			 const TypeType * type,
// 			 DeviceCellListData clist,
// 			 DeviceCellListProperty clistPro,
// 			 ScalorType * forcx,
// 			 ScalorType * forcy,
// 			 ScalorType * forcz,
// 			 bool sharednbForceTable)
// {
//   IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
//   IndexType tid = threadIdx.x;
//   ScalorType fsumx(0.f), fsumy(0.f), fsumz(0.f);
  
//   extern __shared__ volatile char pub_sbuff[];
//   volatile IndexType * targetIndex = (volatile IndexType *) pub_sbuff;
//   CoordType * targetCoord = (CoordType *) &targetIndex[roundUp4(blockDim.x)];
//   volatile TypeType * targetType = (volatile TypeType *) &targetCoord[roundUp4(blockDim.x)];
//   __syncthreads();

//   IndexType ii = get (clist, bid, tid);
//   CoordType ref;
//   TypeType refType;
//   if (ii != MaxIndexValue){
//     ref = tex1Dfetch (global_texRef_interaction_coord, ii);
//     refType = tex1Dfetch(global_texRef_interaction_type, ii);
//   }
  
//   for (unsigned i = 0; i < numNeighborCell(clistPro, bid); ++i){
//     __syncthreads();
//     IndexType targetCellIndex = getTargetCellIndex (clistPro, bid, i);
//     CoordType shift = getShiftValue (clistPro, bid, i);
//     IndexType targetIndex[tid] = get (clist, targetCellIndex, tid);
//     if (targetIndex[tid] != MaxIndexValue){
//       targetCoord[tid] = tex1Dfetch (global_texRef_interaction_coord, targetIndexes[tid]);
//       targetType[tid] = tex1Dfetch (global_texRef_interaction_type, targetIndexes[tid]);
//     }
//     __syncthreads ();
//     if (ii != MaxIndexValue){
//       for (IndexType jj = 0; jj < blockDim.x; ++jj){
// 	if (targetIndex[jj] == MaxIndexValue) continue;
// 	ScalorType diffx = targetCoord[jj].x + shift.x - ref.x;
// 	ScalorType diffy = targetCoord[jj].y + shift.y - ref.y;
// 	ScalorType diffz = targetCoord[jj].z + shift.z - ref.z;
// 	if ((diffx*diffx+diffy*diffy+diffz*diffz) < rlist2 &&
// 	    targetIndex[jj] != ii){
// 	  ForceIndexType fidx;
// 	  if (sharednbForceTable){
// 	    fidx = nonBondedInteractionTableItem (
// 		nonBondedInteractionTable, numAtomType, refType, targetType[jj]);
// 	  }
// 	  else {
// 	    fidx = nonBondedInteractionTableItem (
// 		nonBondedInteractionTable, numAtomType, refType, targetType[jj]);
// 	  }
// 	  ScalorType fx, fy, fz;
// 	  nbforce (nonBondedInteractionType[fidx],
// 		   &nonBondedInteractionParameter
// 		   [nonBondedInteractionParameterPosition[fidx]],
// 		   diffx, diffy, diffz,
// 		   &fx, &fy, &fz);
// 	  fsumx += fx;
// 	  fsumy += fy;
// 	  fsumz += fz;
// 	}
//       }
//     }
//   }

//   if (ii != MaxIndexValue){
//     forcx[ii] += fsumx;
//     forcy[ii] += fsumy;
//     forcz[ii] += fsumz;
//   }
// }

    

__global__ void calNonBondedInteraction (const IndexType numAtom,
					 const CoordType * coord,
					 ScalorType * forcx,
					 ScalorType * forcy, 
					 ScalorType * forcz,
					 const TypeType * type,
					 const RectangularBox box,
					 const DeviceNeighborList nlist,
					 mdError_t * ptr_de,
					 IndexType * errorIndex,
					 ScalorType * errorScalor)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  ScalorType fsumx = 0.0f;
  ScalorType fsumy = 0.0f;
  ScalorType fsumz = 0.0f;
  IndexType ii = tid + bid * blockDim.x;

  // IndexType num = 0;
  // IndexType start = bid * blockDim.x;
  // if (start < numAtom){
  //   if (start + blockDim.x <= numAtom){
  //     num = blockDim.x;
  //   }
  //   else {
  //     num = start + blockDim.x - numAtom;
  //   }
  // }
  // IndexType maxNei = maxVectorBlock (nlist.Nneighbor, start, num);
  // IndexType myNei = nlist.Nneighbor[ii];
  
  if (ii < numAtom) {
#ifdef COMPILE_NO_TEX
    CoordType ref (coord[ii]);
#else
    CoordType ref (tex1Dfetch(global_texRef_interaction_coord, ii));
#endif    
    ScalorType fx(0.f), fy(0.f), fz(0.f);
    for (IndexType jj = 0, nlistPosi = ii;
    	 jj < nlist.Nneighbor[ii];
    	 ++jj, nlistPosi += nlist.stride){

    // for (IndexType jj = 0, nlistPosi = ii;
    // 	 jj < maxNei;
    // 	 ++jj, nlistPosi += nlist.stride){
    //   __syncthreads();
    //   if (jj >= myNei) continue;
      
      IndexType targetIdx ( nlist.data [nlistPosi] );
      ForceIndexType nbForceIndex ( nlist.forceIndex [nlistPosi] );
#ifdef COMPILE_NO_TEX
      CoordType target ( coord[targetIdx] );
#else
      CoordType target ( tex1Dfetch(global_texRef_interaction_coord, targetIdx) );
#endif
      ScalorType diffx ( target.x - ref.x );
      ScalorType diffy ( target.y - ref.y );
      ScalorType diffz ( target.z - ref.z );
      shortestImage (box, &diffx, &diffy, &diffz);
      // ScalorType * forceParam;
      // NBForceSetting::getParam (nbForceIndex, nbForceParam, nbForceParamPosi,
      // 				&forceParam);
      // nbForce (nbForceType[nbForceIndex], forceParam,
      // 	       diffx, diffy, diffz, 
      // 	       &fx, &fy, &fz);
      nbForce (nonBondedInteractionType[nbForceIndex],
	       &nonBondedInteractionParameter
	       [nonBondedInteractionParameterPosition[nbForceIndex]],
      	       diffx, diffy, diffz, 
      	       &fx, &fy, &fz);
      fsumx += fx;
      fsumy += fy;
      fsumz += fz;
    }
    forcx[ii] += fsumx;
    forcy[ii] += fsumy;
    forcz[ii] += fsumz;
  }
}


__global__ void calNonBondedInteraction (const IndexType numAtom,
					 const CoordType * coord,
					 ScalorType * forcx,
					 ScalorType * forcy, 
					 ScalorType * forcz,
					 const TypeType * type,
					 const RectangularBox box,
					 const DeviceNeighborList nlist,
					 ScalorType * statistic_nb_buff0,
					 ScalorType * statistic_nb_buff1,
					 ScalorType * statistic_nb_buff2,
					 ScalorType * statistic_nb_buff3,
					 mdError_t * ptr_de,
					 IndexType * errorIndex,
					 ScalorType * errorScalor)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;

  ScalorType fsumx = 0.0f;
  ScalorType fsumy = 0.0f;
  ScalorType fsumz = 0.0f;
  IndexType ii = tid + bid * blockDim.x;
  ScalorType myPoten = 0.0f, myVxx = 0.0f, myVyy = 0.0f, myVzz = 0.0f;
  if (ii < numAtom) {
    CoordType ref;
#ifdef COMPILE_NO_TEX    
    ref = coord[ii];
#else
    ref = tex1Dfetch(global_texRef_interaction_coord, ii);
#endif
    ScalorType fx(0.f), fy(0.f), fz(0.f);
    ScalorType dp;
    for (IndexType jj = 0, nlistPosi = ii;
	 jj < nlist.Nneighbor[ii];
	 ++jj, nlistPosi += nlist.stride){
      IndexType targetIdx ( nlist.data[nlistPosi] );
      ForceIndexType nbForceIndex ( nlist.forceIndex [nlistPosi] );
#ifdef COMPILE_NO_TEX    
      CoordType target ( coord[targetIdx] );
#else
      CoordType target ( tex1Dfetch(global_texRef_interaction_coord, targetIdx) );
#endif
      ScalorType diffx ( target.x - ref.x );
      ScalorType diffy ( target.y - ref.y );
      ScalorType diffz ( target.z - ref.z );
      shortestImage (box, &diffx, &diffy, &diffz);
      // ScalorType * forceParam;
      // NBForceSetting::getParam (nbForceIndex, nbForceParam, nbForceParamPosi,
      // 				&forceParam);
      // nbForcePoten (nbForceType[nbForceIndex],
      // 		    forceParam,
      // 		    diffx, diffy, diffz, 
      // 		    &fx, &fy, &fz, &dp);
      nbForcePoten (nonBondedInteractionType[nbForceIndex],
		    &nonBondedInteractionParameter
		    [nonBondedInteractionParameterPosition[nbForceIndex]],
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
    forcx[ii] += fsumx;
    forcy[ii] += fsumy;
    forcz[ii] += fsumz;
  }
  
  if (ii < numAtom){
    statistic_nb_buff0[ii] = myPoten * 0.5f;
    statistic_nb_buff1[ii] = myVxx * 0.5f;
    statistic_nb_buff2[ii] = myVyy * 0.5f;
    statistic_nb_buff3[ii] = myVzz * 0.5f;
  }  
}


__global__ void calBondInteraction (const IndexType numAtom,
				    const CoordType * coord,
				    ScalorType * forcx,
				    ScalorType * forcy, 
				    ScalorType * forcz,
				    const RectangularBox box,
				    const DeviceBondList bdlist)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  ScalorType fsumx = 0.0f;
  ScalorType fsumy = 0.0f;
  ScalorType fsumz = 0.0f;
  IndexType ii = tid + bid * blockDim.x;
  
  if (ii >= numAtom) return;
  CoordType ref;
#ifdef COMPILE_NO_TEX
  ref = coord[ii];
#else
  ref = tex1Dfetch(global_texRef_interaction_coord, ii);
#endif  
      
  IndexType myNumBond = bdlist.numBond[ii];
  
  for (IndexType jj = 0; jj < bdlist.maxNumBond; ++jj){
    if (jj == myNumBond) break;
    IndexType targetIdx = bdlist.bondNeighborIndex[jj * bdlist.stride + ii];
    CoordType target;
#ifdef COMPILE_NO_TEX
    target = coord[targetIdx];
#else
    target = tex1Dfetch(global_texRef_interaction_coord, targetIdx);
#endif 
    ScalorType diffx, diffy, diffz;
    diffx = target.x - ref.x;
    diffy = target.y - ref.y;
    diffz = target.z - ref.z;
    shortestImage (box, &diffx, &diffy, &diffz);
    ScalorType fx, fy, fz;
    ForceIndexType bondFindex = bdlist.bondIndex[jj * bdlist.stride + ii];
    bondForce (bondedInteractionType[bondFindex],
	       &bondedInteractionParameter
	       [bondedInteractionParameterPosition[bondFindex]],
	       diffx, diffy, diffz, &fx, &fy, &fz);
    fsumx += fx;
    fsumy += fy;
    fsumz += fz;
  }
  forcx[ii] += fsumx;
  forcy[ii] += fsumy;
  forcz[ii] += fsumz;
}


__global__ void calBondInteraction (const IndexType numAtom,
				    const CoordType * coord,
				    ScalorType * forcx,
				    ScalorType * forcy, 
				    ScalorType * forcz,
				    const RectangularBox box,
				    const DeviceBondList bdlist,
				    ScalorType * statistic_b_buff0,
				    ScalorType * statistic_b_buff1,
				    ScalorType * statistic_b_buff2,
				    ScalorType * statistic_b_buff3,
				    mdError_t * ptr_de)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;

  extern __shared__ volatile ScalorType buff[];
  buff[tid] = 0.f;
  __syncthreads();
  
  ScalorType fsumx = 0.0f;
  ScalorType fsumy = 0.0f;
  ScalorType fsumz = 0.0f;
  IndexType ii = tid + bid * blockDim.x;
  ScalorType myPoten = 0.0f, myVxx = 0.0f, myVyy = 0.0f, myVzz = 0.0f;
  if (ii < numAtom) {
    CoordType ref;
#ifdef COMPILE_NO_TEX
    ref = coord[ii];
#else 
    ref = tex1Dfetch(global_texRef_interaction_coord, ii);
#endif
    IndexType myNumBond = bdlist.numBond[ii];
    for (IndexType jj = 0; jj < bdlist.maxNumBond; ++jj){
      if (jj == myNumBond) break;
      IndexType targetIdx = bdlist.bondNeighborIndex[jj * bdlist.stride + ii];
      CoordType target;
#ifdef COMPILE_NO_TEX
      target = coord[targetIdx];
#else
      target = tex1Dfetch(global_texRef_interaction_coord, targetIdx);
#endif
      ScalorType diffx, diffy, diffz;
      diffx = target.x - ref.x;
      diffy = target.y - ref.y;
      diffz = target.z - ref.z;
      shortestImage (box, &diffx, &diffy, &diffz);
      ScalorType fx, fy, fz;
      ForceIndexType bondFindex = bdlist.bondIndex[jj * bdlist.stride + ii];
      ScalorType dp;
      bondForcePoten (bondedInteractionType[bondFindex],
		      &bondedInteractionParameter
		      [bondedInteractionParameterPosition[bondFindex]],
		      diffx, diffy, diffz, &fx, &fy, &fz, &dp);
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



__global__ void calAngleInteraction (const IndexType numAtom,
				     const CoordType * coord,
				     ScalorType * forcx,
				     ScalorType * forcy, 
				     ScalorType * forcz,
				     const RectangularBox box,
				     const DeviceAngleList anglelist)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  ScalorType fsumx = 0.0f;
  ScalorType fsumy = 0.0f;
  ScalorType fsumz = 0.0f;
  IndexType ii = tid + bid * blockDim.x;
  IndexType myNumAngle;
  
  if (ii < numAtom){
    myNumAngle = anglelist.numAngle[ii];  
  }
  else {
    myNumAngle = 0;
  }
  if (__all(myNumAngle == 0)) return ;
  
  if (ii < numAtom) {
    CoordType ref;
#ifdef COMPILE_NO_TEX
    ref = coord[ii];
#else
    ref = tex1Dfetch(global_texRef_interaction_coord, ii);
#endif
    for (IndexType jj = 0; jj < myNumAngle; ++jj){
      IndexType targetIdx0 =
	  anglelist.angleNeighborIndex[((jj<<1)  ) * anglelist.stride + ii];
      IndexType targetIdx1 =
	  anglelist.angleNeighborIndex[((jj<<1)+1) * anglelist.stride + ii];
      IndexType myPosi     =
	  anglelist.anglePosi[jj * anglelist.stride + ii];
      CoordType target0, target1;
#ifdef COMPILE_NO_TEX
      target0 = coord[targetIdx0];
      target1 = coord[targetIdx1];
#else
      target0 = tex1Dfetch(global_texRef_interaction_coord, targetIdx0);
      target1 = tex1Dfetch(global_texRef_interaction_coord, targetIdx1);
#endif 
      ScalorType diff0x, diff0y, diff0z;
      ScalorType diff1x, diff1y, diff1z;
      bool center (myPosi == 1);
      if (center){
	diff0x = ref.x - target0.x;
	diff0y = ref.y - target0.y;
	diff0z = ref.z - target0.z;
	diff1x = target1.x -  ref.x;
	diff1y = target1.y -  ref.y;
	diff1z = target1.z -  ref.z;
      } else {
	diff0x = target0.x - ref.x;
	diff0y = target0.y - ref.y;
	diff0z = target0.z - ref.z;
	diff1x = target1.x - target0.x;
	diff1y = target1.y - target0.y;
	diff1z = target1.z - target0.z;
      }      
      shortestImage (box, &diff0x, &diff0y, &diff0z);
      shortestImage (box, &diff1x, &diff1y, &diff1z);
      ScalorType fx, fy, fz;
      ForceIndexType angleFindex = anglelist.angleIndex[jj * anglelist.stride + ii];
      angleForce (center,
		  bondedInteractionType[angleFindex],
		  &bondedInteractionParameter
		  [bondedInteractionParameterPosition[angleFindex]],
		  diff0x, diff0y, diff0z,
		  diff1x, diff1y, diff1z,
		  &fx, &fy, &fz);
      fsumx += fx;
      fsumy += fy;
      fsumz += fz;
    }
    forcx[ii] += fsumx;
    forcy[ii] += fsumy;
    forcz[ii] += fsumz;
  }
}


__global__ void calAngleInteraction (const IndexType numAtom,
				     const CoordType * coord,
				     ScalorType * forcx,
				     ScalorType * forcy, 
				     ScalorType * forcz,
				     const RectangularBox box,
				     const DeviceAngleList anglelist,
				     ScalorType * statistic_b_buff0,
				     mdError_t * ptr_de)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  ScalorType fsumx = 0.0f;
  ScalorType fsumy = 0.0f;
  ScalorType fsumz = 0.0f;
  ScalorType myPoten = 0.f;
  IndexType ii = tid + bid * blockDim.x;
  IndexType myNumAngle;
  extern __shared__ volatile ScalorType buff[];
  buff[tid] = 0.f;
  __syncthreads();
  
  if (ii < numAtom) {
    CoordType ref;
#ifdef COMPILE_NO_TEX
    ref = coord[ii];
#else
    ref = tex1Dfetch(global_texRef_interaction_coord, ii);
#endif
    myNumAngle = anglelist.numAngle[ii];  
    for (IndexType jj = 0; jj < myNumAngle; ++jj){
      IndexType targetIdx0 =
	  anglelist.angleNeighborIndex[((jj<<1)  ) * anglelist.stride + ii];
      IndexType targetIdx1 =
	  anglelist.angleNeighborIndex[((jj<<1)+1) * anglelist.stride + ii];
      IndexType myPosi     =
	  anglelist.anglePosi[jj * anglelist.stride + ii];
      CoordType target0, target1;
#ifdef COMPILE_NO_TEX
      target0 = coord[targetIdx0];
      target1 = coord[targetIdx1];
#else
      target0 = tex1Dfetch(global_texRef_interaction_coord, targetIdx0);
      target1 = tex1Dfetch(global_texRef_interaction_coord, targetIdx1);
#endif 
      ScalorType diff0x, diff0y, diff0z;
      ScalorType diff1x, diff1y, diff1z;
      bool center = (myPosi == 1);
      if (center){
	diff0x = ref.x - target0.x;
	diff0y = ref.y - target0.y;
	diff0z = ref.z - target0.z;
	diff1x = target1.x -  ref.x;
	diff1y = target1.y -  ref.y;
	diff1z = target1.z -  ref.z;
      } else {
	diff0x = target0.x - ref.x;
	diff0y = target0.y - ref.y;
	diff0z = target0.z - ref.z;
	diff1x = target1.x - target0.x;
	diff1y = target1.y - target0.y;
	diff1z = target1.z - target0.z;
      }      
      shortestImage (box, &diff0x, &diff0y, &diff0z);
      shortestImage (box, &diff1x, &diff1y, &diff1z);
      ScalorType fx, fy, fz;
      ForceIndexType angleFindex = anglelist.angleIndex[jj * anglelist.stride + ii];
      ScalorType dp;
      angleForcePoten (center,
		       bondedInteractionType[angleFindex],
		       &bondedInteractionParameter
		       [bondedInteractionParameterPosition[angleFindex]],
		       diff0x, diff0y, diff0z,
		       diff1x, diff1y, diff1z,
		       &fx, &fy, &fz, &dp);
      myPoten += dp;
      fsumx += fx;
      fsumy += fy;
      fsumz += fz;
    }
    forcx[ii] += fsumx;
    forcy[ii] += fsumy;
    forcz[ii] += fsumz;
  }

  buff[tid] = myPoten * 0.33333333333333333f;
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) statistic_b_buff0[bid] = buff[0];
}



static __device__ IndexType shiftedD3toD1 (
    DeviceCellList clist,
    RectangularBoxGeometry::RectangularBox box,
    int ix, int iy, int iz,
    ScalorType * shiftx , ScalorType * shifty, ScalorType * shiftz)
{
  int tmp;
  ix += (tmp = -int(floorf(ix * clist.NCelli.x))) * clist.NCell.x;
  *shiftx = tmp * box.size.x;
  iy += (tmp = -int(floorf(iy * clist.NCelli.y))) * clist.NCell.y;
  *shifty = tmp * box.size.y;
  iz += (tmp = -int(floorf(iz * clist.NCelli.z))) * clist.NCell.z;
  *shiftz = tmp * box.size.z;
  return D3toD1 (clist, ix, iy, iz);
}

__global__ void calNonBondedInteraction (
    const IndexType numAtom,
    const CoordType * coord,
    ScalorType * forcx,
    ScalorType * forcy, 
    ScalorType * forcz,
    const TypeType * type,
    const RectangularBox box,
    DeviceCellList clist,
    mdError_t * ptr_de)
{
  // RectangularBoxGeometry::normalizeSystem (box, &ddata);
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType bidx, bidy, bidz;
  D1toD3 (clist, bid, bidx, bidy, bidz);
  
  // load index
  IndexType ii = getDeviceCellListData (clist, bid, tid);
  // load iith coordinate // use texturefetch instead
  CoordType ref;
  TypeType reftype;
  ScalorType fsumx (0.f), fsumy(0.f), fsumz(0.f);
  if (ii != MaxIndexValue){
#ifdef COMPILE_NO_TEX
    ref = coord[ii];
    reftype = type[ii];
#else
    ref = tex1Dfetch (global_texRef_interaction_coord, ii);
    // reftype = tex1Dfetch(global_texRef_interaction_type, ii);
#endif
  }
  ScalorType rlist = clist.rlist;

  // the target index and coordinates are shared

  extern __shared__ volatile char pub_sbuff[];
  
  volatile IndexType * targetIndexes =
      (volatile IndexType *) pub_sbuff;
  CoordType * target =
      (CoordType *) &targetIndexes[roundUp4(blockDim.x)];
  volatile TypeType * targettype =
      (volatile TypeType *) &target[roundUp4(blockDim.x)];
  __syncthreads();

  // bool oneCellX(false), oneCellY(false), oneCellZ(false);
  // if (clist.NCell.x == 1) oneCellX = true;
  // if (clist.NCell.y == 1) oneCellY = true;
  // if (clist.NCell.z == 1) oneCellZ = true;
  // int upperx(1), lowerx(-1);
  // int uppery(1), lowery(-1);
  // int upperz(1), lowerz(-1);
  // if (oneCellX) {lowerx =  0; upperx = 0;}
  // if (oneCellY) {lowery =  0; uppery = 0;}
  // if (oneCellZ) {lowerz =  0; upperz = 0;}
  ScalorType rlist2 = rlist * rlist;
  
  // loop over 27 neighbor cells
#pragma unroll 3
  // for (int nci = bidx + lowerx; nci <= bidx + upperx; ++nci){
  //   for (int ncj = bidy + lowery; ncj <= bidy + uppery; ++ncj){
  //     for (int nck = bidz + lowerz; nck <= bidz + upperz; ++nck){
  for (int nci = int(bidx) - 1; nci <= int(bidx) + 1; ++nci){
    for (int ncj = int(bidy) - 1; ncj <= int(bidy) + 1; ++ncj){
      for (int nck = int(bidz) - 1; nck <= int(bidz) + 1; ++nck){
  // for (int di = lowerx; di <= upperx; ++di){
  //   for (int dj = lowery; dj <= uppery; ++dj){
  //     for (int dk = lowerz; dk <= upperz; ++dk){
	__syncthreads();
	// the shift value of a cell is pre-computed
	ScalorType xshift, yshift, zshift;
	// int nci = di + bidx;
	// int ncj = dj + bidy;
	// int nck = dk + bidz;
	IndexType targetCellIdx = shiftedD3toD1 (clist, box, 
						 nci, ncj, nck, 
						 &xshift, &yshift, &zshift);
	// load target index and coordinates form global memary
	// IndexType tmp = (targetIndexes[tid] = 
	// 		 getDeviceCellListData(clist, targetCellIdx, tid));
	targetIndexes[tid] = getDeviceCellListData(clist, targetCellIdx, tid);
	if (targetIndexes[tid] != MaxIndexValue){
// #ifdef COMPILE_NO_TEX
// 	  target[tid] = coord[tmp];
// 	  // targettype[tid] = type[tmp];
// #else
	  target[tid] = tex1Dfetch(global_texRef_interaction_coord, targetIndexes[tid]);
	  // targettype[tid] = tex1Dfetch(global_texRef_interaction_type, tmp);
// #endif
	}
	__syncthreads();
	// find neighbor
	// if (ii != MaxIndexValue){
	  for (IndexType jj = 0; jj < clist.numbers[targetCellIdx]; ++jj){
	    // if (targetIndexes[jj] == MaxIndexValue) break;
	    ScalorType diffx = target[jj].x - xshift - ref.x;
	    ScalorType diffy = target[jj].y - yshift - ref.y;
	    ScalorType diffz = target[jj].z - zshift - ref.z;
	    // if (oneCellX) shortestImage (box.size.x, box.sizei.x, &diffx);
	    // if (oneCellY) shortestImage (box.size.y, box.sizei.y, &diffy);
	    // if (oneCellZ) shortestImage (box.size.z, box.sizei.z, &diffz);
	    //printf ("%d\t%d\t%f\t%f\n", ii,
	    ScalorType dr2;
	    if ((dr2 = (diffx*diffx+diffy*diffy+diffz*diffz)) < rlist2 &&
		targetIndexes[jj] != ii){
	      ForceIndexType fidx(0);
	      // fidx = AtomNBForceTable::calForceIndex (
	      // 	  nonBondedInteractionTable,
	      // 	  numAtomType[0],
	      // 	  reftype,
	      // 	  targettype[jj]);
	      // if (fidx != mdForceNULL) {
	      ScalorType fx, fy, fz;
	      nbForce (nonBondedInteractionType[fidx],
		       &nonBondedInteractionParameter
		       [nonBondedInteractionParameterPosition[fidx]],
		       diffx, diffy, diffz,
		       dr2,
		       &fx, &fy, &fz);
	      fsumx += fx;
	      fsumy += fy;
	      fsumz += fz;
	      // }
	    }
	  }
	// }
      }
    }
  }
  if (ii != MaxIndexValue){
    forcx[ii] += fsumx;
    forcy[ii] += fsumy;
    forcz[ii] += fsumz;
  }
}






__global__ void calNonBondedInteraction (
    const IndexType numAtom,
    const CoordType * coord,
    ScalorType * forcx,
    ScalorType * forcy, 
    ScalorType * forcz,
    const TypeType * type,
    const RectangularBox box,
    DeviceCellList clist,
    ScalorType * statistic_nb_buff0,
    ScalorType * statistic_nb_buff1,
    ScalorType * statistic_nb_buff2,
    ScalorType * statistic_nb_buff3,
    mdError_t * ptr_de)
{
  // RectangularBoxGeometry::normalizeSystem (box, &ddata);
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType bidx, bidy, bidz;
  D1toD3 (clist, bid, bidx, bidy, bidz);
  
  // load index
  IndexType ii = getDeviceCellListData (clist, bid, tid);
  // load iith coordinate // use texturefetch instead
  CoordType ref;
  TypeType reftype;
  ScalorType fsumx (0.f), fsumy(0.f), fsumz(0.f);
  ScalorType myPoten (0.0f), myVxx (0.0f), myVyy (0.0f), myVzz (0.0f);
  if (ii != MaxIndexValue){
#ifdef COMPILE_NO_TEX
    ref = coord[ii];
    reftype = type[ii];
#else
    ref = tex1Dfetch (global_texRef_interaction_coord, ii);
    reftype = tex1Dfetch(global_texRef_interaction_type, ii);
#endif
  }
  ScalorType rlist = clist.rlist;

  // the target index and coordinates are shared

  extern __shared__ volatile char pub_sbuff[];
  
  volatile IndexType * targetIndexes =
      (volatile IndexType *) pub_sbuff;
  CoordType * target =
      (CoordType *) &targetIndexes[roundUp4(blockDim.x)];
  volatile TypeType * targettype =
      (volatile TypeType *) &target[roundUp4(blockDim.x)];
  __syncthreads();

  bool oneCellX(false), oneCellY(false), oneCellZ(false);
  if (clist.NCell.x == 1) oneCellX = true;
  if (clist.NCell.y == 1) oneCellY = true;
  if (clist.NCell.z == 1) oneCellZ = true;
  int upperx(1), lowerx(-1);
  int uppery(1), lowery(-1);
  int upperz(1), lowerz(-1);
  if (oneCellX) {lowerx =  0; upperx = 0;}
  if (oneCellY) {lowery =  0; uppery = 0;}
  if (oneCellZ) {lowerz =  0; upperz = 0;}
  ScalorType rlist2 = rlist * rlist;
  
  // loop over 27 neighbor cells
#pragma unroll 3
  for (int di = lowerx; di <= upperx; ++di){
    for (int dj = lowery; dj <= uppery; ++dj){
      for (int dk = lowerz; dk <= upperz; ++dk){
	__syncthreads();
	// the shift value of a cell is pre-computed
	ScalorType xshift, yshift, zshift;
	int nci = di + bidx;
	int ncj = dj + bidy;
	int nck = dk + bidz;
	IndexType targetCellIdx = shiftedD3toD1 (clist, box, 
						 nci, ncj, nck, 
						 &xshift, &yshift, &zshift);
	// load target index and coordinates form global memary
	IndexType tmp = (targetIndexes[tid] = 
			 getDeviceCellListData(clist, targetCellIdx, tid));
	if (tmp != MaxIndexValue){
#ifdef COMPILE_NO_TEX
	  target[tid] = coord[tmp];
	  targettype[tid] = type[tmp];
#else
	  target[tid] = tex1Dfetch(global_texRef_interaction_coord, tmp);
	  targettype[tid] = tex1Dfetch(global_texRef_interaction_type, tmp);
#endif
	}
	__syncthreads();
	// find neighbor
	if (ii != MaxIndexValue){
	  for (IndexType jj = 0; jj < blockDim.x; ++jj){
	    if (targetIndexes[jj] == MaxIndexValue) break;
	    ScalorType diffx = target[jj].x - xshift - ref.x;
	    ScalorType diffy = target[jj].y - yshift - ref.y;
	    ScalorType diffz = target[jj].z - zshift - ref.z;
	    if (oneCellX) shortestImage (box.size.x, box.sizei.x, &diffx);
	    if (oneCellY) shortestImage (box.size.y, box.sizei.y, &diffy);
	    if (oneCellZ) shortestImage (box.size.z, box.sizei.z, &diffz);
	    //printf ("%d\t%d\t%f\t%f\n", ii, 
	    if ((diffx*diffx+diffy*diffy+diffz*diffz) < rlist2 &&
		targetIndexes[jj] != ii){
	      ForceIndexType fidx;
	      fidx = AtomNBForceTable::calForceIndex (
		  nonBondedInteractionTable,
		  numAtomType[0],
		  reftype,
		  targettype[jj]);
	      // if (fidx != mdForceNULL) {
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
	      // }
	    }
	  }
	}
      }
    }
  }
  if (ii != MaxIndexValue){
    forcx[ii] += fsumx;
    forcy[ii] += fsumy;
    forcz[ii] += fsumz;
    statistic_nb_buff0[ii] = myPoten * 0.5f;
    statistic_nb_buff1[ii] = myVxx * 0.5f;
    statistic_nb_buff2[ii] = myVyy * 0.5f;
    statistic_nb_buff3[ii] = myVzz * 0.5f;
  }  
}

