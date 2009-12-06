#include "InteractionEngine_interface.h"

texture<ScalorType, 1, cudaReadModeElementType> global_texRef_interaction_coordx;
texture<ScalorType, 1, cudaReadModeElementType> global_texRef_interaction_coordy;
texture<ScalorType, 1, cudaReadModeElementType> global_texRef_interaction_coordz;
// texture<TypeType, 1, cudaReadModeElementType> global_texRef_interaction_type;

__constant__ mdNBInteraction_t nbForceType [MaxNumberNBForce];
__constant__ ScalorType nbForceParam       [MaxNumberNBForceParam];
__constant__ IndexType  nbForceParamPosi   [MaxNumberNBForce];
__constant__ mdBondInteraction_t bondForceType [MaxNumberBondForce];
__constant__ ScalorType bondForceParam         [MaxNumberBondForceParam];
__constant__ IndexType  bondForceParamPosi     [MaxNumberBondForce];

void InteractionEngine_interface::init (const MDSystem  & sys,
					const IndexType & NTread)
{
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

  size_t sizescalor = sizeof(ScalorType)*sys.ddata.numMem;
  // size_t sizetype = sizeof(TypeType)*sys.ddata.numMem;
  cudaBindTexture(0, global_texRef_interaction_coordx, sys.ddata.coordx, sizescalor);
  cudaBindTexture(0, global_texRef_interaction_coordy, sys.ddata.coordy, sizescalor);
  cudaBindTexture(0, global_texRef_interaction_coordz, sys.ddata.coordz, sizescalor);
  // cudaBindTexture(0, global_texRef_interaction_type, sys.ddata.type, sizetype);
  checkCUDAError ("InteractionEngine::init, bind texture");

  // init statistical buff
  // cudaMalloc((void**)&(statistic_nb_buff0), sizeof(ScalorType) * nob);
  // cudaMalloc((void**)&(statistic_nb_buff1), sizeof(ScalorType) * nob);
  // cudaMalloc((void**)&(statistic_nb_buff2), sizeof(ScalorType) * nob);
  // cudaMalloc((void**)&(statistic_nb_buff3), sizeof(ScalorType) * nob);
  // cudaMalloc((void**)&(statistic_b_buff0), sizeof(ScalorType) * nob);
  // cudaMalloc((void**)&(statistic_b_buff1), sizeof(ScalorType) * nob);
  // cudaMalloc((void**)&(statistic_b_buff2), sizeof(ScalorType) * nob);
  // cudaMalloc((void**)&(statistic_b_buff3), sizeof(ScalorType) * nob);
  // checkCUDAError ("InteractionEngine::init, allocate for statistic buff");
  
  // init sum vectors
  sum_nb_p.init (sys.ddata.numAtom, NThreadForSum);
  sum_nb_vxx.init (sys.ddata.numAtom, NThreadForSum);
  sum_nb_vyy.init (sys.ddata.numAtom, NThreadForSum);
  sum_nb_vzz.init (sys.ddata.numAtom, NThreadForSum);
  sum_b_p.init (nob, NThreadForSum);
  sum_b_vxx.init (nob, NThreadForSum);
  sum_b_vyy.init (nob, NThreadForSum);
  sum_b_vzz.init (nob, NThreadForSum);
  for (IndexType i = 0; i < 8; ++i){
    cudaStreamCreate(&sum_stream[i]);
  }
  checkCUDAError ("InteractionEngine::init init sum statistic");
  
  // init nb force param
  if (sys.nbForce.setting.NNBForce > MaxNumberNBForce ){
    throw MDExcptExceedConstantMemLimit ("InteractionEngine::init", "nbForceType",
					 MaxNumberNBForce * sizeof(mdNBInteraction_t));
  }
  if (sys.nbForce.setting.paramLength > MaxNumberNBForceParam ){
    throw MDExcptExceedConstantMemLimit ("InteractionEngine::init", "nbForceParam",
					 MaxNumberNBForceParam * sizeof(ScalorType));
  }

  // cudaMemset (nbForceType, 0,
  // 	      sizeof(mdNBInteraction_t) * sys.nbForce.setting.NNBForce);
  cudaMemcpyToSymbol (nbForceType, sys.nbForce.setting.type, 
  		      sizeof(mdNBInteraction_t) * sys.nbForce.setting.NNBForce);
  // printf ("# cpy nb force type, size %d N %d\n",
  // 	  sizeof(mdNBInteraction_t) * sys.nbForce.setting.NNBForce,
  // 	  sys.nbForce.setting.NNBForce);
  for (IndexType i = 0; i < sys.nbForce.setting.NNBForce; ++i){
    printf ("# %d\n", sys.nbForce.setting.type[i]);
  }
  // cudaMemset (nbForceParam, 0,
  // 	      sizeof(ScalorType) * sys.nbForce.setting.paramLength);
  // cudaThreadSynchronize ();
  cudaMemcpyToSymbol (nbForceParam, sys.nbForce.setting.param,
  		      sizeof(ScalorType) * sys.nbForce.setting.paramLength);
  // cudaThreadSynchronize ();
  // printf ("# cpy nb force param, N%d\n", sys.nbForce.setting.paramLength);
  // for (IndexType i = 0; i < sys.nbForce.setting.paramLength; ++i){
  //   printf("# %02d %f\n", i, sys.nbForce.setting.param[i]);
  // }
  // cudaMemset (nbForceParamPosi, 0,
  // 	      sizeof(IndexType) * sys.nbForce.setting.NNBForce);
  cudaMemcpyToSymbol (nbForceParamPosi, sys.nbForce.setting.paramPosi,
		      sizeof(IndexType) * sys.nbForce.setting.NNBForce);
  // printf ("# cpy nb force param posi\n");
  checkCUDAError ("InteractionEngine::init, init NB force setting");
  // for (IndexType i = 0; i < sys.nbForce.setting.NNBForce; ++i){
  //   printf ("# %d\n", sys.nbForce.setting.paramPosi[i]);
  // }
  
  //init bond force param
  cudaMemcpyToSymbol (bondForceType, sys.bdlist.bondType,
		      sizeof(mdBondInteraction_t) * sys.bdlist.NBondForce);
  cudaMemcpyToSymbol (bondForceParam, sys.bdlist.param,
		      sizeof(ScalorType) * sys.bdlist.paramLength);
  cudaMemcpyToSymbol (bondForceParamPosi, sys.bdlist.paramPosi,
		      sizeof(IndexType) * sys.bdlist.NBondForce);
  checkCUDAError ("InteractionEngine::init, init bond force setting");

  // cal shared buff size
  calBondInteraction_sbuffSize = myBlockDim.x * sizeof(ScalorType);
}

InteractionEngine_interface::~InteractionEngine_interface()
{
  cudaUnbindTexture(global_texRef_interaction_coordx);
  cudaUnbindTexture(global_texRef_interaction_coordy);
  cudaUnbindTexture(global_texRef_interaction_coordz);
  // cudaUnbindTexture(global_texRef_interaction_type);

  // cudaFree(statistic_nb_buff0);
  // cudaFree(statistic_nb_buff1);
  // cudaFree(statistic_nb_buff2);
  // cudaFree(statistic_nb_buff3);  

  // cudaFree(statistic_b_buff0);
  // cudaFree(statistic_b_buff1);
  // cudaFree(statistic_b_buff2);
  // cudaFree(statistic_b_buff3);  
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

	  

void InteractionEngine_interface::applyInteraction (MDSystem & sys,
						    const NeighborList & nlist,
						    MDTimer *timer )
{
  clearInteraction (sys);
  if (timer != NULL) timer->tic(mdTimeNonBondedInteraction);
  calNonBondedInteraction
      <<<atomGridDim, myBlockDim>>> (
	  sys.ddata.numAtom,
	  sys.ddata.coordx, sys.ddata.coordy, sys.ddata.coordz,
	  sys.ddata.forcx,  sys.ddata.forcy,  sys.ddata.forcz,
	  sys.ddata.type, 
	  sys.box, nlist.dnlist,
	  err.ptr_de,
	  err.ptr_dindex,
	  err.ptr_dscalor);
  checkCUDAError ("InteractionEngine::applyInteraction nb");
  err.check ("interaction engine nb");	
  if (timer != NULL) timer->toc(mdTimeNonBondedInteraction);
  if (timer != NULL) timer->tic(mdTimeBondedInteraction);
  calBondInteraction
      <<<atomGridDim, myBlockDim>>> (
  	  sys.ddata.numAtom,
  	  sys.ddata.coordx, sys.ddata.coordy, sys.ddata.coordz,
  	  sys.ddata.forcx,  sys.ddata.forcy,  sys.ddata.forcz,
  	  sys.box, sys.bdlist.dbdlist);
  checkCUDAError ("InteractionEngine::applyInteraction bonded");
  err.check ("interaction engine b");	
  if (timer != NULL) timer->toc(mdTimeBondedInteraction);
}

void InteractionEngine_interface::applyInteraction (MDSystem &sys,
						    const NeighborList & nlist,
						    MDStatistic & st,
						    MDTimer *timer )
{
  clearInteraction (sys);
  if (timer != NULL) timer->tic(mdTimeNBInterStatistic);
  calNonBondedInteraction
      <<<atomGridDim, myBlockDim>>> (
	  sys.ddata.numAtom,
	  sys.ddata.coordx, sys.ddata.coordy, sys.ddata.coordz,
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
  if (timer != NULL) timer->toc(mdTimeNBInterStatistic);
  if (timer != NULL) timer->tic(mdTimeBInterStatistic);
  calBondInteraction
      <<<atomGridDim, myBlockDim,
      calBondInteraction_sbuffSize>>> (
	  sys.ddata.numAtom,
	  sys.ddata.coordx, sys.ddata.coordy, sys.ddata.coordz,
	  sys.ddata.forcx,  sys.ddata.forcy,  sys.ddata.forcz,
	  sys.box, sys.bdlist.dbdlist
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

  if (timer != NULL) timer->tic(mdTimeNBInterStatistic);
  cudaThreadSynchronize();
  sum_nb_p.sumBuffAdd(st.ddata, mdStatisticNonBondedPotential, 0);
  sum_nb_vxx.sumBuffAdd(st.ddata, mdStatisticVirialXX, 1);
  sum_nb_vyy.sumBuffAdd(st.ddata, mdStatisticVirialYY, 2);
  sum_nb_vzz.sumBuffAdd(st.ddata, mdStatisticVirialZZ, 3);
  cudaThreadSynchronize();
  if (timer != NULL) timer->toc(mdTimeNBInterStatistic);
  if (timer != NULL) timer->tic(mdTimeBInterStatistic);
  cudaThreadSynchronize();
  sum_b_p.sumBuffAdd(st.ddata, mdStatisticBondedPotential, 4);
  sum_b_vxx.sumBuffAdd(st.ddata, mdStatisticVirialXX, 5);
  sum_b_vyy.sumBuffAdd(st.ddata, mdStatisticVirialYY, 6);
  sum_b_vzz.sumBuffAdd(st.ddata, mdStatisticVirialZZ, 7);
  cudaThreadSynchronize();
  if (timer != NULL) timer->toc(mdTimeBInterStatistic);
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

__global__ void calNonBondedInteraction (const IndexType numAtom,
					 const ScalorType * coordx,
					 const ScalorType * coordy, 
					 const ScalorType * coordz,
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
  
  if (ii < numAtom) {
    ScalorType refx = 0.0f, refy = 0.0f, refz = 0.0f;
#ifdef COMPILE_NO_TEX
    refx = coordx[ii];
    refy = coordy[ii];
    refz = coordz[ii];
#else
    refx = tex1Dfetch(global_texRef_interaction_coordx, ii);
    refy = tex1Dfetch(global_texRef_interaction_coordy, ii);
    refz = tex1Dfetch(global_texRef_interaction_coordz, ii);
#endif
    IndexType myNumNei = nlist.Nneighbor[ii];
    for (IndexType jj = 0; jj < myNumNei; ++jj){
      IndexType nlistPosi = jj * nlist.stride + ii;
      IndexType targetIdx   = nlist.data [nlistPosi];
      ForceIndexType nbForceIndex = nlist.forceIndex [nlistPosi];
      ScalorType targetx, targety, targetz;
#ifdef COMPILE_NO_TEX
      targetx = coordx[targetIdx];
      targety = coordy[targetIdx];
      targetz = coordz[targetIdx];
#else
      targetx = tex1Dfetch(global_texRef_interaction_coordx, targetIdx);
      targety = tex1Dfetch(global_texRef_interaction_coordy, targetIdx);
      targetz = tex1Dfetch(global_texRef_interaction_coordz, targetIdx);
#endif
      ScalorType diffx, diffy, diffz;
      diffx = targetx - refx;
      diffy = targety - refy;
      diffz = targetz - refz;
      shortestImage (box, &diffx, &diffy, &diffz);
      ScalorType fx(0.f), fy(0.f), fz(0.f);
      ScalorType * forceParam;
      NBForceSetting::getParam (nbForceIndex, nbForceParam, nbForceParamPosi,
				&forceParam);
      nbForce (nbForceType[nbForceIndex], forceParam,
      	       diffx, diffy, diffz, 
      	       &fx, &fy, &fz);
      fsumx -= fx;
      fsumy -= fy;
      fsumz -= fz;
    }
    forcx[ii] += fsumx;
    forcy[ii] += fsumy;
    forcz[ii] += fsumz;
  }
}


__global__ void calNonBondedInteraction (const IndexType numAtom,
					 const ScalorType * coordx,
					 const ScalorType * coordy, 
					 const ScalorType * coordz,
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
    ScalorType refx = 0.0f, refy = 0.0f, refz = 0.0f;
#ifdef COMPILE_NO_TEX    
    refx = coordx[ii];
    refy = coordy[ii];
    refz = coordz[ii];
#else
    refx = tex1Dfetch(global_texRef_interaction_coordx, ii);
    refy = tex1Dfetch(global_texRef_interaction_coordy, ii);
    refz = tex1Dfetch(global_texRef_interaction_coordz, ii);
#endif
    IndexType myNumNei = nlist.Nneighbor[ii];
    for (IndexType jj = 0; jj < myNumNei; ++jj){
      IndexType nlistPosi = jj * nlist.stride + ii;
      IndexType targetIdx = nlist.data[nlistPosi];
      ForceIndexType nbForceIndex = nlist.forceIndex [nlistPosi];
      ScalorType targetx, targety, targetz;
#ifdef COMPILE_NO_TEX    
      targetx = coordx[targetIdx];
      targety = coordy[targetIdx];
      targetz = coordz[targetIdx];
#else
      targetx = tex1Dfetch(global_texRef_interaction_coordx, targetIdx);
      targety = tex1Dfetch(global_texRef_interaction_coordy, targetIdx);
      targetz = tex1Dfetch(global_texRef_interaction_coordz, targetIdx);
#endif
      ScalorType diffx, diffy, diffz;
      diffx = targetx - refx;
      diffy = targety - refy;
      diffz = targetz - refz;
      shortestImage (box, &diffx, &diffy, &diffz);
      ScalorType fx(0.f), fy(0.f), fz(0.f);
      ScalorType * forceParam;
      NBForceSetting::getParam (nbForceIndex, nbForceParam, nbForceParamPosi,
      				&forceParam);
      ScalorType dp;
      // if (nbForceType[nbForceIndex] != 3){
      // 	*ptr_de = mdErrorShortNeighborList;
      // 	break;
      // }
      // dp = CosTail::forcePoten (forceParam,  diffx, diffy, diffz, &fx, &fy, &fz);

      nbForcePoten (nbForceType[nbForceIndex],
      		    forceParam,
      		    diffx, diffy, diffz, 
      		    &fx, &fy, &fz, &dp);
      myPoten += dp;
      myVxx -= fx * diffx;
      myVyy -= fy * diffy;
      myVzz -= fz * diffz;
      fsumx -= fx;
      fsumy -= fy;
      fsumz -= fz;
    }
    forcx[ii] += fsumx;
    forcy[ii] += fsumy;
    forcz[ii] += fsumz;
  }

  // __syncthreads();  
  
  // IndexType num = ((bid+1) * blockDim.x > numAtom) ? 
  //     (blockDim.x - ((bid+1) * blockDim.x) + numAtom) : blockDim.x;
  // ScalorType partialSum;
  // ScalorType tmpP(0.f), tmpVxx(0.f), tmpVyy(0.f), tmpVzz(0.f);
  
  // __syncthreads();
  // extern __shared__ volatile ScalorType buff  [];  
  // buff[tid] = 0.f;
  // __syncthreads();
  
  // buff[tid] = myPoten * 0.5f;
  // sumVectorBlockBuffer_2 (buff);
  // if (threadIdx.x == 0) statistic_nb_buff0[bid] = buff[0];
  // __syncthreads();
  // buff[tid] = myVxx * 0.5f;
  // sumVectorBlockBuffer_2 (buff);
  // if (threadIdx.x == 0) statistic_nb_buff1[bid] = buff[0];
  // __syncthreads();
  // buff[tid] = myVyy * 0.5f;
  // sumVectorBlockBuffer_2 (buff);
  // if (threadIdx.x == 0) statistic_nb_buff2[bid] = buff[0];
  // __syncthreads();
  // buff[tid] = myVzz * 0.5f;
  // sumVectorBlockBuffer_2 (buff);
  // if (threadIdx.x == 0) statistic_nb_buff3[bid] = buff[0];
  // __syncthreads();

  if (ii < numAtom){
    statistic_nb_buff0[ii] = myPoten * 0.5f;
    statistic_nb_buff1[ii] = myVxx * 0.5f;
    statistic_nb_buff2[ii] = myVyy * 0.5f;
    statistic_nb_buff3[ii] = myVzz * 0.5f;
  }
  
  // partialSum = sumVectorBlockBuffer (buff, num);
  // if (sumPartialSum (partialSum, statistic_nb_buff0, &interactionEngine_counter_for_NB,
  // 		     &tmpP) ){
  //   stddata[mdStatisticNonBondedPotential] += tmpP * 0.5f;
  // }
  // __syncthreads();
  // buff[tid] = myVxx;
  // __syncthreads();
  // partialSum = sumVectorBlockBuffer (buff, num);
  // if (sumPartialSum (partialSum, statistic_nb_buff1, &interactionEngine_counter_for_NB_Vxx,
  // 		     &tmpVxx) ){
  //   stddata[mdStatisticVirialXX] += tmpVxx * 0.5f;
  // }
  // __syncthreads();
  // buff[tid] = myVyy;
  // __syncthreads();
  // partialSum = sumVectorBlockBuffer (buff, num);
  // if (sumPartialSum (partialSum, statistic_nb_buff2, &interactionEngine_counter_for_NB_Vyy,
  // 		     &tmpVyy) ){
  //   stddata[mdStatisticVirialYY] += tmpVyy * 0.5f;
  // }
  // __syncthreads();
  // buff[tid] = myVzz;
  // __syncthreads();
  // partialSum = sumVectorBlockBuffer (buff, num);
  // if (sumPartialSum (partialSum, statistic_nb_buff3, &interactionEngine_counter_for_NB_Vzz,
  // 		     &tmpVzz) ){
  //   stddata[mdStatisticVirialZZ] += tmpVzz * 0.5f;
  // }
}










__global__ void calBondInteraction (const IndexType numAtom,
				    const ScalorType * coordx,
				    const ScalorType * coordy, 
				    const ScalorType * coordz,
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
  ScalorType refx = 0.0f, refy = 0.0f, refz = 0.0f;

#ifndef COMPILE_NO_TEX
  refx = coordx[ii];
  refy = coordy[ii];
  refz = coordz[ii];
#else
  refx = tex1Dfetch(global_texRef_interaction_coordx, ii);
  refy = tex1Dfetch(global_texRef_interaction_coordy, ii);
  refz = tex1Dfetch(global_texRef_interaction_coordz, ii);
#endif  
      
  IndexType myNumBond = bdlist.Nbond[ii];
  
  for (IndexType jj = 0; jj < bdlist.listLength; ++jj){
    if (jj == myNumBond) break;
    IndexType targetIdx = bdlist.data[jj * bdlist.stride + ii];
    ScalorType targetx, targety, targetz;

#ifndef COMPILE_NO_TEX
    targetx = coordx[targetIdx];
    targety = coordy[targetIdx];
    targetz = coordz[targetIdx];
#else
    targetx = tex1Dfetch(global_texRef_interaction_coordx, targetIdx);
    targety = tex1Dfetch(global_texRef_interaction_coordy, targetIdx);
    targetz = tex1Dfetch(global_texRef_interaction_coordz, targetIdx);
#endif 
    ScalorType diffx, diffy, diffz;
    diffx = targetx - refx;
    diffy = targety - refy;
    diffz = targetz - refz;
    shortestImage (box, &diffx, &diffy, &diffz);
    ScalorType fx, fy, fz;
    ForceIndexType bondFindex = bdlist.bondIndex[jj * bdlist.stride + ii];
    bondForce (bondForceType[bondFindex],
	       &bondForceParam[bondForceParamPosi[bondFindex]],
	       diffx, diffy, diffz, &fx, &fy, &fz);
    fsumx -= fx;
    fsumy -= fy;
    fsumz -= fz;
  }
  forcx[ii] += fsumx;
  forcy[ii] += fsumy;
  forcz[ii] += fsumz;
}




__global__ void calBondInteraction (const IndexType numAtom,
				    const ScalorType * coordx,
				    const ScalorType * coordy, 
				    const ScalorType * coordz,
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
  // volatile __shared__ ScalorType buff  [MaxThreadsPerBlock];  
  buff[tid] = 0.f;
  // buff[tid+blockDim.x] = 0.f;
  __syncthreads();
  
  ScalorType fsumx = 0.0f;
  ScalorType fsumy = 0.0f;
  ScalorType fsumz = 0.0f;
  IndexType ii = tid + bid * blockDim.x;
  ScalorType myPoten = 0.0f, myVxx = 0.0f, myVyy = 0.0f, myVzz = 0.0f;
  if (ii < numAtom) {
    ScalorType refx = 0.0f, refy = 0.0f, refz = 0.0f;
#ifdef COMPILE_NO_TEX
    refx = coordx[ii];
    refy = coordy[ii];
    refz = coordz[ii];
#else 
    refx = tex1Dfetch(global_texRef_interaction_coordx, ii);
    refy = tex1Dfetch(global_texRef_interaction_coordy, ii);
    refz = tex1Dfetch(global_texRef_interaction_coordz, ii);
#endif
    IndexType myNumBond = bdlist.Nbond[ii];
    for (IndexType jj = 0; jj < bdlist.listLength; ++jj){
      if (jj == myNumBond) break;
      IndexType targetIdx = bdlist.data[jj * bdlist.stride + ii];
      ScalorType targetx, targety, targetz;
#ifdef COMPILE_NO_TEX
      targetx = coordx[targetIdx];
      targety = coordy[targetIdx];
      targetz = coordz[targetIdx];
#else
      targetx = tex1Dfetch(global_texRef_interaction_coordx, targetIdx);
      targety = tex1Dfetch(global_texRef_interaction_coordy, targetIdx);
      targetz = tex1Dfetch(global_texRef_interaction_coordz, targetIdx);
#endif
      ScalorType diffx, diffy, diffz;
      diffx = targetx - refx;
      diffy = targety - refy;
      diffz = targetz - refz;
      shortestImage (box, &diffx, &diffy, &diffz);
      ScalorType fx, fy, fz;
      ForceIndexType bondFindex = bdlist.bondIndex[jj * bdlist.stride + ii];
      ScalorType dp;
      bondForcePoten (bondForceType[bondFindex],
		      &bondForceParam[bondForceParamPosi[bondFindex]],
		      diffx, diffy, diffz, &fx, &fy, &fz, &dp);
      myPoten += dp;
      myVxx -= fx * diffx;
      myVyy -= fy * diffy;
      myVzz -= fz * diffz;
      fsumx -= fx;
      fsumy -= fy;
      fsumz -= fz;
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

  // __syncthreads();
  
  // IndexType num = ((bid+1) * blockDim.x > numAtom) ? 
  //     (blockDim.x - ((bid+1) * blockDim.x) + numAtom) : blockDim.x;
  // ScalorType partialSum;
  // ScalorType tmpP(0.f), tmpVxx(0.f), tmpVyy(0.f), tmpVzz(0.f);
  
  // buff[tid] = myPoten;
  // __syncthreads();
  // partialSum = sumVectorBlockBuffer (buff, num);
  // if (sumPartialSum (partialSum, statistic_b_buff0, &interactionEngine_counter_for_B,
  // 		     &tmpP) ){
  //   stddata[mdStatisticBondedPotential] += tmpP * 0.5f;
  // }
  // __syncthreads();
  // buff[tid] = myVxx;
  // __syncthreads();
  // partialSum = sumVectorBlockBuffer (buff, num);
  // if (sumPartialSum (partialSum, statistic_b_buff1, &interactionEngine_counter_for_B_Vxx,
  // 		     &tmpVxx) ){
  //   stddata[mdStatisticVirialXX] += tmpVxx * 0.5f;
  // }
  // __syncthreads();
  // buff[tid] = myVyy;
  // __syncthreads();
  // partialSum = sumVectorBlockBuffer (buff, num);
  // if (sumPartialSum (partialSum, statistic_b_buff2, &interactionEngine_counter_for_B_Vyy,
  // 		     &tmpVyy) ){
  //   stddata[mdStatisticVirialYY] += tmpVyy * 0.5f;
  // }
  // __syncthreads();
  // buff[tid] = myVzz;
  // __syncthreads();
  // partialSum = sumVectorBlockBuffer (buff, num);
  // if (sumPartialSum (partialSum, statistic_b_buff3, &interactionEngine_counter_for_B_Vzz,
  // 		     &tmpVzz) ){
  //   stddata[mdStatisticVirialZZ] += tmpVzz * 0.5f;
  // }
}
