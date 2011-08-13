#include "SPMERec.h"

bool global_texRef_SPME_QConv_binded;
texture<CoordType,  1, cudaReadModeElementType> global_texRef_SPME_QConv;

SPMERec::
SPMERec ()
    : malloced (false),
      enable_nlist (true)
{
  global_texRef_SPME_QConv_binded = false;

  IndexType * d_arch, cudaArch;
  cudaMalloc ((void**)&d_arch, sizeof(IndexType));
  getArch <<<1,1>>> (d_arch);
  cudaMemcpy (&cudaArch, d_arch, sizeof(IndexType), cudaMemcpyDeviceToHost);
  cudaFree (d_arch);
  printf ("cudaArch is %d\n", cudaArch);
  if (cudaArch == 200) enable_nlist = false;
}

SPMERec::
~SPMERec ()
{
  freeAll();
}

void SPMERec::
freeAll ()
{
  if (malloced){
    cudaFree (vecbx);
    cudaFree (vecby);
    cudaFree (vecbz);
    cudaFree (Q);
    if (enable_nlist){
      cudaFree (nlist_n);
      cudaFree (nlist_list);
    }
    malloced = false;
  }
}

void SPMERec::
calB ()
{
  calBx
      <<<((K.x + meshBlockDim.x - 1) / meshBlockDim.x), meshBlockDim>>> (
	  K,
	  order,
	  vecbx);
  checkCUDAError ("SPMERec::calB x");
  calBy
      <<<((K.y + meshBlockDim.x - 1) / meshBlockDim.x), meshBlockDim>>> (
	  K,
	  order,
	  vecby);
  checkCUDAError ("SPMERec::calB y");
  calBz
      <<<((K.z + meshBlockDim.x - 1) / meshBlockDim.x), meshBlockDim>>> (
	  K,
	  order,
	  vecbz);
  checkCUDAError ("SPMERec::calB z");
}

void  SPMERec::
reinit (const MatrixType & vecA_,
	const IntVectorType & K_,
	const IndexType & order_,
	const ScalorType & beta_,
	const IndexType & natom,
	const IndexType & meshNThread,
	const IndexType & atomNThread)
{
  freeAll();
  vecA = vecA_;
  K = K_;
  order = order_;
  beta = beta_;
  calV();
  calAStar();
  
  IndexType nele = K.x * K.y * K.z;
  IntVectorType N(K);
  N.z = (N.z >> 1) + 1;
  IndexType nelehalf = N.x * N.y * N.z;

  IndexType nob;
  meshBlockDim.x = meshNThread;
  nob = (nele  + meshBlockDim.x - 1) / meshBlockDim.x;
  meshGridDim = toGridDim (nob);
  meshGridDim.y = 8;
  meshGridDim.x /= meshGridDim.y;
  meshGridDim.x ++;
  nob = (nelehalf  + meshBlockDim.x - 1) / meshBlockDim.x;
  meshGridDim_half = toGridDim (nob);
  meshGridDim_half.y = 8;
  meshGridDim_half.x /= meshGridDim_half.y;
  meshGridDim_half.x ++;
  atomBlockDim.x = atomNThread;
  nob = (natom + atomBlockDim.x - 1) / atomBlockDim.x;
  atomGridDim = toGridDim (nob);
  atomGridDim.y = 8;
  atomGridDim.x /= atomGridDim.y;
  atomGridDim.x ++;

  cudaMalloc ((void**)&vecbx, sizeof(ScalorType) * K.x);
  cudaMalloc ((void**)&vecby, sizeof(ScalorType) * K.y);
  cudaMalloc ((void**)&vecbz, sizeof(ScalorType) * K.z);
  calB ();
  cudaMalloc ((void**)&Q, sizeof(cufftReal) * nele);
  checkCUDAError ("SPMERec::reinit malloc");  
  // printf ("vecbx is %f\n", vecbx[0]);
  
  if (enable_nlist) {
    nlist_stride = nele;
    ScalorType rho = natom / volume;
    ScalorType vcell = volume * ScalorType (order * order * order) / ScalorType (K.x * K.y * K.z);
    nlist_length = (rho * vcell * 1.5f + 15.);
    printf ("# spme nlist length is %d\n", nlist_length);
    cudaMalloc ((void**)&nlist_n, sizeof(IndexType) * nlist_stride);
    cudaMalloc ((void**)&nlist_list, sizeof(IndexType) * nlist_stride * nlist_length);
    checkCUDAError ("SPMERec::reinit malloc nlist");
  }
  
  malloced = true;
}


void SPMERec::
calQ (const MDSystem & sys,
      MDTimer * timer)
{
  if (enable_nlist){
//
// fast algorithm of nlist !!!!
//  
    if (timer != NULL) timer->tic (mdTimeSPMERecMeshNeighborList);
    initMeshNeighborList
	<<<meshGridDim, meshBlockDim>>>(
	    K,
	    vecAStar,
	    nlist_n,
	    nlist_list,
	    nlist_stride,
	    nlist_length);
    buildMeshNeighborList
	<<<atomGridDim, atomBlockDim>>> (
	    K,
	    vecAStar,
	    order,
	    sys.ddata.coord,
	    sys.ddata.numAtom,
	    nlist_n,
	    nlist_list,
	    nlist_stride,
	    nlist_length,
	    err.ptr_de);
    checkCUDAError ("SPMERecIk::calQ buildNeighborList");
    err.check ("SPMERecIk::calQ buildNeighborList");
    if (timer != NULL) timer->toc (mdTimeSPMERecMeshNeighborList);
    if (timer != NULL) timer->tic (mdTimeSPMECalQFromNList);
    calQMat
	<<<meshGridDim, meshBlockDim>>> (
	    K,
	    vecAStar,
	    order,
	    nlist_n,
	    nlist_list,
	    nlist_stride,
	    Q);
    checkCUDAError ("SPMERecIk::calQ calQ");
    if (timer != NULL) timer->toc (mdTimeSPMECalQFromNList);
  }
  else {
    setValue
    	<<<meshGridDim, meshBlockDim>>> (
    	    Q,
    	    K.x * K.y * K.z,
    	    ScalorType (0.f));
    calQMat
    	<<<atomGridDim, atomBlockDim>>> (
    	    K, 
    	    vecAStar,
    	    order,
    	    sys.ddata.coord,
    	    sys.ddata.charge,
    	    sys.ddata.numAtom,
    	    Q,
    	    err.ptr_de);
    checkCUDAError ("SPMERecIk::calQ calQ atomicAdd");
  }
  
  // // int nele = K.x * K.y * K.z;
  // // cufftReal * see = (cufftReal * )malloc (sizeof(cufftReal) * nele);
  // cudaMemcpy (see, Q, sizeof(cufftReal) * nele, cudaMemcpyDeviceToHost);
  // cudaMemcpy (see, Q, sizeof(cufftReal) * nele, cudaMemcpyDeviceToHost);
  // free (see);
}

void SPMERec::
calV()
{
  volume =
      vecA.xx * (vecA.yy*vecA.zz - vecA.zy*vecA.yz) - 
      vecA.xy * (vecA.yx*vecA.zz - vecA.zx*vecA.yz) +
      vecA.xz * (vecA.yx*vecA.zy - vecA.zx*vecA.yy);
}
  
void SPMERec::
calAStar ()
{
  ScalorType volumei = ScalorType(1.) / volume;
  vecAStar.xx = ( vecA.yy*vecA.zz - vecA.zy*vecA.yz) * volumei;
  vecAStar.yy = ( vecA.xx*vecA.zz - vecA.zx*vecA.xz) * volumei;
  vecAStar.zz = ( vecA.xx*vecA.yy - vecA.yx*vecA.xy) * volumei;
  vecAStar.yx = (-vecA.yx*vecA.zz + vecA.zx*vecA.yz) * volumei;
  vecAStar.zx = ( vecA.yx*vecA.zy - vecA.zx*vecA.yy) * volumei;
  vecAStar.xy = (-vecA.xy*vecA.zz + vecA.zy*vecA.xz) * volumei;
  vecAStar.zy = (-vecA.xx*vecA.zy + vecA.zx*vecA.xy) * volumei;
  vecAStar.xz = ( vecA.xy*vecA.yz - vecA.yy*vecA.xz) * volumei;
  vecAStar.yz = (-vecA.xx*vecA.yz + vecA.yx*vecA.xz) * volumei;
}


void SPMERecIk::
freeAll ()
{
  if (malloced){
    cudaFree (psiF);
    cudaFree (phiF0);
    cudaFree (phiF1);
    cudaFree (phiF2);
    if (fftOutOfPlace) {
      cudaFree (QF);
      cudaFree (QFxPsiF);
      cudaFree (QFxPhiF0);
      cudaFree (QFxPhiF1);
      cudaFree (QFxPhiF2);
    }
    cudaFree (QConvPsi);
    cudaFree (QConvPhi0);
    cudaFree (QConvPhi1);
    cudaFree (QConvPhi2);
    if (global_texRef_SPME_QConv_binded){
      cudaUnbindTexture(global_texRef_SPME_QConv);
    }
    cudaFree (QConvPhi);
    cufftDestroy (planForward);
    cufftDestroy (planBackward);
    malloced = false;
  }
}

SPMERecIk::
SPMERecIk()
    : malloced (false),
      fftOutOfPlace (true)
{
}

SPMERecIk::
~SPMERecIk ()
{
  freeAll();
}

void  SPMERecIk::
reinit (const MatrixType & vecA_,
	const IntVectorType & K_,
	const IndexType & order_,
	const ScalorType & beta_,
	const IndexType & natom,
	const IndexType & meshNThread,
	const IndexType & atomNThread)
{
  SPMERec::reinit (vecA_, K_, order_, beta_, natom, meshNThread, atomNThread);
  freeAll();
  
  IndexType nele = K.x * K.y * K.z;
  IntVectorType N(K);
  N.z = (N.z >> 1) + 1;
  IndexType nelehalf = N.x * N.y * N.z;
  KPadding = K;
  if (!fftOutOfPlace){
    KPadding.z = (N.z << 1);
  }

  cudaMalloc ((void**)&QF,  sizeof(cufftComplex) * nelehalf);
  checkCUDAError ("SPMERecIk::reinit malloc");
  if (fftOutOfPlace){
    cudaMalloc ((void**)&QFxPsiF,  sizeof(cufftComplex) * nelehalf);
    cudaMalloc ((void**)&QFxPhiF0, sizeof(cufftComplex) * nelehalf);
    cudaMalloc ((void**)&QFxPhiF1, sizeof(cufftComplex) * nelehalf);
    cudaMalloc ((void**)&QFxPhiF2, sizeof(cufftComplex) * nelehalf);
    checkCUDAError ("SPMERecIk::reinit malloc");
    cudaMalloc ((void**)&QConvPsi,  sizeof(cufftReal) * nele);
    cudaMalloc ((void**)&QConvPhi0, sizeof(cufftReal) * nele);
    cudaMalloc ((void**)&QConvPhi1, sizeof(cufftReal) * nele);
    cudaMalloc ((void**)&QConvPhi2, sizeof(cufftReal) * nele);
    checkCUDAError ("SPMERecIk::reinit malloc");
  }
  else {
    int nelePading = KPadding.x * KPadding.y * KPadding.z;
    IndexType nob = (nelePading  + meshBlockDim.x - 1) / meshBlockDim.x;
    dim3 meshGridDim_tmp = toGridDim (nob);
    meshGridDim_tmp.y = 8;
    meshGridDim_tmp.x /= meshGridDim_half.y;
    meshGridDim_tmp.x ++;
    cudaMalloc ((void**)&QConvPsi,  sizeof(cufftReal) * nelePading);
    cudaMalloc ((void**)&QConvPhi0, sizeof(cufftReal) * nelePading);
    cudaMalloc ((void**)&QConvPhi1, sizeof(cufftReal) * nelePading);
    cudaMalloc ((void**)&QConvPhi2, sizeof(cufftReal) * nelePading);    
    setValue <<<meshGridDim_tmp.x, meshBlockDim>>> (QConvPsi, nelePading, ScalorType(0.f));
    setValue <<<meshGridDim_tmp.x, meshBlockDim>>> (QConvPhi0, nelePading, ScalorType(0.f));
    setValue <<<meshGridDim_tmp.x, meshBlockDim>>> (QConvPhi1, nelePading, ScalorType(0.f));
    setValue <<<meshGridDim_tmp.x, meshBlockDim>>> (QConvPhi2, nelePading, ScalorType(0.f));
  }

  cudaMalloc ((void**)&psiF,  sizeof(cufftComplex) * nelehalf);
  cudaMalloc ((void**)&phiF0, sizeof(cufftComplex) * nelehalf);
  cudaMalloc ((void**)&phiF1, sizeof(cufftComplex) * nelehalf);
  cudaMalloc ((void**)&phiF2, sizeof(cufftComplex) * nelehalf);
  checkCUDAError ("SPMERecIk::reinit malloc");
  
  cudaMalloc ((void**)&QConvPhi,  sizeof(CoordType) * nele);
  cudaBindTexture(0, global_texRef_SPME_QConv, QConvPhi,
		  sizeof(CoordType) * nele);
  global_texRef_SPME_QConv_binded = true;
  checkCUDAError ("SPMERecIk::reinit malloc");

  cufftPlan3d (&planForward,  K.x, K.y, K.z, CUFFT_R2C);
  cufftPlan3d (&planBackward, K.x, K.y, K.z, CUFFT_C2R);
  // cufftResult r = cufftSetCompatibilityMode (planForward, CUFFT_COMPATIBILITY_FFTW_PADDING);
  // if (r != CUFFT_SUCCESS) exit(1);
  // cufftSetCompatibilityMode (planBackward, CUFFT_COMPATIBILITY_FFTW_PADDING);
  // cufftSetCompatibilityMode (planForward, CUFFT_COMPATIBILITY_FFTW_ALL);
  // cufftSetCompatibilityMode (planBackward, CUFFT_COMPATIBILITY_FFTW_ALL);
  
  calPsiFPhiF
      <<<meshGridDim_half, meshBlockDim>>> (
	  K,
	  vecAStar,
	  beta,
	  volume,
	  vecbx,
	  vecby,
	  vecbz,
	  psiF,
	  phiF0,
	  phiF1,
	  phiF2);
  checkCUDAError ("SPMERecIk::reinit calPsiFPhiF");

  malloced = true;

  sum_e.  reinit (nele, NThreadForSum);
  // sum_vxx.reinit (nele, NThreadForSum);
  // sum_vyy.reinit (nele, NThreadForSum);
  // sum_vzz.reinit (nele, NThreadForSum);
  checkCUDAError ("EwaldSumRec::reinit reinit sums");
}

void SPMERecIk::
applyInteraction (MDSystem & sys,
		  MDStatistic * pst,
		  MDTimer * timer)
{
  // int nele = K.x * K.y * K.z;
  IntVectorType N(K);
  N.z = (N.z >> 1) + 1;
  IndexType nelehalf = N.x * N.y * N.z;
  // size_t see_size = KPadding.x * KPadding.y * KPadding.z * sizeof(cufftReal);
  // size_t cpl_size = nelehalf * sizeof(cufftComplex);
  // cufftReal* see = (cufftReal *) malloc (see_size);
  // cufftComplex* cpl = (cufftComplex *) malloc (cpl_size);
  if (timer != NULL) timer->tic (mdTimeSPMERecCalQ);
  calQ (sys, timer);
  // cufftHandle planBackward1, planForward1;
  // cufftPlan3d (&planForward1,  KPadding.x, KPadding.y, KPadding.z, CUFFT_R2C);
  // cufftPlan3d (&planBackward1,  KPadding.x, KPadding.y, KPadding.z, CUFFT_C2R);
  // cufftPlan1d (&planForward1,   10, CUFFT_R2C, 1);
  // cufftPlan1d (&planBackward1,   10, CUFFT_C2R, 1);
  // cudaMemcpy (see, Q, see_size, cudaMemcpyDeviceToHost);
  // cudaMemcpy (cpl, phiF0, cpl_size, cudaMemcpyDeviceToHost);
  // cufftExecR2C (planForward,  (cufftReal*)phiF0, phiF0);
  // cudaMemcpy (cpl, phiF0, cpl_size, cudaMemcpyDeviceToHost);
  // cufftExecC2R (planBackward, phiF0, (cufftReal*)phiF0);
  // cudaMemcpy (cpl, phiF0, cpl_size, cudaMemcpyDeviceToHost);
  // cudaMemcpy (cpl, phiF1, cpl_size, cudaMemcpyDeviceToHost);
  // cudaMemcpy (cpl, phiF2, cpl_size, cudaMemcpyDeviceToHost);
  if (timer != NULL) timer->toc (mdTimeSPMERecCalQ);
  
  if (timer != NULL) timer->tic (mdTimeSPMERecFFT);
  cufftExecR2C (planForward, Q, QF);
  checkCUDAError ("SPMERecIk::applyInteraction Q->QF");
  if (timer != NULL) timer->toc (mdTimeSPMERecFFT);
  
  if (timer != NULL) timer->tic (mdTimeSPMERecTimeMatrix);
  ScalorType sizei = 1./(K.x*K.y*K.z);
  if (fftOutOfPlace){
    timeQFPhiF
	<<<meshGridDim_half, meshBlockDim>>> (
	    QF,
	    phiF0,
	    phiF1,
	    phiF2,
	    QFxPhiF0,
	    QFxPhiF1,
	    QFxPhiF2,
	    nelehalf,
	    sizei);
  }
  else {
    timeQFPhiF
	<<<meshGridDim_half, meshBlockDim>>> (
	    QF,
	    phiF0,
	    phiF1,
	    phiF2,
	    QConvPhi0,
	    QConvPhi1,
	    QConvPhi2,
	    nelehalf,
	    sizei);
  }
  checkCUDAError ("SPMERecIk::applyInteraction timeQFPhiF");
  if (timer != NULL) timer->toc (mdTimeSPMERecTimeMatrix);

  if (timer != NULL) timer->tic (mdTimeSPMERecFFT);
  if (fftOutOfPlace){
    cufftExecC2R (planBackward, QFxPhiF0, QConvPhi0);
    cufftExecC2R (planBackward, QFxPhiF1, QConvPhi1);
    cufftExecC2R (planBackward, QFxPhiF2, QConvPhi2);
  }
  else {
    cufftExecC2R (planBackward, (cufftComplex *)QConvPhi0, QConvPhi0);
    cufftExecC2R (planBackward, (cufftComplex *)QConvPhi1, QConvPhi1);
    cufftExecC2R (planBackward, (cufftComplex *)QConvPhi2, QConvPhi2);
  }
  checkCUDAError ("SPMERecIk::applyInteraction QFxPhiF->QConvPhi");
  if (timer != NULL) timer->toc (mdTimeSPMERecFFT);
  // cudaMemcpy (see, QConvPhi0, see_size, cudaMemcpyDeviceToHost);
  // cudaMemcpy (see, QConvPhi1, see_size, cudaMemcpyDeviceToHost);
  // cudaMemcpy (see, QConvPhi2, see_size, cudaMemcpyDeviceToHost);

  if (timer != NULL) timer->tic (mdTimeSPMERecForce);
  if (fftOutOfPlace){
    assembleQConvPhi
	<<<meshGridDim, meshBlockDim>>> (
	    QConvPhi0,
	    QConvPhi1,
	    QConvPhi2,
	    QConvPhi,
	    K.x*K.y*K.z);
  }
  else {
    assembleQConvPhi
	<<<meshGridDim, meshBlockDim>>> (
	    K, KPadding,
	    QConvPhi0,
	    QConvPhi1,
	    QConvPhi2,
	    QConvPhi,
	    K.x*K.y*K.z);
  }
  calForce
      <<<atomGridDim, atomBlockDim>>> (
	  K,
	  vecAStar,
	  order,
	  sys.ddata.coord,
	  sys.ddata.charge,
	  sys.ddata.numAtom,
	  QConvPhi,
	  // QConvPhi0,
	  // QConvPhi1,
	  // QConvPhi2,
	  sys.ddata.forcx,
	  sys.ddata.forcy,
	  sys.ddata.forcz,
	  err.ptr_de);
  checkCUDAError ("SPMERecIk::applyInteraction calForce");
  err.check ("SPMERecIk::applyInteraction calForce");
  if (timer != NULL) timer->toc (mdTimeSPMERecForce);

  if (pst != NULL){
    if (timer != NULL) timer->tic (mdTimeSPMERecTimeMatrix);
    if (fftOutOfPlace){
      timeQFPsiF
	  <<<meshGridDim_half, meshBlockDim>>> (
	      QF,
	      psiF,
	      QFxPsiF,
	      nelehalf,
	      sizei);
    }
    else {
      timeQFPsiF
	  <<<meshGridDim_half, meshBlockDim>>> (
	      QF,
	      psiF,
	      QConvPsi,
	      nelehalf,
	      sizei);
      FILE * fp = fopen ("tmp.tmp", "w");
      for (int i = 0; i < nelehalf; ++i){
	fprintf (fp, "%e %e\n", QConvPsi[i*2], QConvPsi[i*2+1]);
      }
      fclose (fp);
    }
    checkCUDAError ("SPMERecIk::applyInteraction time QF PhiF");
    if (timer != NULL) timer->toc (mdTimeSPMERecTimeMatrix);
    if (timer != NULL) timer->tic (mdTimeSPMERecFFT);
    if (fftOutOfPlace){
      cufftExecC2R (planBackward, QFxPsiF, QConvPsi);
    }
    else {
      cufftExecC2R (planBackward, (cufftComplex*)QConvPsi, QConvPsi);
    }
    checkCUDAError ("SPMERecIk::applyInteraction QFxPsiF->QConvPsi");
    if (timer != NULL) timer->toc (mdTimeSPMERecFFT);
    if (timer != NULL) timer->tic (mdTimeSPMERecEnergy);
    if (fftOutOfPlace){
      calEnergy
	  <<<meshGridDim, meshBlockDim>>> (
	      Q,
	      QConvPsi,
	      sum_e.buff,
	      K.x * K.y * K.z);
    }
    else {
      calEnergy
	  <<<meshGridDim, meshBlockDim>>> (
	      Q,
	      QConvPsi,
	      K, KPadding,
	      sum_e.buff,
	      K.x * K.y * K.z);
    }
    checkCUDAError ("SPMERecIk::applyInteraction cal energy");
    sum_e.  sumBuffAdd(pst->ddata, mdStatisticNonBondedPotential);
    checkCUDAError ("SPMERecIk::applyInteraction sum energy");
    if (timer != NULL) timer->toc (mdTimeSPMERecEnergy);
  }
}




#include "Ewald.h"
#include "CardinalBspline.h"

// atom grid and block
__global__ void
calForce (const IntVectorType K,
	  const MatrixType vecAStar,
	  const IndexType order,
	  const CoordType * coord,
	  const ScalorType * charge,
	  const IndexType natom,
	  const cufftReal * QConvPhi0,
	  const cufftReal * QConvPhi1,
	  const cufftReal * QConvPhi2,
	  ScalorType * forcx,
	  ScalorType * forcy,
	  ScalorType * forcz,
	  mdError_t * ptr_de )
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;

  if (ii < natom){
    CoordType my_coord (coord[ii]);
    VectorType uu;
    uu.x = K.x * (vecAStar.xx * my_coord.x + vecAStar.xy * my_coord.y + vecAStar.xz *my_coord.z);
    uu.y = K.y * (vecAStar.yx * my_coord.x + vecAStar.yy * my_coord.y + vecAStar.yz *my_coord.z);
    uu.z = K.z * (vecAStar.zx * my_coord.x + vecAStar.zy * my_coord.y + vecAStar.zz *my_coord.z);
    if      (uu.x <  0  ) uu.x += K.x;
    else if (uu.x >= K.x) uu.x -= K.x;
    if      (uu.y <  0  ) uu.y += K.y;
    else if (uu.y >= K.y) uu.y -= K.y;
    if      (uu.z <  0  ) uu.z += K.z;
    else if (uu.z >= K.z) uu.z -= K.z;
    IntVectorType meshIdx;
    meshIdx.x = int(uu.x);
    meshIdx.y = int(uu.y);
    meshIdx.z = int(uu.z);
    if (meshIdx.x >= K.x || meshIdx.x < 0 ||
	meshIdx.y >= K.y || meshIdx.y < 0 ||
	meshIdx.z >= K.z || meshIdx.z < 0 ){
      *ptr_de = mdErrorOverFlowMeshIdx;
    }
    else {
      ScalorType Mnx[MaxSPMEOrder];
      ScalorType Mny[MaxSPMEOrder];
      ScalorType Mnz[MaxSPMEOrder];
      {
	VectorType value;
	value.x = uu.x - meshIdx.x;
	value.y = uu.y - meshIdx.y;
	value.z = uu.z - meshIdx.z;
	for (IndexType jj = 0; jj < order; ++jj){
	  Mnx[jj] = BSplineValue (order, double(value.x));
	  Mny[jj] = BSplineValue (order, double(value.y));
	  Mnz[jj] = BSplineValue (order, double(value.z));
	  value.x += 1.;
	  value.y += 1.;
	  value.z += 1.;
	}
      }
      ScalorType mycharge = charge[ii];
      VectorType fsum;
      fsum.x = fsum.y = fsum.z = ScalorType(0.f);
      IntVectorType myMeshIdx;
      IntVectorType dk;
      for (dk.x = 0; dk.x < order; ++dk.x){
	myMeshIdx.x = meshIdx.x - dk.x;
	if (myMeshIdx.x < 0) myMeshIdx.x += K.x;
	for (dk.y = 0; dk.y < order; ++dk.y){
	  myMeshIdx.y = meshIdx.y - dk.y;
	  if (myMeshIdx.y < 0) myMeshIdx.y += K.y;
	  for (dk.z = 0; dk.z < order; ++dk.z){
	    myMeshIdx.z = meshIdx.z - dk.z;
	    if (myMeshIdx.z < 0) myMeshIdx.z += K.z;
	    //
	    // possible improvement!!!!
	    //
	    IndexType index = index3to1 (myMeshIdx, K);
	    ScalorType myP = Mnx[dk.x] * Mny[dk.y] * Mnz[dk.z];
	    //
	    // possible improvement of reading mem!!!!
	    //
	    fsum.x += myP * QConvPhi0[index];
	    fsum.y += myP * QConvPhi1[index];
	    fsum.z += myP * QConvPhi2[index];
	  }
	}
      }
      forcx[ii] += mycharge * fsum.x;
      forcy[ii] += mycharge * fsum.y;
      forcz[ii] += mycharge * fsum.z;
    }
  }
}


// atom grid and block
__global__ void
calForce (const IntVectorType K,
	  const MatrixType vecAStar,
	  const IndexType order,
	  const CoordType * coord,
	  const ScalorType * charge,
	  const IndexType natom,
	  const CoordType * QConvPhi,
	  ScalorType * forcx,
	  ScalorType * forcy,
	  ScalorType * forcz,
	  mdError_t * ptr_de )
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;

  if (ii < natom){
    CoordType my_coord (coord[ii]);
    VectorType uu;
    uu.x = K.x * (vecAStar.xx * my_coord.x + vecAStar.xy * my_coord.y + vecAStar.xz *my_coord.z);
    uu.y = K.y * (vecAStar.yx * my_coord.x + vecAStar.yy * my_coord.y + vecAStar.yz *my_coord.z);
    uu.z = K.z * (vecAStar.zx * my_coord.x + vecAStar.zy * my_coord.y + vecAStar.zz *my_coord.z);
    if      (uu.x <  0  ) uu.x += K.x;
    else if (uu.x >= K.x) uu.x -= K.x;
    if      (uu.y <  0  ) uu.y += K.y;
    else if (uu.y >= K.y) uu.y -= K.y;
    if      (uu.z <  0  ) uu.z += K.z;
    else if (uu.z >= K.z) uu.z -= K.z;
    IntVectorType meshIdx;
    meshIdx.x = int(uu.x);
    meshIdx.y = int(uu.y);
    meshIdx.z = int(uu.z);
    if (meshIdx.x >= K.x || meshIdx.x < 0 ||
	meshIdx.y >= K.y || meshIdx.y < 0 ||
	meshIdx.z >= K.z || meshIdx.z < 0 ){
      *ptr_de = mdErrorOverFlowMeshIdx;
    }
    else {
      ScalorType Mnx[MaxSPMEOrder];
      ScalorType Mny[MaxSPMEOrder];
      ScalorType Mnz[MaxSPMEOrder];
      {
	VectorType value;
	value.x = uu.x - meshIdx.x;
	value.y = uu.y - meshIdx.y;
	value.z = uu.z - meshIdx.z;
	for (IndexType jj = 0; jj < order; ++jj){
	  Mnx[jj] = BSplineValue (order, double(value.x));
	  Mny[jj] = BSplineValue (order, double(value.y));
	  Mnz[jj] = BSplineValue (order, double(value.z));
	  value.x += 1.;
	  value.y += 1.;
	  value.z += 1.;
	}
      }
      ScalorType mycharge = charge[ii];
      VectorType fsum;
      VectorType fsum_small;
      fsum.x = fsum.y = fsum.z = ScalorType(0.f);
      fsum_small.x = fsum_small.y = fsum_small.z = ScalorType(0.f);
      IntVectorType myMeshIdx;
      IntVectorType dk;
      for (dk.x = 0; dk.x < order; ++dk.x){
	myMeshIdx.x = meshIdx.x - dk.x;
	if (myMeshIdx.x < 0) myMeshIdx.x += K.x;
	for (dk.y = 0; dk.y < order; ++dk.y){
	  myMeshIdx.y = meshIdx.y - dk.y;
	  if (myMeshIdx.y < 0) myMeshIdx.y += K.y;
	  for (dk.z = 0; dk.z < order; ++dk.z){
	    myMeshIdx.z = meshIdx.z - dk.z;
	    if (myMeshIdx.z < 0) myMeshIdx.z += K.z;
	    //
	    // possible improvement!!!!
	    //
	    IndexType index = index3to1 (myMeshIdx, K);
	    ScalorType myP = Mnx[dk.x] * Mny[dk.y] * Mnz[dk.z];
	    CoordType tmpQConvPhi = tex1Dfetch (global_texRef_SPME_QConv, index);
	    // CoordType tmpQConvPhi = QConvPhi[index];
	    // fsum.x += myP * tmpQConvPhi.x;
	    // fsum.y += myP * tmpQConvPhi.y;
	    // fsum.z += myP * tmpQConvPhi.z;
	    ScalorType tmp;
	    tmp = myP * tmpQConvPhi.x;
	    if (fabs(tmp) > 1e-2) fsum.x += tmp;
	    else fsum_small.x += tmp;
	    tmp = myP * tmpQConvPhi.y;
	    if (fabs(tmp) > 1e-2) fsum.y += tmp;
	    else fsum_small.y += tmp;
	    tmp = myP * tmpQConvPhi.z;
	    if (fabs(tmp) > 1e-2) fsum.z += tmp;
	    else fsum_small.z += tmp;
	  }
	}
      }
      forcx[ii] += mycharge * (fsum.x + fsum_small.x);
      forcy[ii] += mycharge * (fsum.y + fsum_small.y);
      forcz[ii] += mycharge * (fsum.z + fsum_small.z);
    }
  }
}

