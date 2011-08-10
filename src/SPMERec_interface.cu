#include "SPMERec.h"

texture<CoordType,  1, cudaReadModeElementType> global_texRef_SPME_QConv;

void SPMERecIk::
freeAll ()
{
  if (malloced){
    cudaFree (vecbx);
    cudaFree (vecby);
    cudaFree (vecbz);
    cudaFree (Q);
    cudaFree (psiF);
    cudaFree (phiF0);
    cudaFree (phiF1);
    cudaFree (phiF2);
    cudaFree (QF);
    cudaFree (QFxPsiF);
    cudaFree (QFxPhiF0);
    cudaFree (QFxPhiF1);
    cudaFree (QFxPhiF2);
    cudaFree (QConvPsi);
    cudaFree (QConvPhi0);
    cudaFree (QConvPhi1);
    cudaFree (QConvPhi2);
    cudaUnbindTexture(global_texRef_SPME_QConv);
    cudaFree (QConvPhi);
    cufftDestroy (planForward);
    cufftDestroy (planBackward);
    malloced = false;
  }
}

SPMERecIk::
SPMERecIk()
    : malloced (false)
{
}

SPMERecIk::
~SPMERecIk ()
{
  freeAll();
}

void SPMERecIk::
calB ()
{
  calBx
      <<<((K.x + meshBlockDim.x - 1) / meshBlockDim.x), meshBlockDim>>> (
	  K,
	  order,
	  vecbx);
  checkCUDAError ("SPMERecIk::calB x");
  calBy
      <<<((K.y + meshBlockDim.x - 1) / meshBlockDim.x), meshBlockDim>>> (
	  K,
	  order,
	  vecby);
  checkCUDAError ("SPMERecIk::calB y");
  calBz
      <<<((K.z + meshBlockDim.x - 1) / meshBlockDim.x), meshBlockDim>>> (
	  K,
	  order,
	  vecbz);
  checkCUDAError ("SPMERecIk::calB z");
}


void SPMERecIk::
buildNeighborList (const MDSystem & sys)
{
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
  nob = (nelehalf  + meshBlockDim.x - 1) / meshBlockDim.x;
  meshGridDim_half = toGridDim (nob);
  atomBlockDim.x = atomNThread;
  nob = (natom + atomBlockDim.x - 1) / atomBlockDim.x;
  atomGridDim = toGridDim (nob);

  cudaMalloc ((void**)&vecbx, sizeof(ScalorType) * K.x);
  cudaMalloc ((void**)&vecby, sizeof(ScalorType) * K.y);
  cudaMalloc ((void**)&vecbz, sizeof(ScalorType) * K.z);
  calB ();
  checkCUDAError ("SPMERecIk::reinit malloc");
  cudaMalloc ((void**)&Q, sizeof(cufftReal) * nele);
  checkCUDAError ("SPMERecIk::reinit malloc");
  cudaMalloc ((void**)&psiF,  sizeof(cufftComplex) * nelehalf);
  cudaMalloc ((void**)&phiF0, sizeof(cufftComplex) * nelehalf);
  cudaMalloc ((void**)&phiF1, sizeof(cufftComplex) * nelehalf);
  cudaMalloc ((void**)&phiF2, sizeof(cufftComplex) * nelehalf);
  checkCUDAError ("SPMERecIk::reinit malloc");
  cudaMalloc ((void**)&QF,  sizeof(cufftComplex) * nelehalf);
  cudaMalloc ((void**)&QFxPsiF,  sizeof(cufftComplex) * nelehalf);
  cudaMalloc ((void**)&QFxPhiF0, sizeof(cufftComplex) * nelehalf);
  cudaMalloc ((void**)&QFxPhiF1, sizeof(cufftComplex) * nelehalf);
  cudaMalloc ((void**)&QFxPhiF2, sizeof(cufftComplex) * nelehalf);
  checkCUDAError ("SPMERecIk::reinit malloc");
  cudaMalloc ((void**)&QConvPsi,  sizeof(cufftReal) * nele);
  cudaMalloc ((void**)&QConvPhi0, sizeof(cufftReal) * nele);
  cudaMalloc ((void**)&QConvPhi1, sizeof(cufftReal) * nele);
  cudaMalloc ((void**)&QConvPhi2, sizeof(cufftReal) * nele);
  cudaMalloc ((void**)&QConvPhi,  sizeof(CoordType) * nele);
  cudaBindTexture(0, global_texRef_SPME_QConv, QConvPhi,
		  sizeof(CoordType) * nele);
  checkCUDAError ("SPMERecIk::reinit malloc");

  cufftPlan3d (&planForward,  K.x, K.y, K.z, CUFFT_R2C);
  cufftPlan3d (&planBackward, K.x, K.y, K.z, CUFFT_C2R);

  malloced = true;

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

  nlist_stride = nele;
  nlist_length = order * order * order * 1;
  cudaMalloc ((void**)&nlist_n, sizeof(IndexType) * nlist_stride);
  cudaMalloc ((void**)&nlist_list, sizeof(IndexType) * nlist_stride * nlist_length);
  checkCUDAError ("SPMERecIk::reinit malloc nlist");

  sum_e.  reinit (nele, NThreadForSum);
  sum_vxx.reinit (nele, NThreadForSum);
  sum_vyy.reinit (nele, NThreadForSum);
  sum_vzz.reinit (nele, NThreadForSum);
  checkCUDAError ("EwaldSumRec::reinit reinit sums");
}

void SPMERecIk::
applyInteraction (MDSystem & sys,
		  MDStatistic * pst,
		  MDTimer * timer)
{
  if (timer != NULL) timer->tic (mdTimeSPMERecCalQ);
  calQ (sys);
  if (timer != NULL) timer->toc (mdTimeSPMERecCalQ);
  
  if (timer != NULL) timer->tic (mdTimeSPMERecFFT);
  cufftExecR2C (planForward, Q, QF);
  checkCUDAError ("SPMERecIk::applyInteraction Q->QF");
  if (timer != NULL) timer->toc (mdTimeSPMERecFFT);
  
  if (timer != NULL) timer->tic (mdTimeSPMERecTimeMatrix);
  IntVectorType N(K);
  N.z = (N.z >> 1) + 1;
  IndexType nelehalf = N.x * N.y * N.z;
  ScalorType sizei = 1./(K.x*K.y*K.z);
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
  checkCUDAError ("SPMERecIk::applyInteraction timeQFPhiF");
  if (timer != NULL) timer->toc (mdTimeSPMERecTimeMatrix);

  if (timer != NULL) timer->tic (mdTimeSPMERecFFT);
  cufftExecC2R (planBackward, QFxPhiF0, QConvPhi0);
  cufftExecC2R (planBackward, QFxPhiF1, QConvPhi1);
  cufftExecC2R (planBackward, QFxPhiF2, QConvPhi2);
  checkCUDAError ("SPMERecIk::applyInteraction QFxPhiF->QConvPhi");
  if (timer != NULL) timer->toc (mdTimeSPMERecFFT);

  if (timer != NULL) timer->tic (mdTimeSPMERecForce);
  assembleQConvPhi
      <<<meshGridDim, meshBlockDim>>> (
	  QConvPhi0,
	  QConvPhi1,
	  QConvPhi2,
	  QConvPhi,
	  K.x*K.y*K.z);
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
    timeQFPsiF
	<<<meshGridDim_half, meshBlockDim>>> (
	    QF,
	    psiF,
	    QFxPsiF,
	    nelehalf,
	    sizei);
    checkCUDAError ("SPMERecIk::applyInteraction time QF PhiF");
    if (timer != NULL) timer->toc (mdTimeSPMERecTimeMatrix);
    if (timer != NULL) timer->tic (mdTimeSPMERecFFT);
    cufftExecC2R (planBackward, QFxPsiF, QConvPsi);
    checkCUDAError ("SPMERecIk::applyInteraction QFxPsiF->QConvPsi");
    if (timer != NULL) timer->toc (mdTimeSPMERecFFT);
    if (timer != NULL) timer->tic (mdTimeSPMERecEnergy);
    calEnergy
	<<<meshGridDim, meshBlockDim>>> (
	    Q,
	    QConvPsi,
	    sum_e.buff,
	    K.x * K.y * K.z);
    checkCUDAError ("SPMERecIk::applyInteraction cal energy");
    sum_e.  sumBuffAdd(pst->ddata, mdStatisticNonBondedPotential);
    checkCUDAError ("SPMERecIk::applyInteraction sum energy");
    if (timer != NULL) timer->toc (mdTimeSPMERecEnergy);
  }
}



void SPMERecIk::
calQ (const MDSystem & sys)
{
//
// fast algorithm !!!!
//
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
  // FILE * fp = fopen ("tmpQ.out", "w");
  // for (unsigned i = 0; i < K.x * K.y * K.z; ++i){
  //   fprintf (fp, "%.12e\n", Q[i]);
  // }
  // fclose (fp);
}
  
void SPMERecIk::
calV()
{
  volume =
      vecA.xx * (vecA.yy*vecA.zz - vecA.zy*vecA.yz) - 
      vecA.xy * (vecA.yx*vecA.zz - vecA.zx*vecA.yz) +
      vecA.xz * (vecA.yx*vecA.zy - vecA.zx*vecA.yy);
}
  
void SPMERecIk::
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
	  Mnx[jj] = BSplineValue (order, value.x);
	  Mny[jj] = BSplineValue (order, value.y);
	  Mnz[jj] = BSplineValue (order, value.z);
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
