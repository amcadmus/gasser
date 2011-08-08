#include "SPMERec.h"

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
  cal_Bx
      <<<((K.x + meshBlockDim.x - 1) / meshBlockDim.x), meshBlockDim>>> (
	  K,
	  order,
	  vecbx);
  checkCUDAError ("SPMERecIk::calB x");
  cal_By
      <<<((K.y + meshBlockDim.x - 1) / meshBlockDim.x), meshBlockDim>>> (
	  K,
	  order,
	  vecby);
  checkCUDAError ("SPMERecIk::calB y");
  cal_Bz
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
  checkCUDAError ("SPMERecIk::reinit malloc");

  cufftPlan3d (&planForward,  K.x, K.y, K.z, CUFFT_R2C);
  cufftPlan3d (&planBackward, K.x, K.y, K.z, CUFFT_C2R);

  malloced = true;

  cal_PsiFPhiF
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
  checkCUDAError ("SPMERecIk::reinit cal_PsiFPhiF");

  nlist_stride = nele;
  nlist_length = order * order * order * 2;
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
  calQ (sys);
  
  cufftExecR2C (planForward, Q, QF);
  checkCUDAError ("SPMERecIk::applyInteraction Q->QF");
  
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

  cufftExecC2R (planBackward, QFxPhiF0, QConvPhi0);
  cufftExecC2R (planBackward, QFxPhiF1, QConvPhi1);
  cufftExecC2R (planBackward, QFxPhiF2, QConvPhi2);
  checkCUDAError ("SPMERecIk::applyInteraction QFxPhiF->QConvPhi");

  calForce
      <<<atomGridDim, atomBlockDim>>> (
	  K,
	  vecAStar,
	  order,
	  sys.ddata.coord,
	  sys.ddata.charge,
	  sys.ddata.numAtom,
	  QConvPhi0,
	  QConvPhi1,
	  QConvPhi2,
	  sys.ddata.forcx,
	  sys.ddata.forcy,
	  sys.ddata.forcz,
	  err.ptr_de);
  checkCUDAError ("SPMERecIk::applyInteraction calForce");
  err.check ("SPMERecIk::applyInteraction calForce");

  if (pst != NULL){
    timeQFPsiF
	<<<meshGridDim_half, meshBlockDim>>> (
	    QF,
	    psiF,
	    QFxPsiF,
	    nelehalf,
	    sizei);
    checkCUDAError ("SPMERecIk::applyInteraction time QF PhiF");
    cufftExecC2R (planBackward, QFxPsiF, QConvPsi);
    checkCUDAError ("SPMERecIk::applyInteraction QFxPsiF->QConvPsi");
    calEnergy
	<<<meshGridDim, meshBlockDim>>> (
	    Q,
	    QConvPsi,
	    sum_e.buff,
	    K.x * K.y * K.z);
    checkCUDAError ("SPMERecIk::applyInteraction cal energy");
    sum_e.  sumBuffAdd(pst->ddata, mdStatisticNonBondedPotential);
    checkCUDAError ("SPMERecIk::applyInteraction sum energy");
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
  cal_Q
      <<<meshGridDim, meshBlockDim>>> (
	  K,
	  vecAStar,
	  order,
	  nlist_n,
	  nlist_list,
	  nlist_stride,
	  Q);
  checkCUDAError ("SPMERecIk::calQ cal_Q");
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

