#include "EwaldSumRec.h"

void EwaldSumRec::
applyInteraction (MDSystem & sys,
		  MDStatistic * st,
		  MDTimer * timer)
{
  cal_Sm
      <<<meshGridDim, meshBlockDim,
      sizeof(CoordType) * meshBlockDim.x>>> (
	  K,
	  vecAStar,
	  sys.ddata.coord,
	  sys.ddata.charge,
	  sys.ddata.numAtom,
	  Sm,
	  nele);
  checkCUDAError ("EwaldSumRec::applyInteraction cal_sm");
  applyForce
      <<<atomGridDim, atomBlockDim,
      (sizeof(ScalorType) + sizeof(cufftComplex)) * atomBlockDim.x + 4>>> (
	  K,
	  vecAStar,
	  fm,
	  Sm,
	  nele,
	  sys.ddata.coord,
	  sys.ddata.charge,
	  sys.ddata.forcx,
	  sys.ddata.forcy,
	  sys.ddata.forcz,
	  sys.ddata.numAtom);
  checkCUDAError ("EwaldSumRec::applyInteraction applyForce");
}



void EwaldSumRec::
calV()
{
  volume =
      vecA.xx * (vecA.yy*vecA.zz - vecA.zy*vecA.yz) - 
      vecA.xy * (vecA.yx*vecA.zz - vecA.zx*vecA.yz) +
      vecA.xz * (vecA.yx*vecA.zy - vecA.zx*vecA.yy);
}
  
void EwaldSumRec::
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

void EwaldSumRec::
reinit (const MatrixType & vecA_,
	const IntVectorType & K_,
	const ScalorType & beta_,
	const IndexType & natom,
	const IndexType & meshNThread,
	const IndexType & atomNThread)
{
  vecA = vecA_;
  K = K_;
  beta = beta_;
  calV();
  calAStar();
  
  N.x = (K.x/2) * 2 + 1;
  N.y = (K.y/2) * 2 + 1;
  N.z = (K.z/2) * 2 + 1;
  nele = N.x * N.y * N.z;
    
  IndexType nob;
  meshBlockDim.x = meshNThread;
  nob = (nele  + meshBlockDim.x - 1) / meshBlockDim.x;
  meshGridDim = toGridDim (nob);
  atomBlockDim.x = atomNThread;
  nob = (natom + atomBlockDim.x - 1) / atomBlockDim.x;
  atomGridDim = toGridDim (nob);

  cudaMalloc ((void**)&Sm, sizeof(cufftComplex) * nele);
  cudaMalloc ((void**)&fm, sizeof(ScalorType) * nele);
  checkCUDAError ("EwaldSumRec::reinit malloc");
  malloced = true;
  
  cal_fm
      <<<meshGridDim, meshBlockDim>>> (
	  K,
	  vecAStar,
	  beta,
	  volume,
	  fm,
	  nele);
  checkCUDAError ("EwaldSumRec::reinit cal_fm");
}


void EwaldSumRec::
freeAll ()
{
  if (malloced){
    cudaFree (Sm);
    cudaFree (fm);
    malloced = false;
  }
}

EwaldSumRec::
EwaldSumRec()
    : malloced (false)
{
}

EwaldSumRec::
~EwaldSumRec ()
{
  freeAll();
}


      


