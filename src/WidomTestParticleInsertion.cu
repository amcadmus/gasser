#include "WidomTestParticleInsertion.h"
#include "RandomGenerator.h"

WidomTestParticleInsertion::
WidomTestParticleInsertion ()
    : inited (false), nParticleInserted (0)
{
}

WidomTestParticleInsertion::
~WidomTestParticleInsertion ()
{
  clear ();
}

void WidomTestParticleInsertion::
reinit (const ScalorType & temperature_,
	const IndexType & nParticleInserted_,
	const TypeType & particleType)
{
  clear ();
  
  mytemperature = temperature_;
  nParticleInserted = nParticleInserted_;

  size_t sizec = sizeof(CoordType) * nParticleInserted;
  size_t sizet = sizeof(TypeType)  * nParticleInserted;
//  size_t sizef = sizeof(ScalorType)* nParticleInserted;
  hcoord = (HostCoordType *) malloc (sizec);
  if (hcoord == NULL) throw MDExcptFailedMallocOnHost ("WidomTestParticleInsertion::hcoord", sizec);
  htype  = (TypeType *) malloc (sizet);
  if (htype == NULL) throw MDExcptFailedMallocOnHost ("WidomTestParticleInsertion::htype", sizet);
  cudaMalloc ((void **) &coordTestParticle, sizec);
  checkCUDAError ("WidomTestParticleInsertion::reinit coordTestParticle");
  cudaMalloc ((void **) &typeTestParticle,  sizet);
  checkCUDAError ("WidomTestParticleInsertion::reinit typeTestParticle");
  cudaMalloc ((void **) &dresult, sizeof(ScalorType));
  
  for (unsigned i = 0; i < nParticleInserted; ++i){
    hcoord[i].x = hcoord[i].y = hcoord[i].z = hcoord[i].w = 0.;
    htype[i] = particleType;
  }
  cudaMemcpy (coordTestParticle, hcoord, sizec, cudaMemcpyHostToDevice);
  checkCUDAError ("WidomTestParticleInsertion::reinit copy from hcoord to coordTestParticle");
  cudaMemcpy (typeTestParticle, htype, sizet, cudaMemcpyHostToDevice);
  checkCUDAError ("WidomTestParticleInsertion::reinit copy from htype to typeTestParticle");

  sumExpDeltaU.reinit (numTestParticle(), NThreadForSum);
  
  inited = true;
}

void WidomTestParticleInsertion::
clear ()
{
  nParticleInserted = 0;
  if (inited){
    freeAPointer ((void**)&hcoord);
    freeAPointer ((void**)&htype);
    cudaFree (coordTestParticle);
    checkCUDAError ("WidomTestParticleInsertion::clear free coordTestParticle");
    cudaFree (typeTestParticle);
    checkCUDAError ("WidomTestParticleInsertion::clear free typeTestParticle");
  }
  inited = false;
}

ScalorType WidomTestParticleInsertion::
expMu ()
{
  sumExpDeltaU.sumBuff(dresult, 0);
  cudaMemcpy (&hresult, dresult, sizeof(ScalorType), cudaMemcpyDeviceToHost);
  hresult /= numTestParticle();
  return hresult;
}

void WidomTestParticleInsertion::
generateTestCoords (const MDSystem & sys)
{
  for (unsigned i = 0; i < numTestParticle(); ++i){
    hcoord[i].x = sys.box.size.x * RandomGenerator_MT19937::genrand_real2();
    hcoord[i].y = sys.box.size.y * RandomGenerator_MT19937::genrand_real2();
    hcoord[i].z = sys.box.size.z * RandomGenerator_MT19937::genrand_real2();
  }
  size_t sizec = sizeof(CoordType) * nParticleInserted;
  cudaMemcpy (coordTestParticle, hcoord, sizec, cudaMemcpyHostToDevice);
  checkCUDAError ("WidomTestParticleInsertion::generateTestCoords copy from hcoord to coordTestParticle");
}




