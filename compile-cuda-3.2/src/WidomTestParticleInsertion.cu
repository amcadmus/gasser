#include "WidomTestParticleInsertion.h"
#include "RandomGenerator.h"

__global__ static void
processDeltaEnergy_NVT (const IndexType numTestParticle,
			ScalorType * energy_buff,
			const ScalorType energyCorr,
			const IndexType  natom,
			const ScalorType volume,
			const ScalorType temperature)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;

  if (ii < numTestParticle) {
    energy_buff[ii] = expf (-(energy_buff[ii] + energyCorr) / temperature);
    // energy_buff[ii] = expf (-(energy_buff[ii]) / temperature);
  }
}

__global__ static void
processDeltaEnergy_NPT (const IndexType numTestParticle,
			ScalorType * energy_buff,
			const ScalorType energyCorr,
			const IndexType  natom,
			const ScalorType pressure,
			const ScalorType volume,
			const ScalorType temperature)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;

  if (ii < numTestParticle) {
    // energy_buff[ii] = expf ( - (energy_buff[ii] + energyCorr) / temperature)
    // 	 * pressure * volume /
    // 	(temperature * (natom + 1.));
    energy_buff[ii] = expf ( - (energy_buff[ii] + energyCorr) / temperature)
	 * pressure * volume /
	(temperature * (natom + 1.));
  }
}

__global__ static void
processDeltaEnergy_NPT2 (const IndexType numTestParticle,
			 ScalorType * energy_buff,
			 const ScalorType energyCorr,
			 const IndexType  natom,
			 const ScalorType pressure,
			 const ScalorType volume,
			 const ScalorType temperature)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;

  if (ii < numTestParticle) {
    // energy_buff[ii] = expf ( - (energy_buff[ii] + energyCorr) / temperature)
    // 	 * pressure * volume /
    // 	(temperature * (natom + 1.));
    energy_buff[ii] = expf ( - (energy_buff[ii] + energyCorr) / temperature)
	* pressure /
	(temperature * (natom + 1.));
  }
}


WidomTestParticleInsertion_NVT::
WidomTestParticleInsertion_NVT ()
    : inited (false), nParticleInserted (0)
{
}

WidomTestParticleInsertion_NVT::
~WidomTestParticleInsertion_NVT ()
{
  clear ();
}

void WidomTestParticleInsertion_NVT::
reinit (const ScalorType & temperature_,
	const IndexType & nParticleInserted_,
	const TypeType & particleType,
	const SystemNonBondedInteraction & sysNbInter)
{
  clear ();
  
  mytemperature = temperature_;
  nParticleInserted = nParticleInserted_;

  size_t sizec = sizeof(CoordType) * nParticleInserted;
  size_t sizet = sizeof(TypeType)  * nParticleInserted;
//  size_t sizef = sizeof(ScalorType)* nParticleInserted;
  hcoord = (HostCoordType *) malloc (sizec);
  if (hcoord == NULL) throw MDExcptFailedMallocOnHost ("WidomTestParticleInsertion_NVT::hcoord", sizec);
  htype  = (TypeType *) malloc (sizet);
  if (htype == NULL) throw MDExcptFailedMallocOnHost ("WidomTestParticleInsertion_NVT::htype", sizet);
  cudaMalloc ((void **) &coordTestParticle, sizec);
  checkCUDAError ("WidomTestParticleInsertion_NVT::reinit coordTestParticle");
  cudaMalloc ((void **) &typeTestParticle,  sizet);
  checkCUDAError ("WidomTestParticleInsertion_NVT::reinit typeTestParticle");
  cudaMalloc ((void **) &dresult, sizeof(ScalorType));
  
  for (unsigned i = 0; i < nParticleInserted; ++i){
    hcoord[i].x = hcoord[i].y = hcoord[i].z = hcoord[i].w = 0.;
    htype[i] = particleType;
  }
  cudaMemcpy (coordTestParticle, hcoord, sizec, cudaMemcpyHostToDevice);
  checkCUDAError ("WidomTestParticleInsertion_NVT::reinit copy from hcoord to coordTestParticle");
  cudaMemcpy (typeTestParticle, htype, sizet, cudaMemcpyHostToDevice);
  checkCUDAError ("WidomTestParticleInsertion_NVT::reinit copy from htype to typeTestParticle");

  sumExpDeltaU.reinit (numTestParticle(), NThreadForSum);

  energyCorr = sysNbInter.energyCorrection (particleType) * 2.;
  printf ("# energy correction to widom test particle insertion is %f\n",
	  energyCorr);
  
  inited = true;
}

void WidomTestParticleInsertion_NVT::
clear ()
{
  nParticleInserted = 0;
  if (inited){
    freeAPointer ((void**)&hcoord);
    freeAPointer ((void**)&htype);
    cudaFree (coordTestParticle);
    checkCUDAError ("WidomTestParticleInsertion_NVT::clear free coordTestParticle");
    cudaFree (typeTestParticle);
    checkCUDAError ("WidomTestParticleInsertion_NVT::clear free typeTestParticle");
  }
  inited = false;
}

ScalorType WidomTestParticleInsertion_NVT::
expMu ()
{
  processDeltaEnergy_NVT
      <<<numTestParticle() / DefaultNThreadPerBlock + 1, DefaultNThreadPerBlock>>>
      (numTestParticle(),
       sumExpDeltaU.buff,
       energyCorrection(),
       numAtom,
       volume,
       temperature());
  sumExpDeltaU.sumBuff(dresult, 0);
  cudaMemcpy (&hresult, dresult, sizeof(ScalorType), cudaMemcpyDeviceToHost);
  // printf ("# ntest is %d\n", numTestParticle());
  hresult /= numTestParticle();
  return hresult;
}

void WidomTestParticleInsertion_NVT::
generateTestCoords (const MDSystem & sys)
{
  numAtom = sys.ddata.numAtom;
  volume = (sys.box.size.x * sys.box.size.y * sys.box.size.z);
  for (unsigned i = 0; i < numTestParticle(); ++i){
    hcoord[i].x = sys.box.size.x * RandomGenerator_MT19937::genrand_real2();
    hcoord[i].y = sys.box.size.y * RandomGenerator_MT19937::genrand_real2();
    hcoord[i].z = sys.box.size.z * RandomGenerator_MT19937::genrand_real2();
    // hcoord[i].x = 1;
    // hcoord[i].y = 1.8;
    // hcoord[i].z = 1.45;    
  }
  size_t sizec = sizeof(CoordType) * nParticleInserted;
  cudaMemcpy (coordTestParticle, hcoord, sizec, cudaMemcpyHostToDevice);
  checkCUDAError ("WidomTestParticleInsertion_NVT::generateTestCoords copy from hcoord to coordTestParticle");
  scaledEnergyCorr = energyCorr / volume;
}

////////////////////////////////////////////////////////////////////////////////
// npt system
////////////////////////////////////////////////////////////////////////////////


WidomTestParticleInsertion_NPT::
WidomTestParticleInsertion_NPT ()
    : inited (false), nParticleInserted (0)
{
}

WidomTestParticleInsertion_NPT::
~WidomTestParticleInsertion_NPT ()
{
  clear ();
}

void WidomTestParticleInsertion_NPT::
reinit (const ScalorType & temperature_,
	const ScalorType & pressure_,
	const IndexType & nParticleInserted_,
	const TypeType & particleType,
	const SystemNonBondedInteraction & sysNbInter)
{
  clear ();
  
  mytemperature = temperature_;
  mypressure = pressure_;
  nParticleInserted = nParticleInserted_;

  size_t sizec = sizeof(CoordType) * nParticleInserted;
  size_t sizet = sizeof(TypeType)  * nParticleInserted;
//  size_t sizef = sizeof(ScalorType)* nParticleInserted;
  hcoord = (HostCoordType *) malloc (sizec);
  if (hcoord == NULL) throw MDExcptFailedMallocOnHost ("WidomTestParticleInsertion_NPT::hcoord", sizec);
  htype  = (TypeType *) malloc (sizet);
  if (htype == NULL) throw MDExcptFailedMallocOnHost ("WidomTestParticleInsertion_NPT::htype", sizet);
  cudaMalloc ((void **) &coordTestParticle, sizec);
  checkCUDAError ("WidomTestParticleInsertion_NPT::reinit coordTestParticle");
  cudaMalloc ((void **) &typeTestParticle,  sizet);
  checkCUDAError ("WidomTestParticleInsertion_NPT::reinit typeTestParticle");
  cudaMalloc ((void **) &dresult, sizeof(ScalorType));
  
  for (unsigned i = 0; i < nParticleInserted; ++i){
    hcoord[i].x = hcoord[i].y = hcoord[i].z = hcoord[i].w = 0.;
    htype[i] = particleType;
  }
  cudaMemcpy (coordTestParticle, hcoord, sizec, cudaMemcpyHostToDevice);
  checkCUDAError ("WidomTestParticleInsertion_NPT::reinit copy from hcoord to coordTestParticle");
  cudaMemcpy (typeTestParticle, htype, sizet, cudaMemcpyHostToDevice);
  checkCUDAError ("WidomTestParticleInsertion_NPT::reinit copy from htype to typeTestParticle");

  sumExpDeltaU.reinit (numTestParticle(), NThreadForSum);

  energyCorr = sysNbInter.energyCorrection (particleType) * 2.;
  printf ("# energy correction to widom test particle insertion is %f\n",
	  energyCorr);
  
  inited = true;
}

void WidomTestParticleInsertion_NPT::
clear ()
{
  nParticleInserted = 0;
  if (inited){
    freeAPointer ((void**)&hcoord);
    freeAPointer ((void**)&htype);
    cudaFree (coordTestParticle);
    checkCUDAError ("WidomTestParticleInsertion_NPT::clear free coordTestParticle");
    cudaFree (typeTestParticle);
    checkCUDAError ("WidomTestParticleInsertion_NPT::clear free typeTestParticle");
  }
  inited = false;
}

ScalorType WidomTestParticleInsertion_NPT::
expMu ()
{
  processDeltaEnergy_NPT
      <<<numTestParticle() / DefaultNThreadPerBlock + 1, DefaultNThreadPerBlock>>>
      (numTestParticle(),
       sumExpDeltaU.buff,
       energyCorrection(),
       numAtom,
       pressure(),
       volume,
       temperature());
  sumExpDeltaU.sumBuff(dresult, 0);
  cudaMemcpy (&hresult, dresult, sizeof(ScalorType), cudaMemcpyDeviceToHost);
  // printf ("# ntest is %d\n", numTestParticle());
  hresult /= numTestParticle();
  return hresult;
}

void WidomTestParticleInsertion_NPT::
generateTestCoords (const MDSystem & sys)
{
  numAtom = sys.ddata.numAtom;
  volume = (sys.box.size.x * sys.box.size.y * sys.box.size.z);
  for (unsigned i = 0; i < numTestParticle(); ++i){
    hcoord[i].x = sys.box.size.x * RandomGenerator_MT19937::genrand_real2();
    hcoord[i].y = sys.box.size.y * RandomGenerator_MT19937::genrand_real2();
    hcoord[i].z = sys.box.size.z * RandomGenerator_MT19937::genrand_real2();
  }
  size_t sizec = sizeof(CoordType) * nParticleInserted;
  cudaMemcpy (coordTestParticle, hcoord, sizec, cudaMemcpyHostToDevice);
  checkCUDAError ("WidomTestParticleInsertion_NPT::generateTestCoords copy from hcoord to coordTestParticle");
  scaledEnergyCorr = energyCorr / volume;
}



////////////////////////////////////////////////////////////////////////////////
// nvt 2
////////////////////////////////////////////////////////////////////////////////


WidomTestParticleInsertion_NVT2::
WidomTestParticleInsertion_NVT2 ()
    : inited (false), nParticleInserted (0)
{
}

WidomTestParticleInsertion_NVT2::
~WidomTestParticleInsertion_NVT2 ()
{
  clear ();
}

void WidomTestParticleInsertion_NVT2::
reinit (const ScalorType & temperature_,
	const ScalorType & gridSize,
	const RectangularBox & box,
	const TypeType & particleType,
	const SystemNonBondedInteraction & sysNbInter)
{
  clear ();
  
  mytemperature = temperature_;

  HostCoordType		inteCellSize;
  IndexType		inteCellNumX;
  IndexType		inteCellNumY;
  IndexType		inteCellNumZ;
  inteCellNumX = unsigned((box.size.x - 1e-3) / gridSize) + 1;
  inteCellNumY = unsigned((box.size.y - 1e-3) / gridSize) + 1;
  inteCellNumZ = unsigned((box.size.z - 1e-3) / gridSize) + 1;
  inteCellSize.x = box.size.x / ScalorType(inteCellNumX);
  inteCellSize.y = box.size.y / ScalorType(inteCellNumY);
  inteCellSize.z = box.size.z / ScalorType(inteCellNumZ);
  inteCellVolume = inteCellSize.x * inteCellSize.y * inteCellSize.z;

  nParticleInserted = inteCellNumX * inteCellNumY * inteCellNumZ;
  
  size_t sizec = sizeof(CoordType) * nParticleInserted;
  size_t sizet = sizeof(TypeType)  * nParticleInserted;
//  size_t sizef = sizeof(ScalorType)* nParticleInserted;
  hcoord = (HostCoordType *) malloc (sizec);
  if (hcoord == NULL) throw MDExcptFailedMallocOnHost ("WidomTestParticleInsertion_NVT2::hcoord", sizec);
  htype  = (TypeType *) malloc (sizet);
  if (htype == NULL) throw MDExcptFailedMallocOnHost ("WidomTestParticleInsertion_NVT2::htype", sizet);
  cudaMalloc ((void **) &coordTestParticle, sizec);
  checkCUDAError ("WidomTestParticleInsertion_NVT2::reinit coordTestParticle");
  cudaMalloc ((void **) &typeTestParticle,  sizet);
  checkCUDAError ("WidomTestParticleInsertion_NVT2::reinit typeTestParticle");
  cudaMalloc ((void **) &dresult, sizeof(ScalorType));

  for (unsigned i = 0; i < inteCellNumX; ++i){
    for (unsigned j = 0; j < inteCellNumY; ++j){
      for (unsigned k = 0; k < inteCellNumZ; ++k){
	IndexType idx = k + inteCellNumZ * (j + inteCellNumY * i);
	hcoord[idx].x = inteCellSize.x * i + 0.5 * inteCellSize.x;
	hcoord[idx].y = inteCellSize.y * j + 0.5 * inteCellSize.y;
	hcoord[idx].z = inteCellSize.z * k + 0.5 * inteCellSize.z;
	hcoord[idx].w = 0.;
      }
    }
  }  
  for (unsigned i = 0; i < numTestParticle(); ++i){
    htype[i] = particleType;
  }
  
  cudaMemcpy (coordTestParticle, hcoord, sizec, cudaMemcpyHostToDevice);
  checkCUDAError ("WidomTestParticleInsertion_NVT2::reinit copy from hcoord to coordTestParticle");
  cudaMemcpy (typeTestParticle, htype, sizet, cudaMemcpyHostToDevice);
  checkCUDAError ("WidomTestParticleInsertion_NVT2::reinit copy from htype to typeTestParticle");

  sumExpDeltaU.reinit (numTestParticle(), NThreadForSum);

  volume = (box.size.x * box.size.y * box.size.z);
  scaledEnergyCorr = sysNbInter.energyCorrection (particleType) * 2.;
  scaledEnergyCorr /= volume;
  printf ("# energy correction to widom test particle insertion is %f\n",
	  scaledEnergyCorr);
  
  inited = true;
}

void WidomTestParticleInsertion_NVT2::
clear ()
{
  nParticleInserted = 0;
  if (inited){
    freeAPointer ((void**)&hcoord);
    freeAPointer ((void**)&htype);
    cudaFree (coordTestParticle);
    checkCUDAError ("WidomTestParticleInsertion_NVT2::clear free coordTestParticle");
    cudaFree (typeTestParticle);
    checkCUDAError ("WidomTestParticleInsertion_NVT2::clear free typeTestParticle");
  }
  inited = false;
}

ScalorType WidomTestParticleInsertion_NVT2::
expMu ()
{
  processDeltaEnergy_NVT
      <<<numTestParticle() / DefaultNThreadPerBlock + 1, DefaultNThreadPerBlock>>>
      (numTestParticle(),
       sumExpDeltaU.buff,
       energyCorrection(),
       numAtom,
       volume,
       temperature());
  sumExpDeltaU.sumBuff(dresult, 0);
  cudaMemcpy (&hresult, dresult, sizeof(ScalorType), cudaMemcpyDeviceToHost);
  hresult *= inteCellVolume;
  return hresult;
}

void WidomTestParticleInsertion_NVT2::
generateTestCoords (const MDSystem & sys)
{
  numAtom = sys.ddata.numAtom;
  // for (unsigned i = 0; i < numTestParticle(); ++i){
  //   hcoord[i].x = sys.box.size.x * RandomGenerator_MT19937::genrand_real2();
  //   hcoord[i].y = sys.box.size.y * RandomGenerator_MT19937::genrand_real2();
  //   hcoord[i].z = sys.box.size.z * RandomGenerator_MT19937::genrand_real2();
  // }
  // size_t sizec = sizeof(CoordType) * nParticleInserted;
  // cudaMemcpy (coordTestParticle, hcoord, sizec, cudaMemcpyHostToDevice);
  // checkCUDAError ("WidomTestParticleInsertion_NVT2::generateTestCoords copy from hcoord to coordTestParticle");
}



////////////////////////////////////////////////////////////////////////////////
// for npt system version 2
////////////////////////////////////////////////////////////////////////////////


WidomTestParticleInsertion_NPT2::
WidomTestParticleInsertion_NPT2 ()
    : inited (false), nParticleInserted (0)
{
}

WidomTestParticleInsertion_NPT2::
~WidomTestParticleInsertion_NPT2 ()
{
  clear ();
}

void WidomTestParticleInsertion_NPT2::
reinit (const ScalorType & temperature_,
	const ScalorType & pressure_,
	const ScalorType & gridSize_,
	const TypeType & particleType_,
	const SystemNonBondedInteraction & sysNbInter)
{
  clear ();
  
  mytemperature = temperature_;
  mypressure = pressure_;
  gridSize = gridSize_;
  particleType = particleType_;
  
  energyCorr = sysNbInter.energyCorrection (particleType) * 2.;
  printf ("# energy correction to widom test particle insertion is %f\n",
	  energyCorr);
  
  inited = false;
}

void WidomTestParticleInsertion_NPT2::
clear ()
{
  nParticleInserted = 0;
  if (inited){
    freeAPointer ((void**)&hcoord);
    freeAPointer ((void**)&htype);
    cudaFree (coordTestParticle);
    checkCUDAError ("WidomTestParticleInsertion_NPT2::clear free coordTestParticle");
    cudaFree (typeTestParticle);
    checkCUDAError ("WidomTestParticleInsertion_NPT2::clear free typeTestParticle");
  }
  inited = false;
}

ScalorType WidomTestParticleInsertion_NPT2::
expMu ()
{
  processDeltaEnergy_NPT2
      <<<numTestParticle() / DefaultNThreadPerBlock + 1, DefaultNThreadPerBlock>>>
      (numTestParticle(),
       sumExpDeltaU.buff,
       energyCorrection(),
       numAtom,
       pressure(),
       volume,
       temperature());
  sumExpDeltaU.sumBuff(dresult, 0);
  cudaMemcpy (&hresult, dresult, sizeof(ScalorType), cudaMemcpyDeviceToHost);
  hresult *= inteCellVolume;
  return hresult;
}

void WidomTestParticleInsertion_NPT2::
generateTestCoords (const MDSystem & sys)
{
  clear();
  
  const RectangularBox & box (sys.box);
  IndexType inteCellNumX = unsigned((box.size.x - 1e-3) / gridSize) + 1;
  IndexType inteCellNumY = unsigned((box.size.y - 1e-3) / gridSize) + 1;
  IndexType inteCellNumZ = unsigned((box.size.z - 1e-3) / gridSize) + 1;
  HostCoordType inteCellSize;
  inteCellSize.x = box.size.x / ScalorType(inteCellNumX);
  inteCellSize.y = box.size.y / ScalorType(inteCellNumY);
  inteCellSize.z = box.size.z / ScalorType(inteCellNumZ);
  inteCellVolume = inteCellSize.x * inteCellSize.y * inteCellSize.z;

  nParticleInserted = inteCellNumX * inteCellNumY * inteCellNumZ;

  size_t sizec = sizeof(CoordType) * nParticleInserted;
  size_t sizet = sizeof(TypeType)  * nParticleInserted;
//  size_t sizef = sizeof(ScalorType)* nParticleInserted;
  hcoord = (HostCoordType *) malloc (sizec);
  if (hcoord == NULL) throw MDExcptFailedMallocOnHost ("WidomTestParticleInsertion_NVT2::hcoord", sizec);
  htype  = (TypeType *) malloc (sizet);
  if (htype == NULL) throw MDExcptFailedMallocOnHost ("WidomTestParticleInsertion_NVT2::htype", sizet);
  cudaMalloc ((void **) &coordTestParticle, sizec);
  checkCUDAError ("WidomTestParticleInsertion_NVT2::reinit coordTestParticle");
  cudaMalloc ((void **) &typeTestParticle,  sizet);
  checkCUDAError ("WidomTestParticleInsertion_NVT2::reinit typeTestParticle");
  cudaMalloc ((void **) &dresult, sizeof(ScalorType));

  
  for (unsigned i = 0; i < inteCellNumX; ++i){
    for (unsigned j = 0; j < inteCellNumY; ++j){
      for (unsigned k = 0; k < inteCellNumZ; ++k){
	IndexType idx = k + inteCellNumZ * (j + inteCellNumY * i);
	hcoord[idx].x = inteCellSize.x * i + 0.5 * inteCellSize.x;
	hcoord[idx].y = inteCellSize.y * j + 0.5 * inteCellSize.y;
	hcoord[idx].z = inteCellSize.z * k + 0.5 * inteCellSize.z;
	hcoord[idx].w = 0.;
      }
    }
  }  
  for (unsigned i = 0; i < numTestParticle(); ++i){
    htype[i] = particleType;
  }
  cudaMemcpy (coordTestParticle, hcoord, sizec, cudaMemcpyHostToDevice);
  checkCUDAError ("WidomTestParticleInsertion_NVT2::reinit copy from hcoord to coordTestParticle");
  cudaMemcpy (typeTestParticle, htype, sizet, cudaMemcpyHostToDevice);
  checkCUDAError ("WidomTestParticleInsertion_NVT2::reinit copy from htype to typeTestParticle");

  sumExpDeltaU.reinit (numTestParticle(), NThreadForSum);

  volume = (box.size.x * box.size.y * box.size.z);
  numAtom = sys.ddata.numAtom;
  
  scaledEnergyCorr = energyCorr / volume;
  
  inited = true;
}
