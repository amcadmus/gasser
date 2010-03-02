#include "Parallel_MDData.h"
#include "GromacsFileManager.h"

Parallel::HostMDData::
HostMDData()
    : numAtom_(0), memSize_(0),
      totalMass (0.f), totalMassi (0.f), totalNumFreedom (0),
      coord (NULL),
      coordNoix (NULL), coordNoiy(NULL), coordNoiz(NULL),
      velox (NULL), veloy(NULL), veloz(NULL),
      forcx (NULL), forcy(NULL), forcz(NULL),
      globalIndex (NULL),
      type (NULL), mass (NULL), charge (NULL),
      atomName (NULL), atomIndex (NULL),
      resdName (NULL), resdIndex (NULL)
{
}

Parallel::HostMDData::
~HostMDData ()
{
  clearAll ();
}

void Parallel::HostMDData::
clearAll ()
{
  numAtom_ = 0;
  memSize_ = 0;
  totalNumFreedom = 0;
  totalMass = 0;
  totalMassi = 0;
  
  freeAPointer ((void**)&coord);
  freeAPointer ((void**)&coordNoix);
  freeAPointer ((void**)&coordNoiy);
  freeAPointer ((void**)&coordNoiz);
  freeAPointer ((void**)&velox);
  freeAPointer ((void**)&veloy);
  freeAPointer ((void**)&veloz);
  freeAPointer ((void**)&forcx);
  freeAPointer ((void**)&forcy);
  freeAPointer ((void**)&forcz);
  freeAPointer ((void**)&globalIndex);
  freeAPointer ((void**)&type);
  freeAPointer ((void**)&mass);
  freeAPointer ((void**)&charge);
  freeAPointer ((void**)&atomName);
  freeAPointer ((void**)&atomIndex);
  freeAPointer ((void**)&resdName);
  freeAPointer ((void**)&resdIndex);
}


void Parallel::HostMDData::
reallocAll (const IndexType & memSize__)
{
  reallocCoordinate	(memSize__);
  reallocVelocity	(memSize__);
  reallocForce		(memSize__);
  reallocGlobalIndex	(memSize__);
  reallocTopProperty	(memSize__);
  reallocGroProperty	(memSize__);
}


void Parallel::HostMDData::
reallocCoordinate (const IndexType & memSize__)
{
  if (memSize_ != 0 && memSize_ != memSize__){
    throw MDExcptInconsistentMemorySizeOnHostMDData ();
  }
  if (memSize__ == 0) return ;
  memSize_ = memSize__;

  size_t sizecoord =memSize_ * sizeof(CoordType);
  size_t sizei = memSize_ * sizeof(IntScalorType);

  coord = (CoordType *) realloc (coord, sizecoord);
  if (coord == NULL) throw (MDExcptFailedReallocOnHost("coord", sizecoord));  
  coordNoix = (IntScalorType *) realloc (coordNoix, sizei);
  if (coordNoix == NULL) throw MDExcptFailedReallocOnHost ("coordNoix", sizei);
  coordNoiy = (IntScalorType *) realloc (coordNoiy, sizei);
  if (coordNoiy == NULL) throw MDExcptFailedReallocOnHost ("coordNoiy", sizei);
  coordNoiz = (IntScalorType *) realloc (coordNoiz, sizei);
  if (coordNoiz == NULL) throw MDExcptFailedReallocOnHost ("coordNoiz", sizei);
}

void Parallel::HostMDData::
reallocVelocity (const IndexType & memSize__)
{
  if (memSize_ != 0 && memSize_ != memSize__){
    throw MDExcptInconsistentMemorySizeOnHostMDData ();
  }
  if (memSize__ == 0) return ;
  memSize_ = memSize__;
  
  size_t sizef = memSize_ * sizeof(ScalorType);

  velox = (ScalorType *) realloc (velox, sizef);
  if (velox == NULL) throw MDExcptFailedReallocOnHost("velox", sizef);
  veloy = (ScalorType *) realloc (veloy, sizef);
  if (veloy == NULL) throw MDExcptFailedReallocOnHost("veloy", sizef);
  veloz = (ScalorType *) realloc (veloz, sizef);
  if (veloz == NULL) throw MDExcptFailedReallocOnHost("veloz", sizef);
}

void Parallel::HostMDData::
reallocForce (const IndexType & memSize__)
{
  if (memSize_ != 0 && memSize_ != memSize__){
    throw MDExcptInconsistentMemorySizeOnHostMDData ();
  }
  if (memSize__ == 0) return ;
  memSize_ = memSize__;
  
  size_t sizef = memSize_ * sizeof(ScalorType);
  
  forcx = (ScalorType *) realloc (forcx, sizef);
  if (forcx == NULL) throw MDExcptFailedReallocOnHost("forcx", sizef);
  forcy = (ScalorType *) realloc (forcy, sizef);
  if (forcy == NULL) throw MDExcptFailedReallocOnHost("forcy", sizef);
  forcz = (ScalorType *) realloc (forcz, sizef);
  if (forcz == NULL) throw MDExcptFailedReallocOnHost("forcz", sizef);
}

void Parallel::HostMDData::
reallocGlobalIndex (const IndexType & memSize__)
{
  if (memSize_ != 0 && memSize_ != memSize__){
    throw MDExcptInconsistentMemorySizeOnHostMDData ();
  }
  if (memSize__ == 0) return ;
  memSize_ = memSize__;

  size_t sizeidx = memSize_ * sizeof(IndexType);

  globalIndex = (IndexType *) realloc (globalIndex, sizeidx);
  if (globalIndex == NULL) throw MDExcptFailedReallocOnHost("globalIndex", sizeidx);
}

void Parallel::HostMDData::
reallocTopProperty (const IndexType & memSize__)
{
  if (memSize_ != 0 && memSize_ != memSize__){
    throw MDExcptInconsistentMemorySizeOnHostMDData ();
  }
  if (memSize__ == 0) return ;
  memSize_ = memSize__;

  size_t sizef = memSize_ * sizeof(ScalorType);

  type = (TypeType *) realloc (type, memSize_ * sizeof(TypeType));
  if (type == NULL) throw MDExcptFailedReallocOnHost("type", memSize_ * sizeof(TypeType));
  mass = (ScalorType *) realloc (mass, sizef);
  if (mass == NULL) throw MDExcptFailedReallocOnHost("mass", sizef);
  charge = (ScalorType *) realloc(charge, sizef);
  if (charge == NULL) throw MDExcptFailedReallocOnHost("charge", sizef);
}

void Parallel::HostMDData::
reallocGroProperty (const IndexType & memSize__)
{
  if (memSize_ != 0 && memSize_ != memSize__){
    throw MDExcptInconsistentMemorySizeOnHostMDData ();
  }
  if (memSize__ == 0) return ;
  memSize_ = memSize__;

  size_t sizec = memSize_ * sizeof(char) * StringSize;
  size_t sizeidx = memSize_ * sizeof(IndexType);
  
  atomName = (char *) realloc (atomName, sizec);
  if (atomName == NULL) throw MDExcptFailedReallocOnHost("atomName", sizec);
  atomIndex = (IndexType *) realloc (atomIndex, sizeidx);
  if (atomIndex == NULL) throw MDExcptFailedReallocOnHost("atomIndex", sizeidx);
  resdName = (char *) realloc (resdName, sizec);
  if (resdName == NULL) throw MDExcptFailedReallocOnHost("resdName", sizec);
  resdIndex = (IndexType *) realloc (resdIndex, sizeidx); 
  if (resdIndex == NULL) throw MDExcptFailedReallocOnHost("resdIndex", sizeidx);
}

IndexType Parallel::HostMDData::
numAtomInGroFile (const char * filename)
{
  FILE * fpc = fopen (filename, "r");
  if (fpc == NULL) {
    throw MDExcptCannotOpenFile ("HostMDData::numAtomInGroFile", filename);
  }
  while (fgetc(fpc) != '\n');

  IndexType numAtom;
  if (fscanf (fpc, "%d", &(numAtom)) != 1){
    throw MDExcptWrongFileFormat ("HostMDData::numAtomInGroFile", filename);
  }

  return numAtom;
}


void Parallel::HostMDData::
initCoordinateVelocity_GroFile (const char * filename)
{
  FILE * fpc = fopen (filename, "r");
  if (fpc == NULL) {
    throw MDExcptCannotOpenFile ("HostMDData::initCoordinateVelocity_GroFile", filename);
  }
  while (fgetc(fpc) != '\n');

  if (fscanf (fpc, "%d", &(numAtom_)) != 1){
    throw MDExcptWrongFileFormat ("HostMDData::initCoordinateVelocity_GroFile", filename);
  }

  if (numAtom_ > memSize_) {
    throw MDExcptNumAtomMoreThanMemSize ();
  }

  ScalorType bx, by, bz;
  ScalorType * tmpx, * tmpy, * tmpz;
  tmpx = (ScalorType *)malloc (sizeof(ScalorType) * memSize_);
  if (tmpx == NULL){
    throw MDExcptFailedMallocOnHost ("HostMDData::initCoordinateVelocity_GroFile",
				     "tmpx",
				     sizeof(ScalorType) * memSize_);
  }
  tmpy = (ScalorType *)malloc (sizeof(ScalorType) * memSize_);
  if (tmpy == NULL){
    throw MDExcptFailedMallocOnHost ("HostMDData::initCoordinateVelocity_GroFile",
				     "tmpy",
				     sizeof(ScalorType) * memSize_);
  }
  tmpz = (ScalorType *)malloc (sizeof(ScalorType) * memSize_);
  if (tmpz == NULL){
    throw MDExcptFailedMallocOnHost ("HostMDData::initCoordinateVelocity_GroFile",
				     "tmpz",
				     sizeof(ScalorType) * memSize_);
  }
  GromacsFileManager::readGroFile (filename,
				   resdIndex, resdName, 
				   atomName, atomIndex,
				   tmpx, tmpy, tmpz,
				   velox,  veloy,  veloz,
				   &bx, &by, &bz) ;
  for (IndexType i = 0; i < numAtom_; ++i){
    coord[i].x = tmpx[i];
    coord[i].y = tmpy[i];
    coord[i].z = tmpz[i];
  }
  free (tmpx);
  free (tmpy);
  free (tmpz);

  RectangularBoxGeometry::setBoxSize (bx, by, bz, &box);

  fclose (fpc);  
}


void Parallel::HostMDData::
initMass ()
{
  totalMass = 0;
  for (IndexType i = 0; i < numAtom_; ++i){
    totalMass += mass[i];
  }
  totalMassi = 1./totalMass;  
}


IndexType Parallel::HostMDData::
findMolIndex (const Topology::System & sysTop,
	      const IndexType & globalIndex)
{
  const std::vector<IndexType > & indexShift (sysTop.indexShift);
  if (indexShift.size() < 2){
    throw MDExcptInvalidTopology ();
  }

  IndexType min = 0;
  IndexType max = indexShift.size() - 1;
  while (max - min > 1){
    IndexType mid = (max + min) >> 1;
    if (globalIndex >= indexShift[mid]){
      min = mid;
    }
    else {
      max = mid;
    }
  }
  return min;
}


void Parallel::HostMDData::
initTopology (const Topology::System & sysTop)
{
  for (IndexType i = 0; i < numAtom_; ++i){
    IndexType molIndex = findMolIndex (sysTop, globalIndex[i]);
    IndexType atomIndex = (globalIndex[i] - sysTop.indexShift[molIndex]) %
	(sysTop.molecules[molIndex].size());
    mass[i] = sysTop.molecules[molIndex].atoms[atomIndex].mass;
    charge[i] = sysTop.molecules[molIndex].atoms[atomIndex].charge;
    type[i] = sysTop.molecules[molIndex].atoms[atomIndex].type;
  }

  initMass ();

  totalNumFreedom = numAtom_ * 3;
}


