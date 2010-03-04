#define HOST_CODE
#define CPP_FILE

#include "Parallel_MDData.h"
#include "GromacsFileManager.h"

Parallel::HostMDData::
HostMDData()
    : numAtom_(0), memSize_(0),
      coord (NULL),
      coordNoix (NULL), coordNoiy(NULL), coordNoiz(NULL),
      velox (NULL), veloy(NULL), veloz(NULL),
      forcx (NULL), forcy(NULL), forcz(NULL),
      globalIndex (NULL),
      type (NULL), mass (NULL), charge (NULL)
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
}


void Parallel::HostMDData::
reallocAll (const IndexType & memSize__)
{
  reallocCoordinate	(memSize__);
  reallocVelocity	(memSize__);
  reallocForce		(memSize__);
  reallocGlobalIndex	(memSize__);
  reallocTopProperty	(memSize__);
  // reallocGroProperty	(memSize__);
}


void Parallel::HostMDData::
reallocCoordinate (const IndexType & memSize__)
{
  if (memSize_ != 0 && memSize_ != memSize__){
    throw MDExcptInconsistentMemorySizeOnHostMDData ();
  }
  if (memSize__ == 0) return ;
  memSize_ = memSize__;

  size_t sizecoord =memSize_ * sizeof(HostCoordType);
  size_t sizei = memSize_ * sizeof(IntScalorType);

  coord = (HostCoordType *) realloc (coord, sizecoord);
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


IndexType Parallel::GlobalHostMDData::
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
pushBackAtom  (const HostCoordType & coord_,
	       const IntScalorType & coordNoix_,
	       const IntScalorType & coordNoiy_,
	       const IntScalorType & coordNoiz_,
	       const ScalorType & velox_,
	       const ScalorType & veloy_,
	       const ScalorType & veloz_,
	       const IndexType & globalIndex_,
	       const TypeType & type_,
	       const ScalorType & mass_,
	       const ScalorType & charge_)
{
  if (numAtom_ == memSize_){
    memSize_ ++;
    memSize_ <<= 1;
    reallocAll (memSize_);
  }
  coord[numAtom_] = coord_;
  coordNoix[numAtom_] = coordNoix_;
  coordNoiy[numAtom_] = coordNoiy_;
  coordNoiz[numAtom_] = coordNoiz_;
  velox[numAtom_] = velox_;
  veloy[numAtom_] = veloy_;
  veloz[numAtom_] = veloz_;
  globalIndex[numAtom_] = globalIndex_;
  type[numAtom_] = type_;
  mass[numAtom_] = mass_;
  charge[numAtom_] = charge_;
  numAtom_ ++;
}



void Parallel::GlobalHostMDData::
initConf_GroFile (const char * filename,
		  char * atomName, IndexType * atomIndex,
		  char * resdName, IndexType * resdIndex)
{
  FILE * fpc = fopen (filename, "r");
  if (fpc == NULL) {
    throw MDExcptCannotOpenFile ("HostMDData::initConf_GroFile", filename);
  }
  while (fgetc(fpc) != '\n');

  if (fscanf (fpc, "%d", &(numAtom_)) != 1){
    throw MDExcptWrongFileFormat ("HostMDData::initConf_GroFile", filename);
  }

  if (numAtom_ > memSize_) {
    throw MDExcptNumAtomMoreThanMemSize ();
  }
  

  ScalorType bx, by, bz;
  ScalorType * tmpx, * tmpy, * tmpz;
  tmpx = (ScalorType *)malloc (sizeof(ScalorType) * memSize_);
  if (tmpx == NULL){
    throw MDExcptFailedMallocOnHost ("HostMDData::initConf_GroFile",
				     "tmpx",
				     sizeof(ScalorType) * memSize_);
  }
  tmpy = (ScalorType *)malloc (sizeof(ScalorType) * memSize_);
  if (tmpy == NULL){
    throw MDExcptFailedMallocOnHost ("HostMDData::initConf_GroFile",
				     "tmpy",
				     sizeof(ScalorType) * memSize_);
  }
  tmpz = (ScalorType *)malloc (sizeof(ScalorType) * memSize_);
  if (tmpz == NULL){
    throw MDExcptFailedMallocOnHost ("HostMDData::initConf_GroFile",
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
    globalIndex[i] = i;
    coord[i].x = tmpx[i];
    coord[i].y = tmpy[i];
    coord[i].z = tmpz[i];
    coordNoix[i] = 0;
    coordNoiy[i] = 0;
    coordNoiz[i] = 0;
  }
  free (tmpx);
  free (tmpy);
  free (tmpz);

  RectangularBoxGeometry::setBoxSize (bx, by, bz, &globalBox);

  fclose (fpc);  
}


IndexType Parallel::GlobalHostMDData::
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


void Parallel::GlobalHostMDData::
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
}


void Parallel::HostMDData::
writeData_SimpleFile (const char * filename)
{
  FILE * fp = fopen (filename, "w");
  if (fp == NULL){
    throw MDExcptCannotOpenFile(filename);
  }
  fprintf (fp, "# %d\n", numAtom_);
  for (IndexType i = 0; i < numAtom_; ++i){
    fprintf (fp, "%8.3f %8.3f %8.3f\n", coord[i].x, coord[i].y, coord[i].z);
  }
  fclose (fp);
}

