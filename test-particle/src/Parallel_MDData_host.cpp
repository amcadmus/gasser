#define HOST_CODE

#include "Parallel_MDData.h"
#include "Parallel_TransferEngine.h"
#include "Parallel_Interface.h"
#include "GromacsFileManager.h"

#include "compile_error_mixcode.h"

using namespace Parallel::HostAllocator;

Parallel::HostMDData::
HostMDData(const HostMallocType_t & mallocType_)
    : mallocType (mallocType_),
      _numData(0), _memSize(0),
      coord (NULL),
      coordNoi (NULL),
      velox (NULL), veloy(NULL), veloz(NULL),
      forcx (NULL), forcy(NULL), forcz(NULL),
      globalIndex (NULL),
      type (NULL), mass (NULL), charge (NULL),
      numBond(NULL),
      bondIndex(NULL),
      bondNeighbor_globalIndex(NULL),
      maxNumBond(0),
      numAngle(NULL),
      angleIndex(NULL),
      anglePosi(NULL),
      angleNeighbor_globalIndex(NULL),
      maxNumAngle(0),
      numDihedral(NULL),
      dihedralIndex(NULL),
      dihedralPosi(NULL),
      dihedralNeighbor_globalIndex(NULL),
      maxNumDihedral(0)
{
}

Parallel::HostMDData::
~HostMDData ()
{
  clear ();
}

void Parallel::HostMDData::
clear ()
{
  _numData = 0;
  _memSize = 0;
  
  hostFree ((void**)&coord, mallocType);
  hostFree ((void**)&coordNoi, mallocType);
  hostFree ((void**)&velox, mallocType);
  hostFree ((void**)&veloy, mallocType);
  hostFree ((void**)&veloz, mallocType);
  hostFree ((void**)&forcx, mallocType);
  hostFree ((void**)&forcy, mallocType);
  hostFree ((void**)&forcz, mallocType);
  hostFree ((void**)&globalIndex, mallocType);
  hostFree ((void**)&type, mallocType);
  hostFree ((void**)&mass, mallocType);
  hostFree ((void**)&charge, mallocType);

  hostFree ((void**)&bondNeighbor_globalIndex, mallocType);
  hostFree ((void**)&bondIndex, mallocType);
  hostFree ((void**)&numBond, mallocType);
  maxNumBond = 0;
  hostFree ((void**)&angleNeighbor_globalIndex, mallocType);
  hostFree ((void**)&angleIndex, mallocType);
  hostFree ((void**)&anglePosi, mallocType);
  hostFree ((void**)&numAngle, mallocType);
  maxNumAngle = 0;
  hostFree ((void**)&dihedralNeighbor_globalIndex, mallocType);
  hostFree ((void**)&dihedralIndex, mallocType);
  hostFree ((void**)&dihedralPosi, mallocType);
  hostFree ((void**)&numDihedral, mallocType);
  maxNumDihedral = 0;  
}



void Parallel::HostMDData::
fillZero ()
{
  HostCoordType coord0;
  HostCoordNoiType coordNoi0 ;
  
  for (IndexType i = 0; i < memSize(); ++i){
    coord[i] = coord0;
    coordNoi[i] = coordNoi0;
    velox[i] = 0.f;
    veloy[i] = 0.f;
    veloz[i] = 0.f;
    forcx[i] = 0.f;
    forcy[i] = 0.f;
    forcz[i] = 0.f;
    globalIndex[i] = MaxIndexValue;
    type[i] = 0;
    mass[i] = 0.f;
    charge[i] = 0.f;
  }
  if (maxNumBond != 0){
    for (IndexType i = 0; i < memSize(); ++i){
      numBond[i] = 0;
    }
    for (IndexType j = 0; j < maxNumBond; ++j){
      for (IndexType i = 0; i < memSize(); ++i){
	IndexType index = j * bondTopStride() + i;
	bondNeighbor_globalIndex[index] = MaxIndexValue;
	bondIndex[index] = 0;
      }
    }
  }
  if (maxNumAngle != 0){
    for (IndexType i = 0; i < memSize(); ++i){
      numAngle[i] = 0;
    }
    for (IndexType j = 0; j < maxNumAngle; ++j){
      for (IndexType i = 0; i < memSize(); ++i){
	IndexType index = j * bondTopStride() + i;
	angleIndex[index] = 0;
	anglePosi[index] = 0;
      }
    }
    for (IndexType j = 0; j < 2 * maxNumAngle; ++j){
      for (IndexType i = 0; i < memSize(); ++i){
	IndexType index = j * bondTopStride() + i;
    	angleNeighbor_globalIndex[index] = MaxIndexValue;
      }
    }
  }
  if (maxNumDihedral != 0){
    for (IndexType i = 0; i < memSize(); ++i){
      numDihedral[i] = 0;
    }
    for (IndexType j = 0; j < maxNumDihedral; ++j){
      for (IndexType i = 0; i < memSize(); ++i){
	IndexType index = j * bondTopStride() + i;
	dihedralIndex[index] = 0;
	dihedralPosi[index] = 0;
      }
    }
    for (IndexType j = 0; j < 3 * maxNumDihedral; ++j){
      for (IndexType i = 0; i < memSize(); ++i){
	IndexType index = j * bondTopStride() + i;
    	dihedralNeighbor_globalIndex[index] = MaxIndexValue;
      }
    }
  }
}

void Parallel::HostMDData::
easyMalloc (const IndexType memSize_,
	    const IndexType maxNumBond_,
	    const IndexType maxNumAngle_,
	    const IndexType maxNumDihedral_)
{
  // printf ("# malloc HostMDData\n");
  clear ();

  _memSize = memSize_;
  maxNumBond = maxNumBond_;
  maxNumAngle = maxNumAngle_;
  maxNumDihedral = maxNumDihedral_;

  if (_memSize == 0) return;

  size_t sizecoord = memSize() * sizeof(HostCoordType);
  size_t sizecoordNoi = memSize() * sizeof(HostCoordNoiType);
  size_t sizef = memSize() * sizeof(ScalorType);
  size_t sizeidx = memSize() * sizeof(IndexType);
  size_t sizetype = memSize() * sizeof (TypeType);

  hostMalloc ((void**)&coord, sizecoord, mallocType);
  hostMalloc ((void**)&coordNoi, sizecoordNoi, mallocType);
  hostMalloc ((void**)&velox, sizef, mallocType);
  hostMalloc ((void**)&veloy, sizef, mallocType);
  hostMalloc ((void**)&veloz, sizef, mallocType);
  hostMalloc ((void**)&forcx, sizef, mallocType);
  hostMalloc ((void**)&forcy, sizef, mallocType);
  hostMalloc ((void**)&forcz, sizef, mallocType);
  hostMalloc ((void**)&globalIndex, sizeidx, mallocType);
  hostMalloc ((void**)&type, sizetype, mallocType);
  hostMalloc ((void**)&mass, sizef, mallocType);
  hostMalloc ((void**)&charge, sizef, mallocType);
  
  // coord = (HostCoordType *) malloc (sizecoord);
  // if (coord == NULL) throw MDExcptFailedMallocOnHost ("coord", sizecoord);
  // coordNoi = (HostCoordNoiType *) malloc (sizecoordNoi);
  // if (coordNoi == NULL) throw MDExcptFailedMallocOnHost ("coordNoi", sizecoordNoi);
  // velox = (ScalorType *) malloc (sizef);
  // if (velox == NULL) throw MDExcptFailedMallocOnHost ("velox", sizef);
  // veloy = (ScalorType *) malloc (sizef);
  // if (veloy == NULL) throw MDExcptFailedMallocOnHost ("veloy", sizef);
  // veloz = (ScalorType *) malloc (sizef);
  // if (veloz == NULL) throw MDExcptFailedMallocOnHost ("veloz", sizef);
  // forcx = (ScalorType *) malloc (sizef);
  // if (forcx == NULL) throw MDExcptFailedMallocOnHost ("forcx", sizef);
  // forcy = (ScalorType *) malloc (sizef);
  // if (forcy == NULL) throw MDExcptFailedMallocOnHost ("forcy", sizef);
  // forcz = (ScalorType *) malloc (sizef);
  // if (forcz == NULL) throw MDExcptFailedMallocOnHost ("forcz", sizef);
  // globalIndex = (IndexType *) malloc (sizeidx);
  // if (globalIndex == NULL) throw MDExcptFailedMallocOnHost ("globalIndex", sizeidx);
  // type = (TypeType *) malloc (sizetype);
  // if (type == NULL) throw MDExcptFailedMallocOnHost ("type", sizetype);
  // mass = (ScalorType *) malloc (sizef);
  // if (mass == NULL) throw MDExcptFailedMallocOnHost ("mass", sizef);
  // charge = (ScalorType *) malloc (sizef);
  // if (charge == NULL) throw MDExcptFailedMallocOnHost ("charge", sizef);  

  size_t size0 = sizeof(IndexType) * memSize();

  if (maxNumBond != 0){
    size_t size1 = size0 * maxNumBond;
    hostMalloc ((void**)&numBond, size0, mallocType);
    hostMalloc ((void**)&bondIndex, size1, mallocType);
    hostMalloc ((void**)&bondNeighbor_globalIndex, size1, mallocType);    
    // numBond = (IndexType *) malloc (size0);
    // if (numBond == NULL) throw MDExcptFailedMallocOnHost ("numBond", size0);
    // bondIndex = (IndexType *) malloc (size1);
    // if (bondIndex == NULL) throw MDExcptFailedMallocOnHost ("bondIndex", size1);
    // bondNeighbor_globalIndex = (IndexType *) malloc (size1);
    // if (bondNeighbor_globalIndex == NULL){
    //   throw MDExcptFailedMallocOnHost ("bondNeighbor_globalIndex", size1);
    // }
  }
  if (maxNumAngle != 0){
    size_t size1 = size0 * maxNumAngle;
    size_t size2 = size0 * maxNumAngle * 2;
    hostMalloc ((void**)&numAngle, size0, mallocType);
    hostMalloc ((void**)&angleIndex, size1, mallocType);
    hostMalloc ((void**)&anglePosi , size1, mallocType);
    hostMalloc ((void**)&angleNeighbor_globalIndex, size2, mallocType);
    
    // numAngle = (IndexType *) malloc (size0);
    // if (numAngle == NULL) throw MDExcptFailedMallocOnHost ("numAngle", size0);
    // angleIndex = (IndexType *) malloc (size1);
    // if (angleIndex == NULL) throw MDExcptFailedMallocOnHost ("angleIndex", size1);
    // anglePosi = (IndexType *) malloc (size1);
    // if (anglePosi == NULL) throw MDExcptFailedMallocOnHost ("anglePosi", size1);
    // angleNeighbor_globalIndex = (IndexType *) malloc (size2);
    // if (angleNeighbor_globalIndex == NULL){
    //   throw MDExcptFailedMallocOnHost ("angleNeighbor_globalIndex", size2);
    // }
  }
  if (maxNumDihedral != 0){
    size_t size1 = size0 * maxNumDihedral;
    size_t size2 = size0 * maxNumDihedral * 3;
    hostMalloc ((void**)&numDihedral, size0, mallocType);
    hostMalloc ((void**)&dihedralIndex, size1, mallocType);
    hostMalloc ((void**)&dihedralPosi , size1, mallocType);
    hostMalloc ((void**)&dihedralNeighbor_globalIndex, size2, mallocType);
    
    // numDihedral = (IndexType *) malloc (size0);
    // if (numDihedral == NULL) throw MDExcptFailedMallocOnHost ("numDihedral", size0);
    // dihedralIndex = (IndexType *) malloc (size1);
    // if (dihedralIndex == NULL) throw MDExcptFailedMallocOnHost ("dihedralIndex", size1);
    // dihedralPosi = (IndexType *) malloc (size1);
    // if (dihedralPosi == NULL) throw MDExcptFailedMallocOnHost ("dihedralPosi", size1);
    // dihedralNeighbor_globalIndex = (IndexType *) malloc (size2);
    // if (dihedralNeighbor_globalIndex == NULL){
    //   throw MDExcptFailedMallocOnHost ("dihedralNeighbor_globalIndex", size2);
    // }
  }

  fillZero ();
}


Parallel::HostMDData::
HostMDData (const HostMDData & hdata,
	    const HostMallocType_t & mallocType_)
    : mallocType (mallocType_),
      _numData(0), _memSize(0),
      coord (NULL),
      coordNoi (NULL),
      velox (NULL), veloy(NULL), veloz(NULL),
      forcx (NULL), forcy(NULL), forcz(NULL),
      globalIndex (NULL),
      type (NULL), mass (NULL), charge (NULL),
      numBond(NULL),
      bondIndex(NULL),
      bondNeighbor_globalIndex(NULL),
      maxNumBond(0),
      numAngle(NULL),
      angleIndex(NULL),
      anglePosi(NULL),
      angleNeighbor_globalIndex(NULL),
      maxNumAngle(0),
      numDihedral(NULL),
      dihedralIndex(NULL),
      dihedralPosi(NULL),
      dihedralNeighbor_globalIndex(NULL),
      maxNumDihedral(0)
{
  this->copy (hdata);
}

void Parallel::HostMDData::
copy (const HostMDData & hdata,
      const MDDataItemMask_t mask)
{
  // if (!mask) return;
  IndexType expectedNumBond(0), expectedNumAngle(0), expectedNumDihedral(0);
  bool copyBond = (mask & MDDataItemMask_Bond);
  bool copyAngle = (mask & MDDataItemMask_Angle);
  bool copyDihedral = (mask & MDDataItemMask_Dihedral);  
  if (copyBond){
    expectedNumBond = hdata.getMaxNumBond();
  }
  if (copyAngle){
    expectedNumAngle = hdata.getMaxNumAngle();
  }
  if (copyDihedral){
    expectedNumDihedral = hdata.getMaxNumDihedral();
  }

  if (memSize() < hdata.numData()){
    easyMalloc (hdata.numData(), expectedNumBond, expectedNumAngle, expectedNumDihedral);
  }
  else if ((copyBond && (maxNumBond != hdata.maxNumBond)) ||
	   (copyAngle && (maxNumAngle != hdata.maxNumAngle)) ||
	   (copyDihedral && (maxNumDihedral != hdata.maxNumDihedral)) ){
    easyMalloc (memSize(), expectedNumBond, expectedNumAngle, expectedNumDihedral);
  }

  _numData = hdata.numData();
  setGlobalBox (hdata.getGlobalBox());  
  
  size_t sizecoord = hdata.numData() * sizeof(HostCoordType);
  size_t sizecoordNoi = hdata.numData() * sizeof(HostCoordNoiType);
  size_t sizef = hdata.numData() * sizeof(ScalorType);
  size_t sizeidx = hdata.numData() * sizeof(IndexType);
  size_t sizetype = hdata.numData() * sizeof (TypeType);

  if (mask & MDDataItemMask_Coordinate) {
    memcpy (coord, hdata.coord, sizecoord);
  }
  if (mask & MDDataItemMask_CoordinateNoi){
    memcpy (coordNoi, hdata.coordNoi, sizecoordNoi);
  }
  if (mask & MDDataItemMask_Velocity){
    memcpy (velox, hdata.velox, sizef);
    memcpy (veloy, hdata.veloy, sizef);
    memcpy (veloz, hdata.veloz, sizef);
  }
  if (mask & MDDataItemMask_Force){
    memcpy (forcx, hdata.forcx, sizef);
    memcpy (forcy, hdata.forcy, sizef);
    memcpy (forcz, hdata.forcz, sizef);
  }
  if (mask & MDDataItemMask_GlobalIndex){
    memcpy (globalIndex, hdata.globalIndex, sizetype);
  }
  if (mask & MDDataItemMask_Type){
    memcpy (type, hdata.type, sizetype);
  }
  if (mask & MDDataItemMask_Mass){
    memcpy (mass, hdata.mass, sizef);
  }
  if (mask & MDDataItemMask_Charge){
    memcpy (charge, hdata.charge, sizef);  
  }

  size_t size0 = sizeof(IndexType) * hdata.numData();

  if (copyBond){
    // maxNumBond = hdata.maxNumBond;
    for (IndexType i = 0; i < maxNumBond; ++i){
      if (i == 0) memcpy (numBond, hdata.numBond, size0);
      memcpy (i * bondTopStride() + bondNeighbor_globalIndex,
	      i * bondTopStride() + hdata.bondNeighbor_globalIndex,
	      size0);
      memcpy (i * bondTopStride() + bondIndex,
	      i * bondTopStride() + hdata.bondIndex,
	      size0);
    }
  }

  if (copyAngle ){
    // maxNumAngle = hdata.maxNumAngle;
    for (IndexType i = 0; i < maxNumAngle; ++i){
      if (i == 0) memcpy (numAngle, hdata.numAngle, size0);
      memcpy (i * bondTopStride() + angleIndex,
	      i * bondTopStride() + hdata.angleIndex,
	      size0);
      memcpy (i * bondTopStride() + anglePosi,
	      i * bondTopStride() + hdata.anglePosi,
	      size0);
    }
    for (IndexType i = 0; i < maxNumAngle * 2; ++i){
      memcpy (i * bondTopStride() + angleNeighbor_globalIndex,
	      i * bondTopStride() + hdata.angleNeighbor_globalIndex,
	      size0);
    }
  }

  if (copyDihedral){
    for (IndexType i = 0; i < maxNumDihedral; ++i){
      if (i == 0) memcpy (numDihedral, hdata.numDihedral, size0);
      memcpy (i * size0 + dihedralIndex,
	      i * size0 + hdata.dihedralIndex,
	      size0);
      memcpy (i * size0 + dihedralPosi,
	      i * size0 + hdata.dihedralPosi,
	      size0);
    }
    for (IndexType i = 0; i < maxNumDihedral * 3; ++i){
      memcpy (i * size0 + dihedralNeighbor_globalIndex,
	      i * size0 + hdata.dihedralNeighbor_globalIndex,
	      size0);
    }
  }

}


void Parallel::HostMDData::
pushBackAtom  (const HostCoordType & coord_,
	       const HostCoordNoiType & coordNoi_,
	       const ScalorType & velox_,
	       const ScalorType & veloy_,
	       const ScalorType & veloz_,
	       const IndexType & globalIndex_,
	       const TypeType & type_,
	       const ScalorType & mass_,
	       const ScalorType & charge_)
{
  if (memSize() < numData() + 1){
    MDDataItemMask_t mask = MDDataItemMask_All;
    mask ^= MDDataItemMask_Bond;
    mask ^= MDDataItemMask_Angle;
    mask ^= MDDataItemMask_Dihedral;
    HostMDData backup;
    backup.copy (*this, mask);
    _memSize ++;
    _memSize <<=1;
    easyMalloc (memSize(), maxNumBond, maxNumAngle, maxNumDihedral);
    copy (backup, mask);
  }
  coord[numData()] = coord_;
  coordNoi[numData()].x = coordNoi_.x;
  coordNoi[numData()].y = coordNoi_.y;
  coordNoi[numData()].z = coordNoi_.z;
  velox[numData()] = velox_;
  veloy[numData()] = veloy_;
  veloz[numData()] = veloz_;
  globalIndex[numData()] = globalIndex_;
  type[numData()] = type_;
  mass[numData()] = mass_;
  charge[numData()] = charge_;
  numData() ++;
}


IndexType Parallel::GlobalHostMDData::
numAtomInGroFile (const char * filename)
{
  FILE * fpc = fopen (filename, "r");
  if (fpc == NULL) {
    throw MDExcptCannotOpenFile ("HostMDData::numAtomInGroFile", filename);
  }
  while (fgetc(fpc) != '\n');

  IndexType numData;
  if (fscanf (fpc, "%d", &(numData)) != 1){
    throw MDExcptWrongFileFormat ("HostMDData::numDataInGroFile", filename);
  }

  return numData;
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

  IndexType expectedNumAtom;
  if (fscanf (fpc, "%d", &(expectedNumAtom)) != 1){
    throw MDExcptWrongFileFormat ("HostMDData::initConf_GroFile", filename);
  }

  if (expectedNumAtom > memSize()) {
    easyMalloc (expectedNumAtom);
  }
  numData() = expectedNumAtom;

  ScalorType bx, by, bz;
  ScalorType * tmpx, * tmpy, * tmpz;
  tmpx = (ScalorType *)malloc (sizeof(ScalorType) * numData());
  if (tmpx == NULL){
    throw MDExcptFailedMallocOnHost ("HostMDData::initConf_GroFile",
				     "tmpx",
				     sizeof(ScalorType) * numData());
  }
  tmpy = (ScalorType *)malloc (sizeof(ScalorType) * numData());
  if (tmpy == NULL){
    throw MDExcptFailedMallocOnHost ("HostMDData::initConf_GroFile",
				     "tmpy",
				     sizeof(ScalorType) * numData());
  }
  tmpz = (ScalorType *)malloc (sizeof(ScalorType) * numData());
  if (tmpz == NULL){
    throw MDExcptFailedMallocOnHost ("HostMDData::initConf_GroFile",
				     "tmpz",
				     sizeof(ScalorType) * numData());
  }
  GromacsFileManager::readGroFile (filename,
				   resdIndex, resdName, 
				   atomName, atomIndex,
				   tmpx, tmpy, tmpz,
				   velox,  veloy,  veloz,
				   &bx, &by, &bz) ;
  for (IndexType i = 0; i < numData(); ++i){
    globalIndex[i] = i;
    coord[i].x = tmpx[i];
    coord[i].y = tmpy[i];
    coord[i].z = tmpz[i];
    coordNoi[i].x = 0;
    coordNoi[i].y = 0;
    coordNoi[i].z = 0;
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
initTopology_property (const Topology::System & sysTop)
{
  for (IndexType i = 0; i < numData(); ++i){
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
  fprintf (fp, "# %d\n", numData());
  for (IndexType i = 0; i < numData(); ++i){
    fprintf (fp, "%8.3f %8.3f %8.3f\n", coord[i].x, coord[i].y, coord[i].z);
  }
  fclose (fp);
}


void Parallel::
distributeGlobalMDData (const GlobalHostMDData & gdata,
			HostMDData & ldata)
{
  Parallel::TransferEngine transSend ;
  Parallel::TransferEngine transRecv ;
  
  int Nx, Ny, Nz;
  Parallel::Interface::numProcDim (Nx, Ny, Nz);
  ScalorType hx(0), hy(0), hz(0);

  HostMDData sdata;
  int myRank = Parallel::Interface::myRank();
  
  if (myRank == 0){
    // IndexType naiveLocalMemSize  = (gdata.numData() /
    // 				    Parallel::Interface::numProc() * 2);
    IndexType naiveLocalMemSize = gdata.numData();
    if (naiveLocalMemSize < 100) naiveLocalMemSize = 100 ;
    sdata.easyMalloc   (naiveLocalMemSize);
    sdata.setGlobalBox (gdata.getGlobalBox());
    
    hx = (gdata.getGlobalBox().size.x) / Nx;
    hy = (gdata.getGlobalBox().size.y) / Ny;
    hz = (gdata.getGlobalBox().size.z) / Nz;
  }
  
  for (int ix = 0; ix < Nx; ++ix){
    for (int iy = 0; iy < Ny; ++iy){
      for (int iz = 0; iz < Nz; ++iz){

	if (myRank == 0){
	  sdata.clearData();
	  ScalorType lx = ix * hx;
	  ScalorType ux = lx + hx;
	  ScalorType ly = iy * hy;
	  ScalorType uy = ly + hy;
	  ScalorType lz = iz * hz;
	  ScalorType uz = lz + hz;
	  for (IndexType i = 0; i < gdata.numData(); ++i){
	    if (gdata.cptr_coordinate()[i].x >= lx &&
		gdata.cptr_coordinate()[i].x <  ux &&
		gdata.cptr_coordinate()[i].y >= ly &&
		gdata.cptr_coordinate()[i].y <  uy &&
		gdata.cptr_coordinate()[i].z >= lz &&
		gdata.cptr_coordinate()[i].z <  uz){
	      sdata.pushBackAtom (
		  gdata.cptr_coordinate()[i],
		  gdata.cptr_coordinateNoi()[i],
		  gdata.cptr_velocityX()[i],
		  gdata.cptr_velocityY()[i],
		  gdata.cptr_velocityZ()[i],
		  gdata.cptr_globalIndex()[i],
		  gdata.cptr_type()[i],
		  gdata.cptr_mass()[i],
		  gdata.cptr_charge()[i]);
	    }
	  }
	}

	int workingRank;
	Parallel::Interface::cartCoordToRank (ix, iy, iz, workingRank);
	int recvNumAtom, recvMemSize;
	int sendNumAtom, sendMemSize;
	
	if (myRank == 0){
	  sendNumAtom = sdata.numData();
	  sendMemSize = sdata.memSize();
	  transSend.clearRegistered ();
	  transSend.registerBuff (&sendNumAtom, sizeof(IndexType));
	  transSend.registerBuff (&sendMemSize, sizeof(IndexType));
	  transSend.build ();
	  transSend.Isend (workingRank, 0);
	  // transSend.wait ();
	}
	if (myRank == workingRank){
	  transRecv.clearRegistered();
	  transRecv.registerBuff (&recvNumAtom, sizeof(IndexType));
	  transRecv.registerBuff (&recvMemSize, sizeof(IndexType));
	  transRecv.build ();
	  transRecv.Irecv (0, 0);
	  transRecv.wait ();
	  ldata.easyMalloc (recvMemSize);
	}
	if (myRank == 0){
	  transSend.wait ();
	}
	if (myRank == 0){
	  transSend.clearRegistered ();
	  transSend.registerBuff (sdata.cptr_coordinate(),
				  sizeof(HostCoordType) * sdata.numData());
	  transSend.registerBuff (sdata.cptr_coordinateNoi(),
				  sizeof(HostCoordNoiType) * sdata.numData());
	  transSend.registerBuff (sdata.cptr_velocityX(),
				  sizeof(ScalorType) * sdata.numData());
	  transSend.registerBuff (sdata.cptr_velocityY(),
				  sizeof(ScalorType) * sdata.numData());
	  transSend.registerBuff (sdata.cptr_velocityZ(),
				  sizeof(ScalorType) * sdata.numData());
	  transSend.registerBuff (sdata.cptr_globalIndex(),
				  sizeof(IndexType) * sdata.numData());
	  transSend.registerBuff (sdata.cptr_type(),
				  sizeof(TypeType) * sdata.numData());
	  transSend.registerBuff (sdata.cptr_mass(),
				  sizeof(ScalorType) * sdata.numData());
	  transSend.registerBuff (sdata.cptr_charge(),
				  sizeof(ScalorType) * sdata.numData());
	  ScalorType tmpx, tmpy, tmpz;
	  tmpx = sdata.getGlobalBox().size.x;
	  tmpy = sdata.getGlobalBox().size.y;
	  tmpz = sdata.getGlobalBox().size.z;
	  transSend.registerBuff (&tmpx,
				  sizeof(ScalorType) * 1);
	  transSend.registerBuff (&tmpy,
				  sizeof(ScalorType) * 1);
	  transSend.registerBuff (&tmpz,
				  sizeof(ScalorType) * 1);
	  transSend.build ();
	  transSend.Isend (workingRank, 0);
	}
	if (myRank == workingRank){
	  transRecv.clearRegistered();
	  transRecv.registerBuff (ldata.cptr_coordinate(),
				  sizeof(HostCoordType) * recvNumAtom);
	  transRecv.registerBuff (ldata.cptr_coordinateNoi(),
				  sizeof(HostCoordNoiType) * recvNumAtom);
	  transRecv.registerBuff (ldata.cptr_velocityX(),
				  sizeof(ScalorType) * recvNumAtom);
	  transRecv.registerBuff (ldata.cptr_velocityY(),
				  sizeof(ScalorType) * recvNumAtom);
	  transRecv.registerBuff (ldata.cptr_velocityZ(),
				  sizeof(ScalorType) * recvNumAtom);
	  transRecv.registerBuff (ldata.cptr_globalIndex(),
				  sizeof(IndexType) * recvNumAtom);
	  transRecv.registerBuff (ldata.cptr_type(),
				  sizeof(TypeType) * recvNumAtom);
	  transRecv.registerBuff (ldata.cptr_mass(),
				  sizeof(ScalorType) * recvNumAtom);
	  transRecv.registerBuff (ldata.cptr_charge(),
				  sizeof(ScalorType) * recvNumAtom);
	  ScalorType tmpx, tmpy, tmpz;
	  transRecv.registerBuff (&tmpx,
				  sizeof(ScalorType) * 1);
	  transRecv.registerBuff (&tmpy,
				  sizeof(ScalorType) * 1);
	  transRecv.registerBuff (&tmpz,
				  sizeof(ScalorType) * 1);
	  transRecv.build ();
	  transRecv.Irecv (0, 0);
	  transRecv.wait ();
	  ldata.setGlobalBox (tmpx, tmpy, tmpz);
	  ldata.numData () = recvNumAtom;
	}
	if (myRank == 0){
	  transSend.wait ();
	}
      }
    }
  }

  (Parallel::Interface::barrier());
}


  
void Parallel::GlobalHostMDData::
writeData_GroFile (const char * filename,
		   const char * atomName, const IndexType * atomIndex,
		   const char * resdName, const IndexType * resdIndex)
{
  FILE * fp = fopen (filename, "w");
  if (fp == NULL){
    throw MDExcptCannotOpenFile (filename);
  }
  // fprintf (fp, "# at time = %f, step = %d", time, step);
  ScalorType * tmpx, * tmpy, * tmpz;
  tmpx = (ScalorType *)malloc (sizeof(ScalorType) * numData());
  tmpy = (ScalorType *)malloc (sizeof(ScalorType) * numData());
  tmpz = (ScalorType *)malloc (sizeof(ScalorType) * numData());
  for (IndexType i = 0; i < numData(); ++i){
    tmpx[i] = coord[i].x;
    tmpy[i] = coord[i].y;
    tmpz[i] = coord[i].z;
  }
  GromacsFileManager::writeGroFile (fp,
				    numData(),
				    resdIndex, resdName, 
				    atomName, atomIndex,
				    tmpx,  tmpy,  tmpz,
				    velox, veloy, veloz,
				    getGlobalBox().size.x,
				    getGlobalBox().size.y,
				    getGlobalBox().size.z) ;
  free (tmpx);
  free (tmpy);
  free (tmpz);
  fclose (fp);
}

void Parallel::HostMDData::
mallocFromHost (const HostMDData & hdata)
{
  setGlobalBox (hdata.getGlobalBox());
  easyMalloc (hdata.memSize(), hdata.getMaxNumBond(), hdata.getMaxNumAngle(),
	      hdata.getMaxNumDihedral());
  _numData = 0;
}

