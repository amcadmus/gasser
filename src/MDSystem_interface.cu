#define DEVICE_CODE
#define CPLUSPLUS
#include "MDSystem_interface.h"
#include "GromacsFileManager.h"
#include "FileManager.h"
#include "MDException.h"
#include "Reshuffle_interface.h"

MDSystem::MDSystem()
{
  xdfile = NULL;
  xdx = NULL;
  tmpNAtomType = 0;
  // setNULL (&hdata);
  // setNULL (&ddata);
}

void MDSystem::
normalizeDeviceData (MDTimer * timer)
{
  if (timer != NULL) timer->tic(mdTimeNormalizeSys);
  
  IndexType nob = (ddata.numAtom + DefaultNThreadPerBlock - 1) / DefaultNThreadPerBlock;
  dim3 atomGridDim = toGridDim (nob);

  normalizeSystem
      <<<atomGridDim, DefaultNThreadPerBlock>>> (
	  box,
	  ddata.numAtom,
	  ddata.coord,
	  ddata.coordNoix,
	  ddata.coordNoiy,
	  ddata.coordNoiz);
  checkCUDAError ("NeighborList::rebuild, normalize System");
  if (timer != NULL) timer->toc(mdTimeNormalizeSys);
}


void MDSystem::initConfig (const char * configfile,
			   // const char * mapfile,
			   const IndexType & maxNumAtom)
{
  FILE * fpc = fopen (configfile, "r");
  if (fpc == NULL) {
    throw MDExcptCannotOpenFile ("MDSystem::initConfig:", configfile);
  }
  while (fgetc(fpc) != '\n');

  IndexType numAtom, numMem;
  if (fscanf (fpc, "%d", &(numAtom)) != 1){
    throw MDExcptWrongFileFormat ("MDSystem::initConfig", configfile);
  }
  if (maxNumAtom != 0) {
    numMem = maxNumAtom;
  }
  else {
    numMem = numAtom;
  }
  mallocHostMDData (numAtom, numMem, &hdata);

  IndexType * tmpatomIndex = (IndexType * )malloc(sizeof(IndexType) * numMem);
  if (tmpatomIndex == NULL){
    throw MDExcptFailedMallocOnHost ("MDSystem::initConfig", "tmpatomIndex",
				     sizeof(IndexType) * numMem);
  }
  ScalorType bx, by, bz;
#ifdef COORD_IN_ONE_VEC
  ScalorType * tmpx, * tmpy, * tmpz;
  tmpx = (ScalorType *)malloc (sizeof(ScalorType) * numMem);
  tmpy = (ScalorType *)malloc (sizeof(ScalorType) * numMem);
  tmpz = (ScalorType *)malloc (sizeof(ScalorType) * numMem);
#endif
  GromacsFileManager::readGroFile (configfile,
				   hdata.resdIndex, hdata.resdName, 
				   hdata.atomName, hdata.atomIndex,
#ifndef COORD_IN_ONE_VEC
				   hdata.coordx, hdata.coordy, hdata.coordz,
#else
				   tmpx, tmpy, tmpz,
#endif
				   hdata.velox,  hdata.veloy,  hdata.veloz,
				   &bx, &by, &bz) ;
#ifdef COORD_IN_ONE_VEC
  for (IndexType i = 0; i < numAtom; ++i){
    hdata.coord[i].x = tmpx[i];
    hdata.coord[i].y = tmpy[i];
    hdata.coord[i].z = tmpz[i];
  }
  free (tmpx);
  free (tmpy);
  free (tmpz);
#endif
  freeAPointer ((void**)&tmpatomIndex);
  RectangularBoxGeometry::setBoxSize (bx, by, bz, &box);
  
  // tmpNAtomType = readAtomNameMapFile (mapfile, hdata.numAtom, hdata.atomName,
  // 				      hdata.type, hdata.mass, hdata.charge) ;
  // initMass (&hdata);

  printf ("# total %d atoms found, %d types are presented in mapping file\n",
	  hdata.numAtom, tmpNAtomType);

  for (IndexType i = 0; i < hdata.numAtom; ++i){
    hdata.forcx[i] = 0.f;
    hdata.forcy[i] = 0.f;
    hdata.forcz[i] = 0.f;
  }
  
  hdata.NFreedom = hdata.numAtom * 3;

  fclose (fpc);  
}


void MDSystem::
initTopology (const Topology::System & sysTop)
{
  if (hdata.numAtom != sysTop.indexShift.back()){
    throw MDExcptWrongNumberAtomDataTopology ();
  }
  unsigned shift = 0;
  for (unsigned i = 0; i < sysTop.molecules.size(); ++i){
    for (unsigned j = 0; j < sysTop.numbers[i]; ++j){
      for (unsigned k = 0; k < sysTop.molecules[i].size(); ++k){
	hdata.mass[shift] = sysTop.molecules[i].atoms[k].mass;
	hdata.charge[shift] = sysTop.molecules[i].atoms[k].charge;
	hdata.type[shift] = sysTop.molecules[i].atoms[k].type;
	shift++;
      }
    }
  }
  initMass (&hdata);
}	       

static __global__ void
init_backMapTable (const IndexType numAtom,
		   IndexType * backMapTable)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  if (ii < numAtom){
    backMapTable[ii] = ii;
  }
}

void MDSystem::
initDeviceData ()
{ 
  ////////////////////////////////////////////////////////////
  // init device system
  ////////////////////////////////////////////////////////////
  initDeviceMDData (&hdata, &ddata);
  initDeviceMDData (&hdata, &recoveredDdata);
  initDeviceMDData (&hdata, &bkDdata);
  cudaMalloc ((void**)&backMapTable, sizeof(IndexType) * hdata.numAtom);
  cudaMalloc ((void**)&backMapTableBuff, sizeof(IndexType) * hdata.numAtom);
  checkCUDAError ("MDSystem::initDeviceMDData, malloc back map table");
  
  dim3 myBlockDim, atomGridDim;
  myBlockDim.x = DefaultNThreadPerBlock;
  IndexType nob;
  if (hdata.numAtom % myBlockDim.x == 0){
    nob = hdata.numAtom / myBlockDim.x;
  } else {
    nob = hdata.numAtom / myBlockDim.x + 1;
  }
  atomGridDim = toGridDim (nob);
  init_backMapTable
      <<<atomGridDim, myBlockDim>>> (
	  bkDdata.numAtom, backMapTable);
}


void MDSystem::writeHostDataGro (const char * filename,
				 int step,
				 float time,
				 MDTimer *timer)
{
  if (timer != NULL) timer->tic(mdTimeDataIO);
  FILE * fp = fopen (filename, "w");
  if (fp == NULL){
    throw MDExcptCannotOpenFile (filename);
  }
  fprintf (fp, "# at time = %f, step = %d", time, step);
#ifdef COORD_IN_ONE_VEC
  ScalorType * tmpx, * tmpy, * tmpz;
  tmpx = (ScalorType *)malloc (sizeof(ScalorType) * hdata.numMem);
  tmpy = (ScalorType *)malloc (sizeof(ScalorType) * hdata.numMem);
  tmpz = (ScalorType *)malloc (sizeof(ScalorType) * hdata.numMem);
  for (IndexType i = 0; i < hdata.numAtom; ++i){
    tmpx[i] = hdata.coord[i].x;
    tmpy[i] = hdata.coord[i].y;
    tmpz[i] = hdata.coord[i].z;
  }
#endif
  GromacsFileManager::writeGroFile (fp,
				    hdata.numAtom,
				    hdata.resdIndex, hdata.resdName, 
				    hdata.atomName, hdata.atomIndex,
#ifndef COORD_IN_ONE_VEC
				    hdata.coordx, hdata.coordy, hdata.coordz,
#else
				    tmpx, tmpy, tmpz,
#endif
				    hdata.velox,  hdata.veloy,  hdata.veloz,
				    box.size.x, box.size.y, box.size.z) ;
#ifdef COORD_IN_ONE_VEC
  free (tmpx);
  free (tmpy);
  free (tmpz);
#endif
  fclose (fp);
  if (timer != NULL) timer->toc(mdTimeDataIO);
}


MDSystem::~MDSystem()
{
  freeAPointer ((void **)&xdx);
}

  
void MDSystem::updateHostFromDevice (MDTimer *timer)
{
  if (timer != NULL) timer->tic(mdTimeDataTransfer);
  cpyDeviceMDDataToHost (&ddata, &hdata);
  if (timer != NULL) timer->toc(mdTimeDataTransfer);
}

void MDSystem::updateHostFromRecovered (MDTimer *timer)
{
  if (timer != NULL) timer->tic(mdTimeDataTransfer);
  cpyDeviceMDDataToHost (&recoveredDdata, &hdata);
  if (timer != NULL) timer->toc(mdTimeDataTransfer);
}

void MDSystem::initWriteXtc (const char * filename, float prec)
{
  xdfile = NULL;
  xdfile = xdrfile_open (filename, "w");
  if (xdfile == NULL){
    MDExcptCannotOpenFile ("MDSystem::initWriteXtc", filename);
  }
  for (unsigned i = 0; i < 3; ++i){
    for (unsigned j = 0; j < 3; ++j){
      xdbox[i][j] = 0.f;
    }	      
  }
  xdbox[0][0] = box.size.x;
  xdbox[1][1] = box.size.y;
  xdbox[2][2] = box.size.z;
  xdx = (rvec *) malloc (sizeof(rvec) * hdata.numMem);
  if (xdx == NULL){
    MDExcptFailedMallocOnHost ("MDSystem::initWriteXtc", "xdx", sizeof(rvec) * hdata.numMem);
  }
  xdprec = prec;
}

void MDSystem::writeHostDataXtc (int step, float time, MDTimer *timer)
{
  if (timer != NULL) timer->tic(mdTimeDataIO);
  for (IndexType i = 0; i < hdata.numAtom; ++i){
#ifndef COORD_IN_ONE_VEC
    xdx[i][0] = hdata.coordx[i];
    xdx[i][1] = hdata.coordy[i];
    xdx[i][2] = hdata.coordz[i];
#else
    xdx[i][0] = hdata.coord[i].x;
    xdx[i][1] = hdata.coord[i].y;
    xdx[i][2] = hdata.coord[i].z;
#endif
  }
  xdbox[0][0] = box.size.x;
  xdbox[1][1] = box.size.y;
  xdbox[2][2] = box.size.z;
  write_xtc (xdfile, hdata.numAtom, step, time, xdbox, xdx, xdprec);
  if (timer != NULL) timer->tic(mdTimeDataIO);
}

void MDSystem::endWriteXtc()
{
  xdrfile_close(xdfile);
}


void MDSystem::
recoverDeviceData (MDTimer * timer)
{
  if (timer != NULL) timer->tic(mdTimeDataTransfer);

  dim3 myBlockDim, atomGridDim;
  myBlockDim.x = DefaultNThreadPerBlock;
  IndexType nob;
  if (hdata.numAtom % myBlockDim.x == 0){
    nob = hdata.numAtom / myBlockDim.x;
  } else {
    nob = hdata.numAtom / myBlockDim.x + 1;
  }
  atomGridDim = toGridDim (nob);

  Reshuffle_reshuffleArray
      <<<atomGridDim, myBlockDim>>>
      (ddata.coord, ddata.numAtom, backMapTable, recoveredDdata.coord);
  Reshuffle_reshuffleArray
      <<<atomGridDim, myBlockDim>>>
      (ddata.coordNoix, ddata.numAtom, backMapTable, recoveredDdata.coordNoix);
  Reshuffle_reshuffleArray
      <<<atomGridDim, myBlockDim>>>
      (ddata.coordNoiy, ddata.numAtom, backMapTable, recoveredDdata.coordNoiy);
  Reshuffle_reshuffleArray
      <<<atomGridDim, myBlockDim>>>
      (ddata.coordNoiz, ddata.numAtom, backMapTable, recoveredDdata.coordNoiz);
  Reshuffle_reshuffleArray
      <<<atomGridDim, myBlockDim>>>
      (ddata.velox, ddata.numAtom, backMapTable, recoveredDdata.velox);
  Reshuffle_reshuffleArray
      <<<atomGridDim, myBlockDim>>>
      (ddata.veloy, ddata.numAtom, backMapTable, recoveredDdata.veloy);
  Reshuffle_reshuffleArray
      <<<atomGridDim, myBlockDim>>>
      (ddata.veloz, ddata.numAtom, backMapTable, recoveredDdata.veloz);
  Reshuffle_reshuffleArray
      <<<atomGridDim, myBlockDim>>>
      (ddata.forcx, ddata.numAtom, backMapTable, recoveredDdata.forcx);
  Reshuffle_reshuffleArray
      <<<atomGridDim, myBlockDim>>>
      (ddata.forcy, ddata.numAtom, backMapTable, recoveredDdata.forcy);
  Reshuffle_reshuffleArray
      <<<atomGridDim, myBlockDim>>>
      (ddata.forcz, ddata.numAtom, backMapTable, recoveredDdata.forcz);
  Reshuffle_reshuffleArray
      <<<atomGridDim, myBlockDim>>>
      (ddata.type, ddata.numAtom, backMapTable, recoveredDdata.type);
  Reshuffle_reshuffleArray
      <<<atomGridDim, myBlockDim>>>
      (ddata.mass, ddata.numAtom, backMapTable, recoveredDdata.mass);
  Reshuffle_reshuffleArray
      <<<atomGridDim, myBlockDim>>>
      (ddata.massi, ddata.numAtom, backMapTable, recoveredDdata.massi);
  Reshuffle_reshuffleArray
      <<<atomGridDim, myBlockDim>>>
      (ddata.charge, ddata.numAtom, backMapTable, recoveredDdata.charge);  
  if (timer != NULL) timer->tic(mdTimeDataTransfer);
}


static __global__ void
Reshuffle_calBackMapTable (const IndexType numAtom,
			   const IndexType * backMapTableBuff,
			   const IndexType * idxTable,
			   IndexType *backMapTable)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  if (ii < numAtom){
    backMapTable[idxTable[ii]] = backMapTableBuff[ii];
  }
}

void MDSystem::
reshuffle (const IndexType * indexTable,
	   const IndexType & numAtom,
	   MDTimer * timer)
{
  if (timer != NULL) timer->tic(mdTimeReshuffleSystem);

  cpyDeviceMDDataToDevice (&ddata, &bkDdata);
  
  IndexType nob;
  dim3 myBlockDim;
  myBlockDim.x = DefaultNThreadPerBlock;
  if (ddata.numAtom % myBlockDim.x == 0){
    nob = ddata.numAtom / myBlockDim.x;
  } else {
    nob = ddata.numAtom / myBlockDim.x + 1;
  }
  dim3 atomGridDim = toGridDim (nob);
  Reshuffle_reshuffleArray
      <<<atomGridDim, myBlockDim>>>
      (bkDdata.coord, ddata.numAtom, indexTable, ddata.coord);
  Reshuffle_reshuffleArray
      <<<atomGridDim, myBlockDim>>>
      (bkDdata.coordNoix, ddata.numAtom, indexTable, ddata.coordNoix);
  Reshuffle_reshuffleArray
      <<<atomGridDim, myBlockDim>>>
      (bkDdata.coordNoiy, ddata.numAtom, indexTable, ddata.coordNoiy);
  Reshuffle_reshuffleArray
      <<<atomGridDim, myBlockDim>>>
      (bkDdata.coordNoiz, ddata.numAtom, indexTable, ddata.coordNoiz);
  Reshuffle_reshuffleArray
      <<<atomGridDim, myBlockDim>>>
      (bkDdata.velox, ddata.numAtom, indexTable, ddata.velox);
  Reshuffle_reshuffleArray
      <<<atomGridDim, myBlockDim>>>
      (bkDdata.veloy, ddata.numAtom, indexTable, ddata.veloy);
  Reshuffle_reshuffleArray
      <<<atomGridDim, myBlockDim>>>
      (bkDdata.veloz, ddata.numAtom, indexTable, ddata.veloz);
  Reshuffle_reshuffleArray
      <<<atomGridDim, myBlockDim>>>
      (bkDdata.forcx, ddata.numAtom, indexTable, ddata.forcx);
  Reshuffle_reshuffleArray
      <<<atomGridDim, myBlockDim>>>
      (bkDdata.forcy, ddata.numAtom, indexTable, ddata.forcy);
  Reshuffle_reshuffleArray
      <<<atomGridDim, myBlockDim>>>
      (bkDdata.forcz, ddata.numAtom, indexTable, ddata.forcz);
  Reshuffle_reshuffleArray
      <<<atomGridDim, myBlockDim>>>
      (bkDdata.type, ddata.numAtom, indexTable, ddata.type);
  Reshuffle_reshuffleArray
      <<<atomGridDim, myBlockDim>>>
      (bkDdata.mass, ddata.numAtom, indexTable, ddata.mass);
  Reshuffle_reshuffleArray
      <<<atomGridDim, myBlockDim>>>
      (bkDdata.massi, ddata.numAtom, indexTable, ddata.massi);
  Reshuffle_reshuffleArray
      <<<atomGridDim, myBlockDim>>>
      (bkDdata.charge, ddata.numAtom, indexTable, ddata.charge);  

  cudaMemcpy (backMapTableBuff, backMapTable, sizeof(IndexType) * ddata.numAtom,
	      cudaMemcpyDeviceToDevice);
  Reshuffle_calBackMapTable
      <<<atomGridDim, myBlockDim>>> (
	  ddata.numAtom,
	  backMapTableBuff,
	  indexTable,
	  backMapTable);

  if (timer != NULL) timer->toc(mdTimeReshuffleSystem);
}

