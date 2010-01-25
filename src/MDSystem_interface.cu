#define CPLUSPLUS
#include "MDSystem_interface.h"
#include "GromacsFileManager.h"
#include "FileManager.h"
#include "MDException.h"

MDSystem::MDSystem()
{
  xdfile = NULL;
  xdx = NULL;
  tmpNAtomType = 0;
  hasBond = false;
  hasAngle = false;
  // setNULL (&hdata);
  // setNULL (&ddata);
}


void MDSystem::initConfig (const char * configfile, const char * mapfile,
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
  
  tmpNAtomType = readAtomNameMapFile (mapfile, hdata.numAtom, hdata.atomName,
				      hdata.type, hdata.mass, hdata.charge) ;
  initMass (&hdata);

  printf ("# total %d atoms found, %d types are presented in mapping file\n",
	  hdata.numAtom, tmpNAtomType);

  for (IndexType i = 0; i < hdata.numAtom; ++i){
    hdata.forcx[i] = 0.f;
    hdata.forcy[i] = 0.f;
    hdata.forcz[i] = 0.f;
  }
  
  hdata.NFreedom = hdata.numAtom * 3;

  fclose (fpc);
  
  ////////////////////////////////////////////////////////////
  // init device system
  ////////////////////////////////////////////////////////////

  initDeviceMDData (&hdata, &ddata);
  initDeviceMDData (&hdata, &recoveredDdata);
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


ScalorType MDSystem::
calMaxNBRcut()
{
  return nbInter.maxRcut();
}

void MDSystem::
addNonBondedInteraction(const TypeType &i,
				const TypeType &j,
			const NonBondedInteractionParameter & param)
{
  nbInter.add (i, j, param);
}

void MDSystem::
buildNonBondedInteraction ()
{
  nbInter.build();
}



MDSystem::~MDSystem()
{
  freeAPointer ((void **)&xdx);
}

  
void MDSystem::initBond ()
{
  hasBond = true;
  bdlist.init (ddata);
}

void MDSystem::addBond (const IndexType & ii,
			const IndexType & jj,
			const BondInteractionParameter & param)
{
  bdlist.addBond(ii, jj, param);
}

void MDSystem::buildBond ()
{
  bdlist.build();
}

void MDSystem::initAngle ()
{
  hasAngle = true;
  anglelist.init (ddata);
}

void MDSystem::addAngle (const IndexType & ii,
			 const IndexType & jj,
			 const IndexType & kk,
			 const AngleInteractionParameter & param)
{
  anglelist.addAngle(ii, jj, kk, param);
}

void MDSystem::buildAngle ()
{
  anglelist.build();
}

void MDSystem::updateHost(MDTimer *timer)
{
  if (timer != NULL) timer->tic(mdTimeDataTransfer);
  cpyDeviceMDDataToHost (&ddata, &hdata);
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
