#ifndef __MDSystem_h_wanghan__
#define __MDSystem_h_wanghan__

#include "common.h"

class HostMDData
{
public:
  IndexType numAtom;
  IndexType numMem;
  IndexType NFreedom;
#ifndef COORD_IN_ONE_VEC
  ScalorType * coordx;
  ScalorType * coordy;
  ScalorType * coordz;
#else
  CoordType * coord;
#endif
  IntScalorType * coordNoix;
  IntScalorType * coordNoiy;
  IntScalorType * coordNoiz;
  ScalorType * velox;
  ScalorType * veloy;
  ScalorType * veloz;
  ScalorType * forcx;
  ScalorType * forcy;
  ScalorType * forcz;
  TypeType   * type;
  ScalorType * mass;
  ScalorType * massi;
  ScalorType totalMass;
  ScalorType totalMassi;
  ScalorType * charge;
  char * atomName;
  IndexType * atomIndex;
  char * resdName;
  IndexType * resdIndex;
public:
  HostMDData() ;
  ~HostMDData() ;
}
    ;



class DeviceMDData
{
public:
  bool malloced ;
  IndexType numAtom;
  IndexType numMem;
  IndexType NFreedom;
#ifndef COORD_IN_ONE_VEC
  ScalorType * coordx;
  ScalorType * coordy;
  ScalorType * coordz;
#else
  CoordType * coord;
#endif
  IntScalorType * coordNoix;
  IntScalorType * coordNoiy;
  IntScalorType * coordNoiz;
  ScalorType * velox;
  ScalorType * veloy;
  ScalorType * veloz;
  ScalorType * forcx;
  ScalorType * forcy;
  ScalorType * forcz;
  TypeType   * type;
  ScalorType * mass;
  ScalorType * massi;
  ScalorType totalMass;
  ScalorType totalMassi;
  ScalorType * charge;
  // texture<TypeType, 1, cudaReadModeElementType> texReftype;
public:
  DeviceMDData () ;
  ~DeviceMDData () ;
}
    ;


__host__ void cpyHostMDDataToDevice (const HostMDData * hdata, DeviceMDData * ddata);
__host__ void cpyDeviceMDDataToHost (const DeviceMDData * ddata, HostMDData * hdata);
#ifdef COORD_IN_ONE_VEC
__global__ void deviceCpyTypeToCoordW (CoordType * coord,
				       const TypeType * type,
				       const IndexType N);
#endif

__host__ void mallocHostMDData (IndexType numAtom, IndexType expectedMaxNumAtom,
				HostMDData * hdata);
__host__ void lazyInitHostMDData (HostMDData * hdata);
__host__ void initMass (HostMDData * hdata);
__host__ void destroyHostMDData (HostMDData * hdata);

__host__ void initDeviceMDData (const HostMDData * hdata, DeviceMDData * ddata);
__host__ void destroyDeviceMDData (DeviceMDData * ddata);
// __host__ void bindTextureOnDeviceMDData (DeviceMDData * data);

__device__ void cpyDeviceMDDataElement (const DeviceMDData * ddata1,
					const IndexType indx1,
					DeviceMDData * ddata2,
					const IndexType indx2);


////////////////////////////////////////////////////////////
// implementation
// ////////////////////////////////////////////////////////

__device__ void cpyDeviceMDDataElement (const DeviceMDData * ddata1,
					const IndexType indx1,
					DeviceMDData * ddata2,
					const IndexType indx2)
{
#ifndef COORD_IN_ONE_VEC
  ddata2->coordx[indx2] = ddata1->coordx[indx1];
  ddata2->coordy[indx2] = ddata1->coordy[indx1];
  ddata2->coordz[indx2] = ddata1->coordz[indx1];
#else
  (ddata2->coord[indx2]) = (ddata1->coord[indx1]);
#endif
  
  ddata2->coordNoix[indx2] = ddata1->coordNoix[indx1];
  ddata2->coordNoiy[indx2] = ddata1->coordNoiy[indx1];
  ddata2->coordNoiz[indx2] = ddata1->coordNoiz[indx1];

  ddata2->velox[indx2] = ddata1->velox[indx1];
  ddata2->veloy[indx2] = ddata1->veloy[indx1];
  ddata2->veloz[indx2] = ddata1->veloz[indx1];

  ddata2->forcx[indx2] = ddata1->forcx[indx1];
  ddata2->forcy[indx2] = ddata1->forcy[indx1];
  ddata2->forcz[indx2] = ddata1->forcz[indx1];

  ddata2->type[indx2] = ddata1->type[indx1];
  ddata2->mass[indx2] = ddata1->mass[indx1];
  ddata2->massi[indx2] = ddata1->massi[indx1];
  ddata2->charge[indx2] = ddata1->charge[indx1];

  // IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  // IndexType tid = threadIdx.x;

  // if (bid + tid == 0){
  //   ddata2->numAtom = ddata1->numAtom;
  //   ddata2->numMem  = ddata1->numMem;
  //   ddata2->totalMass = ddata1->totalMass;
  //   ddata2->totalMassi = ddata1->totalMassi;
  // }
}





#endif
