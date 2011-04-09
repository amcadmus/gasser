#ifndef __MDSystem_h_wanghan__
#define __MDSystem_h_wanghan__

#define DEVICE_CODE

#include "systemDefines.h"
// #include "common.h"

class HostMDData
{
public:
  IndexType numAtom;
  IndexType numMem;
  IndexType NFreedom;
  CoordType * coord;
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


void cpyHostMDDataToDevice (const HostMDData * hdata,
			    DeviceMDData * ddata);
void cpyDeviceMDDataToHost (const DeviceMDData * ddata,
			    HostMDData * hdata);
void cpyDeviceMDDataToDevice (const DeviceMDData * ddata1,
			      DeviceMDData * ddata);

__global__ void deviceCpyTypeToCoordW (CoordType * coord,
				       const TypeType * type,
				       const IndexType N);

__host__ void mallocHostMDData (IndexType numAtom,
				IndexType expectedMaxNumAtom,
				HostMDData * hdata);
__host__ void lazyInitHostMDData (HostMDData * hdata);
__host__ void initMass (HostMDData * hdata);
__host__ void destroyHostMDData (HostMDData * hdata);

__host__ void initDeviceMDData (const HostMDData * hdata,
				DeviceMDData * ddata);
__host__ void destroyDeviceMDData (DeviceMDData * ddata);




// __host__ void bindTextureOnDeviceMDData (DeviceMDData * data);

__device__ void
cpyDeviceMDDataElement (const DeviceMDData * ddata1,
			const IndexType indx1,
			DeviceMDData * ddata2,
			const IndexType indx2);


#endif
