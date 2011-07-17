#ifndef __MDSystem_h_wanghan__
#define __MDSystem_h_wanghan__

#define DEVICE_CODE

#include "systemDefines.h"
// #include "common.h"

/// MD data on host.

class HostMDData
{
public:
  IndexType numAtom;		/**< Number of atoms. */
  IndexType numMem;		/**< Size of the memory (counted by
				 * possible number of atoms). */
  IndexType NFreedom;		/**< Degre of freedom. */
  CoordType * coord;		/**< The coordinates of atoms */
  IntScalorType * coordNoix;	/**< Number of images on x. */
  IntScalorType * coordNoiy;	/**< Number of images on y. */
  IntScalorType * coordNoiz;	/**< Number of images on z. */
  ScalorType * velox;		/**< x component of velocity. */
  ScalorType * veloy;		/**< y component of velocity. */
  ScalorType * veloz;		/**< z component of velocity. */
  ScalorType * forcx;		/**< x component of force. */
  ScalorType * forcy;		/**< y component of force. */
  ScalorType * forcz;		/**< z component of force. */
  TypeType   * type;		/**< Type of atom. */
  ScalorType * mass;		/**< Mass of atom. */
  ScalorType * massi;		/**< Inverse mass of atom. */
  ScalorType totalMass;		/**< Total mass of the system. */
  ScalorType totalMassi;	/**< Inverse total mass. */
  ScalorType * charge;		/**< Charge. */
  char * atomName;		/**< Name of the atom. The atom name
				 * is defined in the .gro
				 * configuration file. */
  IndexType * atomIndex;	/**< Index of the atom. Also defined
				 * in the .gro file. */
  char * resdName;		/**< Name of the residue. The residue name
				 * is defined in the .gro
				 * configuration file. */
  IndexType * resdIndex;	/**< Index of the residue. Also
				 * defined in the .gro file. */
public:
  HostMDData() ;
  ~HostMDData() ;
}
    ;

/// System MD data on device.

/**
 * All pointers are public because of the restriction of CUDA.
 * 
 */


class DeviceMDData
{
public:
  bool malloced ;		/**< Tell us if malloced. */
  IndexType numAtom;		/**< Number of atoms. */
  IndexType numMem;		/**< Size of the memory (counted by
				 * possible number of atoms). */
  IndexType NFreedom;		/**< Degree of freedom. */
  CoordType * coord;		/**< The coordinates of atoms. */
  IntScalorType * coordNoix;	/**< Number of images on x. */
  IntScalorType * coordNoiy;	/**< Number of images on y. */
  IntScalorType * coordNoiz;	/**< Number of images on z. */
  ScalorType * velox;		/**< x component of velocity. */
  ScalorType * veloy;		/**< y component of velocity. */
  ScalorType * veloz;		/**< z component of velocity. */
  ScalorType * forcx;		/**< x component of force. */
  ScalorType * forcy;		/**< y component of force. */
  ScalorType * forcz;		/**< z component of force. */
  TypeType   * type;		/**< Type of atom. */
  ScalorType * mass;		/**< Mass of atom. */
  ScalorType * massi;		/**< Inverse mass of atom. */
  ScalorType totalMass;		/**< Total mass of the system. */
  ScalorType totalMassi;	/**< Inverse total mass */
  ScalorType * charge;		/**< Charge */
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
