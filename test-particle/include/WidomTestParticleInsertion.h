#ifndef __WidomTestParticleInsertion_h_wanghan__
#define __WidomTestParticleInsertion_h_wanghan__

#define DEVICE_CODE

#include "common.h"
#include "MDException.h"
#include "MDSystem_interface.h"
#include "NeighborList_interface.h"

class WidomTestParticleInsertion
{
  bool			inited;
  IndexType		nParticleInserted;
  HostCoordType *	hcoord;
  TypeType *		htype;
  ScalorType 		hresult;
  ScalorType		mytemperature;
public:
  CoordType *		coordTestParticle;
  TypeType *		typeTestParticle;
  SumVector<ScalorType> sumExpDeltaU;
  ScalorType *		dresult;
public:
  WidomTestParticleInsertion  ();
  ~WidomTestParticleInsertion ();
public:
  void reinit (const ScalorType & temperature,
	       const IndexType & nParticleInserted,
	       const TypeType & particleType);
  void clear ();
public:
  const IndexType &  numTestParticle () const {return nParticleInserted;}
  const ScalorType & temperature     () const {return mytemperature;}
  void generateTestCoords (const MDSystem & sys);
  ScalorType expMu ();  
};

#endif
