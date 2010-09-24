#ifndef __WidomTestParticleInsertion_h_wanghan__
#define __WidomTestParticleInsertion_h_wanghan__

#define DEVICE_CODE

#include "common.h"
#include "MDException.h"
#include "MDSystem_interface.h"
#include "NeighborList_interface.h"
#include "SystemNonBondedInteraction.h"

class WidomTestParticleInsertion
{
  bool			inited;
  IndexType		nParticleInserted;
  HostCoordType *	hcoord;
  TypeType *		htype;
  ScalorType 		hresult;
  ScalorType		mytemperature;
  ScalorType		energyCorr;
  ScalorType		scaledEnergyCorr;
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
	       const TypeType & particleType,
	       const SystemNonBondedInteraction & sysNbInter);
  void clear ();
public:
  const IndexType &  numTestParticle () const {return nParticleInserted;}
  const ScalorType & temperature     () const {return mytemperature;}
  // void setEnergyCorr (const ScalorType & c) {energCorr = c;}
  ScalorType energyCorrection () const {return scaledEnergyCorr;}
  void generateTestCoords (const MDSystem & sys);
  ScalorType expMu ();  
};

#endif
