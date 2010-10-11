#ifndef __WidomTestParticleInsertion_h_wanghan__
#define __WidomTestParticleInsertion_h_wanghan__

#define DEVICE_CODE

#include "common.h"
#include "MDException.h"
#include "MDSystem_interface.h"
#include "NeighborList_interface.h"
#include "SystemNonBondedInteraction.h"
#include "BoxGeometry.h"

using namespace RectangularBoxGeometry;

class WidomTestParticleInsertion_NVT
{
  bool			inited;
  IndexType		nParticleInserted;
  HostCoordType *	hcoord;
  TypeType *		htype;
  ScalorType 		hresult;
  ScalorType		mytemperature;
  ScalorType		energyCorr;
  ScalorType		scaledEnergyCorr;
  IndexType		numAtom;
  ScalorType		volume;
public:
  CoordType *		coordTestParticle;
  TypeType *		typeTestParticle;
  SumVector<ScalorType> sumExpDeltaU;
  ScalorType *		dresult;
public:
  WidomTestParticleInsertion_NVT  ();
  ~WidomTestParticleInsertion_NVT ();
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

class WidomTestParticleInsertion_NPT
{
  bool			inited;
  IndexType		nParticleInserted;
  HostCoordType *	hcoord;
  TypeType *		htype;
  ScalorType 		hresult;
  ScalorType		mytemperature;
  ScalorType		mypressure;
  ScalorType		energyCorr;
  ScalorType		scaledEnergyCorr;
  IndexType		numAtom;
  ScalorType		volume;
public:
  CoordType *		coordTestParticle;
  TypeType *		typeTestParticle;
  SumVector<ScalorType> sumExpDeltaU;
  ScalorType *		dresult;
public:
  WidomTestParticleInsertion_NPT  ();
  ~WidomTestParticleInsertion_NPT ();
public:
  void reinit (const ScalorType & temperature,
	       const ScalorType & pressure,
	       const IndexType & nParticleInserted,
	       const TypeType & particleType,
	       const SystemNonBondedInteraction & sysNbInter);
  void clear ();
public:
  const IndexType &  numTestParticle () const {return nParticleInserted;}
  const ScalorType & temperature     () const {return mytemperature;}
  const ScalorType & pressure	     () const {return mypressure;}
  // void setEnergyCorr (const ScalorType & c) {energCorr = c;}
  ScalorType energyCorrection () const {return scaledEnergyCorr;}
  void generateTestCoords (const MDSystem & sys);
  ScalorType expMu ();  
};


// use center point integral formula rather than Monte Carlo
class WidomTestParticleInsertion_NVT2
{
  bool			inited;
  IndexType		nParticleInserted;
  HostCoordType *	hcoord;
  TypeType *		htype;
  ScalorType 		hresult;
  ScalorType		mytemperature;
  ScalorType		scaledEnergyCorr;
  IndexType		numAtom;
  ScalorType		volume;
  ScalorType		inteCellVolume;
public:
  CoordType *		coordTestParticle;
  TypeType *		typeTestParticle;
  SumVector<ScalorType> sumExpDeltaU;
  ScalorType *		dresult;
public:
  WidomTestParticleInsertion_NVT2  ();
  ~WidomTestParticleInsertion_NVT2 ();
public:
  void reinit (const ScalorType & temperature_,
	       const ScalorType & gridSize,
	       const RectangularBox & box,
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


class WidomTestParticleInsertion_NPT2
{
  bool			inited;
  IndexType		nParticleInserted;
  HostCoordType *	hcoord;
  TypeType *		htype;
  TypeType		particleType;
  ScalorType 		hresult;
  ScalorType		mytemperature;
  ScalorType		mypressure;
  IndexType		numAtom;
  ScalorType		volume;
  ScalorType		energyCorr;
  ScalorType		scaledEnergyCorr;
  ScalorType		gridSize;
  ScalorType		inteCellVolume;
public:
  CoordType *		coordTestParticle;
  TypeType *		typeTestParticle;
  SumVector<ScalorType> sumExpDeltaU;
  ScalorType *		dresult;
public:
  WidomTestParticleInsertion_NPT2  ();
  ~WidomTestParticleInsertion_NPT2 ();
public:
  void reinit (const ScalorType & temperature,
	       const ScalorType & pressure,
	       const ScalorType & gridSize,
	       const TypeType & particleType,
	       const SystemNonBondedInteraction & sysNbInter);
  void clear ();
public:
  const IndexType &  numTestParticle () const {return nParticleInserted;}
  const ScalorType & temperature     () const {return mytemperature;}
  const ScalorType & pressure	     () const {return mypressure;}
  // void setEnergyCorr (const ScalorType & c) {energCorr = c;}
  ScalorType energyCorrection () const {return scaledEnergyCorr;}
  void generateTestCoords (const MDSystem & sys);
  ScalorType expMu ();  
};


#endif
