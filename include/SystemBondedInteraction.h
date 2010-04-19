#ifndef __SystemBondedInteraction_h_wanghan__
#define __SystemBondedInteraction_h_wanghan__

#include "common.h"
#include "Topology.h"
#include <vector>

class SystemBondedInteraction
{
  struct sortUnit 
  {
    TypeType mytype;
    IndexType myindex;
    sortUnit (TypeType type, IndexType index);
    bool operator < (const sortUnit & su) const;
  };
public:
  std::vector<std::vector<std::vector<IndexType > > > bondIndex;
  std::vector<std::vector<std::vector<IndexType > > > bondNeighborIndex;
  std::vector<std::vector<std::vector<IndexType > > > angleIndex;
  std::vector<std::vector<std::vector<IndexType > > > angleNeighborIndex;
  std::vector<std::vector<std::vector<IndexType > > > anglePosi;
  // HostBondList   hdblist;
  // DeviceBondList dbdlist;
  unsigned angleInteractionShift;
public:
  const std::vector<IndexType > &
  getTopBondNeighborIndex (const IndexType & topMolIndex,
			   const IndexType & topAtomIndex) const
      {return bondNeighborIndex[topMolIndex][topAtomIndex];}
  const std::vector<IndexType > &
  getTopBondIndex (const IndexType & topMolIndex,
		   const IndexType & topAtomIndex) const
      {return bondIndex[topMolIndex][topAtomIndex];}
private:
  bool hasBond_;
  bool hasAngle_;
  bool findBond (const Topology::Bond & bd, IndexType & i);
  bool findAngle (const Topology::Angle & ag, IndexType & i);
  void reallocInteraction ();
  void reallocParameter ();
  void printEverything ();
  void sortBond ();
  void sortAngle ();
public:
  IndexType numBondedInteraction;
  IndexType memBondedInteraction;
  InteractionType * type;
  IndexType  * bondedParameterPosition;
  ScalorType * bondedParameter;
  IndexType numBondedParameter;
  IndexType memBondedParameter;
public:
  SystemBondedInteraction ();
  ~SystemBondedInteraction ();
  void reinit (const Topology::System & sysTop);
  void clear ();
public:
  bool hasBond () const {return hasBond_;}
  bool hasAngle () const {return hasAngle_;}
  IndexType numberOfInteraction () const
      {return numBondedInteraction;}
  IndexType numberOfParameter () const
      {return numBondedParameter;}
  const InteractionType * interactionType () const
      {return type;}
  const ScalorType * interactionParameter () const
      {return bondedParameter;}
  const IndexType * interactionParameterPosition () const
      {return bondedParameterPosition;}
};

#endif
