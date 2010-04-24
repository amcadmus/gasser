#ifndef __SystemNonBondedInteraction_h_wanghan__
#define __SystemNonBondedInteraction_h_wanghan__

#include "common.h"
#include <vector>
#include "Interaction.h"
#include "Topology.h"

class SystemNonBondedInteraction
{
private:  
  std::vector<std::vector<std::vector<ScalorType> > > paramMat;
  std::vector<std::vector<InteractionType > > typeMat;
private:
  void resizeMem (IndexType size);
  bool checkIntegrity ();
public:
  IndexType numAtomTypes;
private:
  bool isBuilt;
  IndexType numInteractionItems;
  InteractionType * types;
  ScalorType * parameters;
  IndexType * positions;
  IndexType numParameters;
private:
  IndexType * interactionTable_;
  IndexType interactionTableLength;
private:
  ScalorType maxrc;
  void freeBuilt ();
  IndexType * interactionTableItem (TypeType atom0,
				    TypeType atom1);
  IndexType * interactionTableItem (IndexType atom0,
				    IndexType atom1);
  void calInteractionTableLength ();
  void add (const TypeType &i,
  	    const TypeType &j,
  	    const NonBondedInteractionParameter & param);
  void add (const TypeType &i,
  	    const TypeType &j,
  	    const InteractionType & type,
	    const std::vector<ScalorType > & paramArray,
	    const ScalorType & rcut);
  void build ();
public:
  SystemNonBondedInteraction();
  SystemNonBondedInteraction(const Topology::System & sysTop);
  ~SystemNonBondedInteraction() { clear(); freeAPointer((void**)&positions);}
  void reinit (const Topology::System & sysTop);
  void clear ();
  ScalorType maxRcut() {return maxrc;}
public:
  bool beBuilt () const {return isBuilt;}
  IndexType numberOfAtomTypes () const
      {return numAtomTypes;}
  IndexType numberOfInteraction () const
      {return numInteractionItems;}
  IndexType numberOfParameter () const
      {return numParameters;}
  const InteractionType * interactionType () const
      {return types;}
  const ScalorType * interactionParameter () const
      {return parameters;}
  const IndexType * interactionParameterPosition () const
      {return positions;}
  IndexType interactionTableSize () const
      {return interactionTableLength;}
  const IndexType * interactionTable () const
      {return interactionTable_;}
};


#ifdef DEVICE_CODE
__device__ IndexType
nonBondedInteractionTableItem (IndexType * interactionTable,
			       IndexType numAtomType,
			       TypeType atom0,
			       TypeType atom1)
{
  IndexType i, j;
  atom0 <= atom1 ?
      (i = IndexType(atom0), j = IndexType(atom1)) :
      (i = IndexType(atom1), j = IndexType(atom0)) ;
  return interactionTable[i * numAtomType + j - ((i*(i+1)) >> 1)];
}

__device__ IndexType
calNonBondedInteractionTableLength (IndexType numAtomType)
{
  return ((numAtomType * (numAtomType + 1)) >> 1);
}
#endif

#endif

