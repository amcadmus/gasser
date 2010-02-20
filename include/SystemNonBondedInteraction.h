#ifndef __SystemNonBondedInteraction_h_wanghan__
#define __SystemNonBondedInteraction_h_wanghan__

#include "common.h"
#include <vector>
#include "Interaction.h"

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
public:
  SystemNonBondedInteraction();
  ~SystemNonBondedInteraction() { clear(); freeAPointer((void**)&positions);}
  void clear ();
  ScalorType maxRcut() {return maxrc;}
  void add (const TypeType &i,
	    const TypeType &j,
	    const NonBondedInteractionParameter & param);
  void build ();
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


#endif

