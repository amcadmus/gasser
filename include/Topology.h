#ifndef __Topology_h_wanghan__
#define __Topology_h_wanghan__

#include "common.h"
#include "Interaction.h"
#include <vector> 

namespace Topology {
    
    struct Atom 
    {
      char name[8];
      ScalorType mass;
      ScalorType charge;
      TypeType type;
  public:
      Atom ();
      Atom (const Atom & a);
      Atom (const ScalorType & mass,
	    const ScalorType & charge,
	    const TypeType & type);
      void setProperty  (const ScalorType & mass,
			 const ScalorType & charge,
			 const TypeType & type);
      void copy (const Atom & a);
    };

    struct NonBondedInteraction
    {
      TypeType atomType0;
      TypeType atomType1;
      InteractionType type;
      std::vector<ScalorType > paramArray;
      ScalorType rcut;
  public:
      NonBondedInteraction (const TypeType & atomType0,
			    const TypeType & atomType1,
			    const NonBondedInteractionParameter & p);
      void specifyInteraction (const TypeType & atomType0,
			       const TypeType & atomType1,
			       const NonBondedInteractionParameter & p);
    };
    
    struct Bond
    {
      IndexType atom0;
      IndexType atom1;
      InteractionType type;
      std::vector<ScalorType > paramArray;
  public:
      Bond (const IndexType & atom0,
	    const IndexType & atom1,
	    const BondInteractionParameter & p);
      void specifyInteraction (const IndexType & atom0,
			       const IndexType & atom1,
			       const BondInteractionParameter & p);
    };

    struct Angle
    {
      IndexType edge0;
      IndexType edge1;
      IndexType center;
      InteractionType type;
      std::vector<ScalorType > paramArray;
  public:
      Angle (const IndexType & edge0,
	     const IndexType & center,
	     const IndexType & edge1,
	     const AngleInteractionParameter & p);
      void specifyInteraction (const IndexType & edge0,
			       const IndexType & center,
			       const IndexType & edge1,
			       const AngleInteractionParameter & p);
    };
    
    struct Molecule 
    {
      char name[8];
      std::vector<Atom > atoms;
      std::vector<NonBondedInteraction > nonBondedInteractions;
      std::vector<Bond > bonds;
      std::vector<Angle > angles;
  public:
      Molecule ();
      void clear();
      void pushAtom (const Atom & a);
      IndexType size () const {return atoms.size();}
      Atom & operator [] (const IndexType & i) {return atoms[i];}
      const Atom & operator [] (const IndexType & i) const {return atoms[i];}
      void addBond (const Bond & bd);
      void addAngle (const Angle & ag);
      void addNonBondedInteraction (const NonBondedInteraction & nb);
    };


    struct System
    {
      char name[8];
      std::vector<Molecule > molecules;
      std::vector<IndexType > numbers;
      std::vector<IndexType > indexShift;
  public:
      System();
      void addMolecules (const Molecule & mol,
			 const IndexType & number);
      void clear();
    };
}


#endif
