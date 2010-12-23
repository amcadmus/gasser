#ifndef __Topology_h_wanghan__
#define __Topology_h_wanghan__

#include "common.h"
#include "Interaction.h"
#include <vector> 
#include "MDException.h"

namespace Topology {

  class MDExcptTopology : public MDException {
    char message[MaxExceptionMsgLength];
public:
    MDExcptTopology (const char * description) 
	{strncpy (message, description, MaxExceptionMsgLength);}
    virtual const char* what() const throw()
	{return message;}
  };

  
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

  //   struct NonBondedInteraction
  //   {
  //     TypeType atomType0;
  //     TypeType atomType1;
  //     InteractionType type;
  //     std::vector<ScalorType > paramArray;
  //     ScalorType rcut;
  //     NonBondedInteractionParameter * ptr_param;
  // public:
  //     NonBondedInteraction (const TypeType & atomType0,
  // 			    const TypeType & atomType1,
  // 			    const NonBondedInteractionParameter & p);
  //     void specifyInteraction (const TypeType & atomType0,
  // 			       const TypeType & atomType1,
  // 			       const NonBondedInteractionParameter & p);
  //   };

  class NonBondedInteraction
  {
private:
    TypeType atomType0;
    TypeType atomType1;
    const NonBondedInteractionParameter * ptr_nbInter;
public:
    NonBondedInteraction (const TypeType & atomType0,
			  const TypeType & atomType1,
			  const NonBondedInteractionParameter & p);
    void specifyInteraction (const TypeType & atomType0,
			     const TypeType & atomType1,
			     const NonBondedInteractionParameter & p);
    const TypeType &	typeOfAtom0 ()	const {return atomType0;}
    const TypeType &	typeOfAtom1 ()	const {return atomType1;}
    InteractionType	typeOfInter ()	const {return ptr_nbInter->type();}
    unsigned		numParam ()	const {return ptr_nbInter->numParam();}
    const ScalorType* ptr_param ()	const {return ptr_nbInter->c_ptr();}
    ScalorType	rcut ()		const {return ptr_nbInter->rcut();}
    ScalorType	shiftAtCut ()	const {return ptr_nbInter->shiftAtCut();}
    ScalorType	energyCorrection   (const ScalorType & rcut) const {return ptr_nbInter->energyCorrection(rcut);}
    ScalorType	pressureCorrection (const ScalorType & rcut) const {return ptr_nbInter->pressureCorrection(rcut);}
  };

  struct Exclusion 
  {
    IndexType atom0;
    IndexType atom1;
public:
    Exclusion (const IndexType & i,
	       const IndexType & j);
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
    std::vector<Bond > bonds;
    std::vector<Angle > angles;
    std::vector<Exclusion > exclusions;
public:
    Molecule ();
    void clear();
    void pushAtom (const Atom & a);
    IndexType size () const {return atoms.size();}
    Atom & operator [] (const IndexType & i) {return atoms[i];}
    const Atom & operator [] (const IndexType & i) const {return atoms[i];}
    void addBond (const Bond & bd);
    void addAngle (const Angle & ag);
    void addExclusion (const Exclusion & exclusion);
  };


  struct System
  {
    char name[8];
    IndexType numFreedom;
    std::vector<NonBondedInteraction > nonBondedInteractions;
    std::vector<Molecule  > molecules;
    std::vector<IndexType > numbers;
    std::vector<IndexType > indexShift;
public:
    System();
    void addNonBondedInteraction (const NonBondedInteraction & nb);
    void addMolecules (const Molecule & mol,
		       const IndexType & number);
    void clear();
public:
    const IndexType & numAtom () const {return indexShift.back();}
    const IndexType & getNumFreedom () const {return numFreedom;}
    void calMolTopPosition (const IndexType & globalIndex,
			    IndexType & top_molIndex,
			    IndexType & top_atomIndex) const;
    const Atom & getAtom (const IndexType & top_molIndex,
			  const IndexType & top_atomIndex) const
	{return molecules[top_molIndex].atoms[top_atomIndex];}
  };
}


#endif
