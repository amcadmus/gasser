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


  /// The atom topology.
  struct Atom 
  {
    char name[8];		/**< Name of the atom. */
    ScalorType mass;		/**< Mass of the atom. */
    ScalorType charge;		/**< Charge of the atom. */
    TypeType type;		/**< Type of the atom. */
public:
    Atom ();
    Atom (const Atom & a);
    /** 
     * Constructor
     * 
     * @param mass Mass of the atom.
     * @param charge Charge of the atom.
     * @param type Type of the atom.
     */
    Atom (const ScalorType & mass,
	  const ScalorType & charge,
	  const TypeType & type);
    /** 
     * Set property of the atom.
     * 
     * @param mass Mass of the atom.
     * @param charge Charge of the atom.
     * @param type Type of the atom.
     */
    void setProperty  (const ScalorType & mass,
		       const ScalorType & charge,
		       const TypeType & type);
    /** 
     * Copy
     * 
     * @param a An other atom.
     */
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


  /// Non-bonded interaction.
  class NonBondedInteraction
  {
private:
    TypeType atomType0;
    TypeType atomType1;
    const NonBondedInteractionParameter * ptr_nbInter;
public:
    /** 
     * Constructor, initialize a non-bonded interaction.
     * 
     * @param atomType0 The first atom type
     * @param atomType1 The second atom type, which can be the same as the first one.
     * @param p Parameters of the non-bonded interaction.
     */
    NonBondedInteraction (const TypeType & atomType0,
			  const TypeType & atomType1,
			  const NonBondedInteractionParameter & p);
    /** 
     * Specify the non-bonded interaction.
     * 
     * @param atomType0 The first atom type
     * @param atomType1 The second atom type, which can be the same as the first one.
     * @param p Parameters of the non-bonded interaction.
     */
    void specifyInteraction (const TypeType & atomType0,
			     const TypeType & atomType1,
			     const NonBondedInteractionParameter & p);
    /** 
     * Get the type of the first atom.
     * 
     * @return The type of the first atom.
     */
    const TypeType &	typeOfAtom0 ()	const {return atomType0;}
    /** 
     * Get the type of the second atom.
     * 
     * @return The type of the second atom.
     */
    const TypeType &	typeOfAtom1 ()	const {return atomType1;}
    /** 
     * Get the type of the interaction.
     * 
     * @return The type of the interaction.
     */
    InteractionType	typeOfInter ()	const {return ptr_nbInter->type();}
    /** 
     * Get the number of non-bonded interaction parameters.
     * 
     * @return The number of non-bonded interaction parameters.
     */
    unsigned		numParam ()	const {return ptr_nbInter->numParam();}
    /** 
     * Get the C pointer to the array of non-bonded interaction parameters.
     * 
     * @return The C pointer to the array of non-bonded interaction parameters.
     */
    const ScalorType*   ptr_param ()	const {return ptr_nbInter->c_ptr();}
    /** 
     * Get cut-off radius.
     * 
     * @return The cut-off radius.
     */
    ScalorType	rcut ()		const {return ptr_nbInter->rcut();}
    /** 
     * Get the shift at cut-off.
     * 
     * @return The shift at cut-off.
     */
    ScalorType	shiftAtCut ()	const {return ptr_nbInter->shiftAtCut();}
    /** 
     * Get the energy correction of the non-bonded interaction.
     * 
     * @param rcut The cut-off radius.
     * 
     * @return The energy correction.
     */
    ScalorType	energyCorrection   (const ScalorType & rcut) const {return ptr_nbInter->energyCorrection(rcut);}
    /** 
     * Get the pressure correction of the non-bonded interaction.
     * 
     * @param rcut The cut-off radius.
     * 
     * @return The pressure correction.
     */
    ScalorType	pressureCorrection (const ScalorType & rcut) const {return ptr_nbInter->pressureCorrection(rcut);}
  };


  /// The exclusion list.
  /**
   * Exclude non-bonded interaction between certain atoms within a
   * molecule.
   */
  struct Exclusion 
  {
    IndexType atom0;
    IndexType atom1;
public:
    /** 
     * Constructor.
     * 
     * @param i The index of the first atom in the molecule.
     * @param j The index of the second atom in the molecule.
     */
    Exclusion (const IndexType & i,
	       const IndexType & j);
  };
  
  /// Bond topology.
  struct Bond
  {
    IndexType atom0;
    IndexType atom1;
    InteractionType type;
    std::vector<ScalorType > paramArray;
public:
    /** 
     * Constructor.
     * 
     * @param atom0 The index of the first atom bonded.
     * @param atom1 The index of the second atom bonded.
     * @param p The parameters of the bond interaction.
     */
    Bond (const IndexType & atom0,
	  const IndexType & atom1,
	  const BondInteractionParameter & p);
    /** 
     * Specify a bond.
     * 
     * @param atom0 The index of the first atom bonded.
     * @param atom1 The index of the second atom bonded.
     * @param p The parameters of the bond interaction.
     */
    void specifyInteraction (const IndexType & atom0,
			     const IndexType & atom1,
			     const BondInteractionParameter & p);
  };

  /// Angle topology.
  struct Angle
  {
    IndexType edge0;
    IndexType edge1;
    IndexType center;
    InteractionType type;
    std::vector<ScalorType > paramArray;
public:
    /** 
     * Constructor.
     * 
     * @param edge0 The index of the first atom on edge.
     * @param center The index of the atom on the angle center.
     * @param edge1 The index of the second atom on edge.
     * @param p The parameters of the angle interaction.
     */
    Angle (const IndexType & edge0,
	   const IndexType & center,
	   const IndexType & edge1,
	   const AngleInteractionParameter & p);
    /** 
     * Specify the angle interaction.
     * 
     * @param edge0 The index of the first atom on edge.
     * @param center The index of the atom on the angle center.
     * @param edge1 The index of the second atom on edge.
     * @param p The parameters of the angle interaction.
     */
    void specifyInteraction (const IndexType & edge0,
			     const IndexType & center,
			     const IndexType & edge1,
			     const AngleInteractionParameter & p);
  };

  /// Molecule topology.
  struct Molecule 
  {
    char name[8];
    std::vector<Atom > atoms;
    std::vector<Bond > bonds;
    std::vector<Angle > angles;
    std::vector<Exclusion > exclusions;
public:
    Molecule ();
    /** 
     * Clear everything in the molecule.
     * 
     */
    void clear();
    /** 
     * Push an atom into the molecule.
     * 
     * @param a The atom.
     */
    void pushAtom (const Atom & a);
    /** 
     * Get the size of the molecule.
     * 
     * @return The size.
     */
    IndexType size () const {return atoms.size();}
    /** 
     * Get the i-th atom.
     * 
     * @param i The index of the atom.
     * 
     * @return The atom.
     */
    Atom & operator [] (const IndexType & i) {return atoms[i];}
    /** 
     * Get the i-th atom.
     * 
     * @param i The index of the atom.
     * 
     * @return The atom.
     */
    const Atom & operator [] (const IndexType & i) const {return atoms[i];}
    /** 
     * Add a bond to the molecule.
     * 
     * @param bd The bond.
     */
    void addBond (const Bond & bd);
    /** 
     * Add an angle to the molecule.
     * 
     * @param ag The angle.
     */
    void addAngle (const Angle & ag);
    /** 
     * Add an exclusion to the molecule.
     * 
     * @param exclusion The exclusion.
     */
    void addExclusion (const Exclusion & exclusion);
  };

  /// The system topology.
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
    /** 
     * Add a non-bonded interaction to the system.
     * 
     * @param nb The non-bonded interaction.
     */
    void addNonBondedInteraction (const NonBondedInteraction & nb);
    /** 
     * Add molecules to the system.
     * 
     * @param mol The molecule.
     * @param number Number of this molecule to be added.
     */
    void addMolecules (const Molecule & mol,
		       const IndexType & number);
    /** 
     * Clear everything.
     * 
     */
    void clear();
public:
    /** 
     * Get number of atoms in the system.
     * 
     * @return The number of atoms.
     */
    const IndexType & numAtom () const {return indexShift.back();}
    /** 
     * Get the number of freedom.
     * 
     * @return The number of freedom.
     */
    const IndexType & getNumFreedom () const {return numFreedom;}
    /** 
     * Given the atomic global index, calculate the index of molecule
     * that the atom belongs to and the index of this atom in the
     * molecule.
     * 
     * @param globalIndex The global atomic index.
     * @param top_molIndex The index of the molecule.
     * @param top_atomIndex The index of the atom in the molecule.
     */
    void calMolTopPosition (const IndexType & globalIndex,
			    IndexType & top_molIndex,
			    IndexType & top_atomIndex) const;
    /** 
     * Get atom of certain molecule.
     * 
     * @param top_molIndex The index of molecule.
     * @param top_atomIndex The index of atom in the molecule.
     * 
     * @return The atom.
     */
    const Atom & getAtom (const IndexType & top_molIndex,
			  const IndexType & top_atomIndex) const
	{return molecules[top_molIndex].atoms[top_atomIndex];}
  };
}


#endif
