#ifndef __Interaction_h_wanghan__
#define __Interaction_h_wanghan__


typedef int InteractionType;

/// Type and parameter of a interaction.

class InteractionParamter 
{
public:
  /** 
   * Get type of the interaction.
   * 
   * @return Type of the interaction.
   */
  virtual InteractionType type () const= 0;
  /** 
   * Get number of parameters of the interaction.
   * 
   * @return Number of parameters.
   */
  virtual unsigned numParam () const = 0;
  /** 
   * Get C pointer of the array of parameters.
   * 
   * @return C pointer.
   */  
  virtual ScalorType * c_ptr () = 0;
  /** 
   * Get constant C pointer of the array of parameters.
   * 
   * @return Constant C pointer.
   */  
  virtual const ScalorType * c_ptr () const = 0;
};

enum mdInteraction {
  mdForceNULL				= 00000,
  mdForceLennardJones6_12		= 00001,
  mdForceLennardJones6_12B		= 00005,
  mdForceLennardJones6_12_cap		= 00002,
  mdForceLennardJones6_12B_cap		= 00006,
  mdForceCosTail			= 00003,
  mdForceCosTail_cap			= 00004,
  mdForceHarmonicSpring			= 10000,
  mdForceFENE				= 10001,
  mdForceFENE2				= 10002,
  mdForceAngleHarmonic			= 20000,
  mdForceCosAngle0			= 20001
};

/// Type and parameters of a non-bonded interaction.

class NonBondedInteractionParameter : public InteractionParamter
{
public:
  NonBondedInteractionParameter () {}
  // NonBondedInteractionParameter (const NonBondedInteractionParameter & p1);
  // const NonBondedInteractionParameter & copy (const NonBondedInteractionParameter & p1);
  /** 
   * Test if two non-bonded interactions are the same.
   * 
   * @param f1 The other non-bonded interaction.
   * 
   * @return true if *this and f1 are the same.
   */
  bool same (const NonBondedInteractionParameter & f1) const ;
  /** 
   * Test if two non-bonded interactions are the same.
   * 
   * @param f1 The other non-bonded interaction.
   * 
   * @return true if *this and f1 are the same.
   */
  bool operator == (const NonBondedInteractionParameter & f1) const;
public:
  /** 
   * Get cut-off radius of the short-range non-bonded interaction.
   * 
   * @return The cut-off radius.
   */
  virtual ScalorType rcut () const = 0;
public:
  /** 
   * Get the shift at the cut-off
   * 
   * @return The shift.
   */
  virtual ScalorType shiftAtCut () const;
  /** 
   * Get the energy correction of the short-range non-bonded interaction.
   * 
   * @param rcut Cut-off radius.
   * 
   * @return The energy correction.
   */
  virtual ScalorType energyCorrection   (const ScalorType & rcut) const;
  /** 
   * Get the pressure correction of the short-range non-bonded interaction.
   * 
   * @param rcut Cut-off radius.
   * 
   * @return The pressure correction.
   */
  virtual ScalorType pressureCorrection (const ScalorType & rcut) const;  
};

/// Type and parameters of a bond interaction. 

class BondInteractionParameter : public InteractionParamter
{
public:
  BondInteractionParameter () {}
  // BondInteractionParameter (const BondInteractionParameter & p1);
  // const BondInteractionParameter & copy (const BondInteractionParameter & p1);  
  /** 
   * Test if two bond interactions are the same.
   * 
   * @param f1 The other bond interaction.
   * 
   * @return true if *this and f1 are the same.
   */
  bool same (const BondInteractionParameter & f1) const ;
  /** 
   * Test if two bond interactions are the same.
   * 
   * @param f1 The other bond interaction.
   * 
   * @return true if *this and f1 are the same.
   */
  bool operator == (const BondInteractionParameter & f1) const;
};

/// Type and parameters of an angle interaction.

class AngleInteractionParameter : public InteractionParamter
{
public:
  AngleInteractionParameter () {}
  // AngleInteractionParameter (const AngleInteractionParameter & p1);
  // const AngleInteractionParameter & copy (const AngleInteractionParameter & p1);  
  /** 
   * Test if two angle interactions are the same.
   * 
   * @param f1 The other angle interaction.
   * 
   * @return true if *this and f1 are the same.
   */
  bool same (const AngleInteractionParameter & f1) const ;
  /** 
   * Test if two angle interactions are the same.
   * 
   * @param f1 The other angle interaction.
   * 
   * @return true if *this and f1 are the same.
   */
  bool operator == (const AngleInteractionParameter & f1) const;
};

#endif
