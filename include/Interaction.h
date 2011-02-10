#ifndef __Interaction_h_wanghan__
#define __Interaction_h_wanghan__


typedef int InteractionType;

class InteractionParamter 
{
public:
  virtual InteractionType type () const= 0;
  virtual unsigned numParam () const = 0;
  virtual ScalorType * c_ptr () = 0;
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


class NonBondedInteractionParameter : public InteractionParamter
{
public:
  NonBondedInteractionParameter () {}
  // NonBondedInteractionParameter (const NonBondedInteractionParameter & p1);
  // const NonBondedInteractionParameter & copy (const NonBondedInteractionParameter & p1);  
  bool same (const NonBondedInteractionParameter & f1) const ;
  bool operator == (const NonBondedInteractionParameter & f1) const;
public:
  virtual ScalorType rcut () const = 0;
public:
  virtual ScalorType shiftAtCut () const;
  virtual ScalorType energyCorrection   (const ScalorType & rcut) const;
  virtual ScalorType pressureCorrection (const ScalorType & rcut) const;  
};

class BondInteractionParameter : public InteractionParamter
{
public:
  BondInteractionParameter () {}
  // BondInteractionParameter (const BondInteractionParameter & p1);
  // const BondInteractionParameter & copy (const BondInteractionParameter & p1);  
  bool same (const BondInteractionParameter & f1) const ;
  bool operator == (const BondInteractionParameter & f1) const;
};

class AngleInteractionParameter : public InteractionParamter
{
public:
  AngleInteractionParameter () {}
  // AngleInteractionParameter (const AngleInteractionParameter & p1);
  // const AngleInteractionParameter & copy (const AngleInteractionParameter & p1);  
  bool same (const AngleInteractionParameter & f1) const ;
  bool operator == (const AngleInteractionParameter & f1) const;
};

#endif
