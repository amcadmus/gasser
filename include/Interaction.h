#ifndef __Interaction_h_wanghan__
#define __Interaction_h_wanghan__


typedef int InteractionType;

class InteractionParamter 
{
public:
  virtual InteractionType type () const= 0;
  virtual unsigned numParam () const = 0;
  virtual const ScalorType * c_ptr () const = 0;
};

enum mdInteraction {
  mdForceNULL				= 00000,
  mdForceLennardJones6_12		= 00001,
  mdForceLennardJones6_12_cap		= 00002,
  mdForceCosTail			= 00003,
  mdForceCosTail_cap			= 00004,
  mdForceHarmonicSpring			= 10000,
  mdForceFENE				= 10001,
  mdForceAngleHarmonic			= 20000
};




#endif
