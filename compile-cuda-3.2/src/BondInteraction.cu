#include "BondInteraction.h"

bool BondInteractionParameter::
same (const BondInteractionParameter & f1) const
{
  if (f1.type() != this->type()){
    return false;
  }
  const ScalorType * myparam = this->c_ptr();
  const ScalorType * fparam  = f1.c_ptr();
  for (unsigned i = 0; i < this->numParam(); ++i){
    if (myparam[i] != fparam[i]){
      return false;
    }
  }
  return true;
}

bool BondInteractionParameter::
operator == (const BondInteractionParameter & f1) const
{
  return this->same(f1);
}

// BondInteractionParameter::
// BondInteractionParameter(const BondInteractionParameter & p1)
// {
//   this->copy(p1);
// }

// const BondInteractionParameter & BondInteractionParameter::
// copy (const BondInteractionParameter & p1)
// {
//   ScalorType * myparam = this->c_ptr();
//   const ScalorType * fparam  = p1.c_ptr();
//   for (unsigned i = 0; i < p1.numParam(); ++i){
//     myparam[i] = fparam[i];
//   }
//   return *this;
// }


InteractionType HarmonicSpringParameter::
type () const 
{
  return mdForceHarmonicSpring;
}

unsigned HarmonicSpringParameter::
numParam () const 
{
  return mdForceNParamHarmonicSpring;
}

const ScalorType * HarmonicSpringParameter::
c_ptr () const 
{
  return param;
}

ScalorType * HarmonicSpringParameter::
c_ptr ()  
{
  return param;
}

void HarmonicSpringParameter::
reinit (ScalorType k, ScalorType r0)
{
  HarmonicSpring::initParameter (param, k, r0);
}

HarmonicSpringParameter::
HarmonicSpringParameter (ScalorType k, ScalorType r0)
{
  HarmonicSpring::initParameter (param, k, r0);
}


InteractionType FENEParameter::
type () const 
{
  return mdForceFENE;
}

unsigned FENEParameter::
numParam () const 
{
  return mdForceNParamFENE;
}

const ScalorType * FENEParameter::
c_ptr () const 
{
  return param;
}

ScalorType * FENEParameter::
c_ptr () 
{
  return param;
}

void FENEParameter::
reinit (ScalorType k, ScalorType rinf)
{
  FENE::initParameter (param, k, rinf);
}


FENEParameter::
FENEParameter (ScalorType k, ScalorType rinf)
{
  FENE::initParameter (param, k, rinf);
}




InteractionType FENE2Parameter::
type () const 
{
  return mdForceFENE2;
}

unsigned FENE2Parameter::
numParam () const 
{
  return mdForceNParamFENE2;
}

const ScalorType * FENE2Parameter::
c_ptr () const 
{
  return param;
}

ScalorType * FENE2Parameter::
c_ptr () 
{
  return param;
}

void FENE2Parameter::
reinit (ScalorType k, ScalorType rs, ScalorType r0)
{
  FENE2::initParameter (param, k, rs, r0);
}


FENE2Parameter::
FENE2Parameter (ScalorType k, ScalorType rs, ScalorType r0)
{
  FENE2::initParameter (param, k, rs, r0);
}



