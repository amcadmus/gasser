#include "AngleInteraction.h"

bool AngleInteractionParameter::
same (const AngleInteractionParameter & f1) const
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

bool AngleInteractionParameter::
operator == (const AngleInteractionParameter & f1) const
{
  return this->same(f1);
}

// AngleInteractionParameter::
// AngleInteractionParameter(const AngleInteractionParameter & p1)
// {
//   this->copy(p1);
// }

// const AngleInteractionParameter & AngleInteractionParameter::
// copy (const AngleInteractionParameter & p1)
// {
//   ScalorType * myparam = this->c_ptr();
//   const ScalorType * fparam  = p1.c_ptr();
//   for (unsigned i = 0; i < p1.numParam(); ++i){
//     myparam[i] = fparam[i];
//   }
//   return *this;
// }

    
InteractionType AngleHarmonicParameter::
type () const 
{
  return mdForceAngleHarmonic;
}

unsigned AngleHarmonicParameter::
numParam () const 
{
  return mdForceNParamAngleHarmonic;
}

const ScalorType * AngleHarmonicParameter::
c_ptr () const 
{
  return param;
}

ScalorType * AngleHarmonicParameter::
c_ptr ()  
{
  return param;
}

void AngleHarmonicParameter::
reinit (ScalorType k, ScalorType theta0)
{
  AngleHarmonic::initParameter (param, k, theta0);
}

AngleHarmonicParameter::
AngleHarmonicParameter (ScalorType k, ScalorType theta0)
{
  AngleHarmonic::initParameter (param, k, theta0);
}






    
InteractionType CosAngle0Parameter::
type () const 
{
  return mdForceCosAngle0;
}

unsigned CosAngle0Parameter::
numParam () const 
{
  return mdForceNParamCosAngle0;
}

const ScalorType * CosAngle0Parameter::
c_ptr () const 
{
  return param;
}

ScalorType * CosAngle0Parameter::
c_ptr ()  
{
  return param;
}

void CosAngle0Parameter::
reinit (ScalorType k)
{
  CosAngle0::initParameter (param, k);
}

CosAngle0Parameter::
CosAngle0Parameter (ScalorType k)
{
  CosAngle0::initParameter (param, k);
}




