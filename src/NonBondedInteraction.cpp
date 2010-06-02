#define HOST_CODE
#define CPP_FILE

#include "NonBondedInteraction.h"

bool NonBondedInteractionParameter::
same (const NonBondedInteractionParameter & f1) const
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

bool NonBondedInteractionParameter::
operator == (const NonBondedInteractionParameter & f1) const
{
  return this->same(f1);
}

InteractionType LennardJones6_12Parameter::
type () const
{
  return mdForceLennardJones6_12;
}

unsigned LennardJones6_12Parameter::
numParam () const 
{
  return mdForceNParamLennardJones6_12;
}

const ScalorType * LennardJones6_12Parameter::
c_ptr () const
{
  return param;
}

ScalorType * LennardJones6_12Parameter::
c_ptr () 
{
  return param;
}

ScalorType LennardJones6_12Parameter::
rcut() const
{
  return LennardJones6_12::calRcut (param);
}

void LennardJones6_12Parameter::
reinit (ScalorType epsilon,
	ScalorType sigma,
	ScalorType shift,
	ScalorType rcut)
{
  LennardJones6_12::initParameter (param, epsilon, sigma, shift, rcut);
}

LennardJones6_12Parameter::
LennardJones6_12Parameter (ScalorType epsilon,
			   ScalorType sigma,
			   ScalorType shift,
			   ScalorType rcut)
{
  LennardJones6_12::initParameter (param, epsilon, sigma, shift, rcut);
}

InteractionType LennardJones6_12CapParameter::
type () const 
{
  return mdForceLennardJones6_12_cap;
}

unsigned LennardJones6_12CapParameter::
numParam () const 
{
  return mdForceNParamLennardJones6_12_cap;
}

const ScalorType * LennardJones6_12CapParameter::
c_ptr() const 
{
  return param;
}

ScalorType * LennardJones6_12CapParameter::
c_ptr() 
{
  return param;
}

ScalorType LennardJones6_12CapParameter::
rcut () const
{
  return LennardJones6_12_cap::calRcut (param);
}

void LennardJones6_12CapParameter::
reinit (ScalorType epsilon,
      ScalorType sigma,
      ScalorType shift,
      ScalorType rcut,
      ScalorType cap)
{
  LennardJones6_12_cap::initParameter
      (param, epsilon, sigma, shift, rcut, cap);
}

LennardJones6_12CapParameter::
LennardJones6_12CapParameter (ScalorType epsilon,
			      ScalorType sigma,
			      ScalorType shift,
			      ScalorType rcut,
			      ScalorType cap)
{
  LennardJones6_12_cap::initParameter
      (param, epsilon, sigma, shift, rcut, cap);
}

InteractionType CosTailParameter::
type () const 
{
  return mdForceCosTail;
}

unsigned CosTailParameter::
numParam () const 
{
  return mdForceNParamCosTail;
}

const ScalorType * CosTailParameter::
c_ptr () const 
{
  return param;
}

ScalorType * CosTailParameter::
c_ptr () 
{
  return param;
}

ScalorType CosTailParameter::
rcut () const 
{
  return CosTail::calRcut (param);
}

void  CosTailParameter::
reinit (ScalorType epsilon,
	ScalorType b,
	ScalorType wc)
{
  CosTail::initParameter (param, epsilon, b, wc);
}

CosTailParameter::
CosTailParameter (ScalorType epsilon,
		  ScalorType b,
		  ScalorType wc)
{
  CosTail::initParameter (param, epsilon, b, wc);
}


InteractionType CosTailCapParameter::
type () const 
{
  return mdForceCosTail_cap;
}

unsigned CosTailCapParameter::
numParam () const 
{
  return mdForceNParamCosTail_cap;
}

const ScalorType * CosTailCapParameter::
c_ptr () const 
{
  return param;
}

ScalorType * CosTailCapParameter::
c_ptr () 
{
  return param;
}

ScalorType CosTailCapParameter::
rcut () const 
{
  return CosTail_cap::calRcut (param);
}

void  CosTailCapParameter::
reinit (ScalorType epsilon,
	ScalorType b,
	ScalorType wc,
	ScalorType cap)
{
  CosTail_cap::initParameter (param, epsilon, b, wc, cap);
}

CosTailCapParameter::
CosTailCapParameter (ScalorType epsilon,
		     ScalorType b,
		     ScalorType wc,
		     ScalorType cap)
{
  CosTail_cap::initParameter (param, epsilon, b, wc, cap);
}

  
  
