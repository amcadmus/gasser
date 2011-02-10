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


ScalorType LennardJones6_12Parameter::
shiftAtCut () const 
{    
  ScalorType dr2 = param[LennardJones6_12::rcut] * param[LennardJones6_12::rcut];
  ScalorType ri2 = 1.f/dr2;
  ScalorType sri2 = param[LennardJones6_12::sigma] * param[LennardJones6_12::sigma] * ri2;
  ScalorType sri6 = sri2*sri2*sri2;
  return 4.f * param[LennardJones6_12::epsilon] * (sri6*sri6 - sri6) + param[LennardJones6_12::shift];    
  // return 0.;
}

ScalorType LennardJones6_12Parameter::
energyCorrection (const ScalorType & r) const
{
  ScalorType sri = param[LennardJones6_12::sigma] / r;
  ScalorType sri2 = sri * sri;
  ScalorType sri8 = sri2 * sri2;
  sri8 = sri8 * sri8;
  ScalorType s3 = param[LennardJones6_12::sigma];
  s3 = s3 * s3 * s3;	
  return 8.f * M_PI * param[LennardJones6_12::epsilon] * s3 / 3.
      * (sri8 / 3. - sri2 ) * sri;	
  // return 0.;
}  

ScalorType LennardJones6_12Parameter::
pressureCorrection (const ScalorType & r) const
{
  ScalorType sri = param[LennardJones6_12::sigma] / r;
  ScalorType sri2 = sri * sri;
  ScalorType sri8 = sri2 * sri2;
  sri8 = sri8 * sri8;
  ScalorType s3 = param[LennardJones6_12::sigma];
  s3 = s3 * s3 * s3;
  return 16.f * M_PI * param[LennardJones6_12::epsilon] * s3 / 3.
      * (sri8 * 2. / 3. - sri2 ) * sri;    
  // return 0.;
}



InteractionType LennardJones6_12BParameter::
type () const
{
  return mdForceLennardJones6_12B;
}

unsigned LennardJones6_12BParameter::
numParam () const 
{
  return mdForceNParamLennardJones6_12B;
}

const ScalorType * LennardJones6_12BParameter::
c_ptr () const
{
  return param;
}

ScalorType * LennardJones6_12BParameter::
c_ptr () 
{
  return param;
}

ScalorType LennardJones6_12BParameter::
rcut() const
{
  return LennardJones6_12B::calRcut (param);
}

void LennardJones6_12BParameter::
reinit (ScalorType epsilon,
	ScalorType sigma,
	ScalorType shift,
	ScalorType rcut)
{
  LennardJones6_12B::initParameter (param, epsilon, sigma, shift, rcut);
}

LennardJones6_12BParameter::
LennardJones6_12BParameter (ScalorType epsilon,
			    ScalorType sigma,
			    ScalorType shift,
			    ScalorType rcut)
{
  LennardJones6_12B::initParameter (param, epsilon, sigma, shift, rcut);
}

ScalorType LennardJones6_12BParameter::
shiftAtCut () const 
{    
  ScalorType dr2 = param[LennardJones6_12B::rcut] * param[LennardJones6_12B::rcut];
  ScalorType ri2 = 1.f/dr2;
  ScalorType sri2 = param[LennardJones6_12B::sigma] * param[LennardJones6_12B::sigma] * ri2;
  ScalorType sri6 = sri2*sri2*sri2;
  return param[LennardJones6_12B::epsilon] * (sri6*sri6 - 2.f * sri6 + param[LennardJones6_12B::shift]);
  // return 0.;
}

ScalorType LennardJones6_12BParameter::
energyCorrection (const ScalorType & r) const
{
  return 0.;
}  

ScalorType LennardJones6_12BParameter::
pressureCorrection (const ScalorType & r) const
{
  return 0.;
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

ScalorType LennardJones6_12CapParameter::
shiftAtCut () const 
{    
  ScalorType dr2 = param[LennardJones6_12_cap::rcut] * param[LennardJones6_12_cap::rcut];
  ScalorType ri2 = 1.f/dr2;
  ScalorType sri2 = param[LennardJones6_12_cap::sigma] * param[LennardJones6_12_cap::sigma] * ri2;
  ScalorType sri6 = sri2*sri2*sri2;
  return 4.f * param[LennardJones6_12_cap::epsilon] * (sri6*sri6 - sri6) + param[LennardJones6_12_cap::shift];    
}

ScalorType LennardJones6_12CapParameter::
energyCorr (const ScalorType & r) const
{
  ScalorType sri = param[LennardJones6_12_cap::sigma] / r;
  ScalorType sri2 = sri * sri;
  ScalorType sri8 = sri2 * sri2;
  sri8 = sri8 * sri8;
  ScalorType s3 = param[LennardJones6_12_cap::sigma];
  s3 = s3 * s3 * s3;	
  return 8.f * M_PI * param[LennardJones6_12_cap::epsilon] * s3 / 3.
      * (sri8 / 3. - sri2 ) * sri;	
}  

ScalorType LennardJones6_12CapParameter::
pressureCorr (const ScalorType & r) const
{
  ScalorType sri = param[LennardJones6_12_cap::sigma] / r;
  ScalorType sri2 = sri * sri;
  ScalorType sri8 = sri2 * sri2;
  sri8 = sri8 * sri8;
  ScalorType s3 = param[LennardJones6_12_cap::sigma];
  s3 = s3 * s3 * s3;
  return 16.f * M_PI * param[LennardJones6_12_cap::epsilon] * s3 / 3.
      * (sri8 * 2. / 3. - sri2 ) * sri;    
}


InteractionType LennardJones6_12BCapParameter::
type () const 
{
  return mdForceLennardJones6_12B_cap;
}

unsigned LennardJones6_12BCapParameter::
numParam () const 
{
  return mdForceNParamLennardJones6_12B_cap;
}

const ScalorType * LennardJones6_12BCapParameter::
c_ptr() const 
{
  return param;
}

ScalorType * LennardJones6_12BCapParameter::
c_ptr() 
{
  return param;
}

ScalorType LennardJones6_12BCapParameter::
rcut () const
{
  return LennardJones6_12B_cap::calRcut (param);
}

void LennardJones6_12BCapParameter::
reinit (ScalorType epsilon,
	ScalorType sigma,
	ScalorType shift,
	ScalorType rcut,
	ScalorType cap)
{
  LennardJones6_12B_cap::initParameter
      (param, epsilon, sigma, shift, rcut, cap);
}

LennardJones6_12BCapParameter::
LennardJones6_12BCapParameter (ScalorType epsilon,
			      ScalorType sigma,
			      ScalorType shift,
			      ScalorType rcut,
			      ScalorType cap)
{
  LennardJones6_12B_cap::initParameter
      (param, epsilon, sigma, shift, rcut, cap);
}

ScalorType LennardJones6_12BCapParameter::
shiftAtCut () const 
{    
  ScalorType dr2 = param[LennardJones6_12B_cap::rcut] * param[LennardJones6_12B_cap::rcut];
  ScalorType ri2 = 1.f/dr2;
  ScalorType sri2 = param[LennardJones6_12B_cap::sigma] * param[LennardJones6_12B_cap::sigma] * ri2;
  ScalorType sri6 = sri2*sri2*sri2;
  return param[LennardJones6_12B_cap::epsilon] * (sri6*sri6 - 2.f * sri6 + param[LennardJones6_12B_cap::shift]);    
}

ScalorType LennardJones6_12BCapParameter::
energyCorr (const ScalorType & r) const
{
  return 0.f;
}  

ScalorType LennardJones6_12BCapParameter::
pressureCorr (const ScalorType & r) const
{
  return 0.f;
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

  
  
