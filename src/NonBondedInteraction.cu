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

ScalorType LennardJones6_12Parameter::
rcut() const
{
  return LennardJones6_12::calRcut (param);
}

void LennardJones6_12Parameter::
init (ScalorType epsilon,
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

ScalorType LennardJones6_12CapParameter::
rcut () const
{
  return LennardJones6_12_cap::calRcut (param);
}

void LennardJones6_12CapParameter::
init (ScalorType epsilon,
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

ScalorType CosTailParameter::
rcut () const 
{
  return CosTail::calRcut (param);
}

void  CosTailParameter::
init (ScalorType epsilon,
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

ScalorType CosTailCapParameter::
rcut () const 
{
  return CosTail_cap::calRcut (param);
}

void  CosTailCapParameter::
init (ScalorType epsilon,
      ScalorType b,
      ScalorType wc,
      ScalorType cap)
{
  CosTail_cap::initParameter (param, epsilon, b, wc, cap);
}




void SystemNonBondedInteraction::
resizeMem (IndexType size)
{
  paramMat.resize (unsigned(size));
  typeMat.resize (unsigned(size));
  for (IndexType i = 0; i < unsigned(size); ++i){
    paramMat[i].resize(unsigned(size));
    typeMat[i].resize(unsigned(size), mdForceNULL);
  }
}

SystemNonBondedInteraction::
SystemNonBondedInteraction ()
{
  numParameters = 0;
  numInteractionItems = 0;
  numAtomTypes = 0;
  interactionTableLength = 0;
  types = NULL;
  parameters = NULL;
  positions = (IndexType *) malloc (sizeof(IndexType));
  if (positions == NULL){
    throw MDExcptFailedMallocOnHost ("SystemNonBondedInteraction::SystemNonBondedInteraction",
				     "positions",
				     sizeof(IndexType) );
  }
  positions[0] = 0;
  interactionTable = NULL;
  isBuilt = false;
  maxrc = 0;
}

void SystemNonBondedInteraction::
freeBuilt ()
{
  freeAPointer ((void**)&types);
  freeAPointer ((void**)&parameters);
  freeAPointer ((void**)&positions);
  positions = (IndexType *) malloc (sizeof(IndexType));
  if (positions == NULL){
    throw MDExcptFailedMallocOnHost ("SystemNonBondedInteraction::SystemNonBondedInteraction",
				     "positions",
				     sizeof(IndexType) );
  }
  positions[0] = 0;
  freeAPointer ((void**)&interactionTable);
  numParameters = 0;
  numInteractionItems = 0;
  interactionTableLength = 0;
}

void SystemNonBondedInteraction::
clear ()
{
  paramMat.clear();
  typeMat.clear();
  numAtomTypes = 0;
  if (isBuilt){
    freeBuilt();
    isBuilt = false;
  }
  maxrc = 0;
}

void SystemNonBondedInteraction::
add (const TypeType &i,
     const TypeType &j,
     const NonBondedInteractionParameter & param)
{
  IndexType ii, jj;
  (i < j) ?
      (ii = IndexType(i), jj = IndexType(j)) :
      (ii = IndexType(j), jj = IndexType(i)) ;
  if (jj+1 > numAtomTypes){
    numAtomTypes = jj+1;
    resizeMem (numAtomTypes);
  }
  if (typeMat[ii][jj] != mdForceNULL){
    printf ("## warning: SystemNonBondedInteraction::add: atome type %d and %d, overwriting existing non-bonded interaction\n", ii, jj);
  }
  typeMat[ii][jj] = param.type();
  paramMat[ii][jj].resize (param.numParam());
  const ScalorType * ptr = param.c_ptr();
  for (unsigned k = 0; k < param.numParam(); ++k){
    paramMat[ii][jj][k] = ptr[k];
  }

  if (param.rcut() > maxrc){
    maxrc = param.rcut();
  }
}

IndexType * SystemNonBondedInteraction::
interactionTableItem (TypeType atom0, TypeType atom1)
{
  IndexType i, j;
  atom0 <= atom1 ?
      (i = IndexType(atom0), j = IndexType(atom1)) :
      (i = IndexType(atom1), j = IndexType(atom0)) ;
  return &interactionTable[i * numAtomTypes + j - ((i*(i+1)) >> 1)];
}
IndexType * SystemNonBondedInteraction::
interactionTableItem (IndexType atom0, IndexType atom1)
{
  IndexType i, j;
  atom0 <= atom1 ?
      (i = IndexType(atom0), j = IndexType(atom1)) :
      (i = IndexType(atom1), j = IndexType(atom0)) ;
  return &interactionTable[i * numAtomTypes + j - ((i*(i+1)) >> 1)];
}

void SystemNonBondedInteraction::
calInteractionTableLength ()
{
  interactionTableLength = ((numAtomTypes * (numAtomTypes + 1)) >> 1);
}



void SystemNonBondedInteraction::
build ()
{
  if (isBuilt){
    freeBuilt();
    isBuilt = false;
  }
  calInteractionTableLength();
  interactionTable = (IndexType *) malloc (sizeof(IndexType) *
					   interactionTableLength);
  if (interactionTable == NULL){
    throw MDExcptFailedMallocOnHost ("SystemNonBondedInteraction::build",
				     "interactionTable",
				     sizeof(IndexType) * interactionTableLength);
  }
  for (unsigned i = 0; i < unsigned(numAtomTypes); ++i){
    for (unsigned j = i; j < unsigned(numAtomTypes); ++j){
      bool exist = false;
      unsigned findex = 0;
      for (unsigned k = 0; k < numInteractionItems; ++k){
	if (typeMat[i][j] == types[k]){
	  if (types[k] != mdForceNULL){
	    bool isMatch = true;
	    for (unsigned m = positions[k]; m < positions[k+1]; ++m){
	      if (paramMat[i][j][m-positions[k]] != parameters[m]){
		isMatch = false;
		break;
	      }
	    }
	    if (isMatch){
	      exist = true;
	      findex = k;
	      break;
	    }
	  }	  
	  else {
	    exist = true;
	    findex = k;
	    break;
	  }
	}
      }
      if (exist){
	(*interactionTableItem(i,j)) = findex;
      }
      else if (typeMat[i][j] != mdForceNULL) {
	numInteractionItems ++;
	// printf ("add new\n");

	types = (InteractionType *) realloc (
	    (void*)types, sizeof(InteractionType) * numInteractionItems);
	if (types == NULL){
	  throw MDExcptFailedReallocOnHost ("SystemNonBondedInteraction::build",
					    "types",
					    sizeof(InteractionType) *
					    numInteractionItems);
	}
	positions = (IndexType *) realloc (
	    (void*)positions, sizeof(IndexType) * (numInteractionItems+1));
	if (positions == NULL){
	  throw MDExcptFailedReallocOnHost ("SystemNonBondedInteraction::build",
					    "positions",
					    sizeof(IndexType) *
					    numInteractionItems+1);
	}
	numParameters += paramMat[i][j].size();
	parameters = (ScalorType *) realloc (
	    (void*)parameters, sizeof(ScalorType) * numParameters);
	if (parameters == NULL){
	  throw MDExcptFailedReallocOnHost ("SystemNonBondedInteraction::build",
					    "parameters",
					    sizeof(ScalorType) *
					    numParameters);
	}
	types[numInteractionItems-1] = typeMat[i][j];
	positions[numInteractionItems] = numParameters;
	for (unsigned k = positions[numInteractionItems-1];
	     k < positions[numInteractionItems]; ++k){
	  parameters[k] = paramMat[i][j][k-positions[numInteractionItems-1]];
	}
	(*interactionTableItem(i,j)) = numInteractionItems-1;
      }
      else {
	numInteractionItems ++;
	// printf ("add new null\n");
	
	types = (InteractionType *) realloc (
	    (void*)types, sizeof(InteractionType) * numInteractionItems);
	if (types == NULL){
	  throw MDExcptFailedReallocOnHost ("SystemNonBondedInteraction::build",
					    "types",
					    sizeof(InteractionType) *
					    numInteractionItems);
	}
	positions = (IndexType *) realloc (
	    (void*)positions, sizeof(IndexType) * (numInteractionItems+1));
	if (positions == NULL){
	  throw MDExcptFailedReallocOnHost ("SystemNonBondedInteraction::build",
					    "positions",
					    sizeof(IndexType) *
					    numInteractionItems+1);
	}
	numParameters += 1;
	parameters = (ScalorType *) realloc (
	    (void*)parameters, sizeof(ScalorType) * numParameters);
	if (parameters == NULL){
	  throw MDExcptFailedReallocOnHost ("SystemNonBondedInteraction::build",
					    "parameters",
					    sizeof(ScalorType) *
					    numParameters);
	}
	types[numInteractionItems-1] = typeMat[i][j];
	positions[numInteractionItems] = numParameters;
	for (unsigned k = positions[numInteractionItems-1];
	     k < positions[numInteractionItems]; ++k){
	  parameters[k] = 0;
	}
	(*interactionTableItem(i,j)) = numInteractionItems-1;
      }
    }
  }
  isBuilt = true;
}




  
  
