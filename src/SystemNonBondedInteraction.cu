#include "SystemNonBondedInteraction.h"


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
  interactionTable_ = NULL;
  isBuilt = false;
  maxrc = 0;
}

SystemNonBondedInteraction::
SystemNonBondedInteraction (const Topology::System & sysTop)
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
  interactionTable_ = NULL;
  isBuilt = false;
  maxrc = 0;

  reinit (sysTop);
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
  freeAPointer ((void**)&interactionTable_);
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

void SystemNonBondedInteraction::
reinit (const Topology::System & sysTop)
{
  clear();
  TypeType maxType = 0;
  // first, cal max type
  for (unsigned i = 0; i < sysTop.molecules.size(); ++i){
    for (unsigned j = 0; j < sysTop.molecules[i].size(); ++j){
      if (sysTop.molecules[i].atoms[j].type > maxType){
	maxType = sysTop.molecules[i].atoms[j].type;
      }
    }
  }
  std::vector<ScalorType > empty;
  add (maxType, maxType, mdForceNULL, empty, 0);
  for (unsigned i = 0; i < sysTop.molecules.size(); ++i){
    for (unsigned j = 0; j < sysTop.molecules[i].nonBondedInteractions.size(); ++j){
      add (sysTop.molecules[i].nonBondedInteractions[j].atomType0,
	   sysTop.molecules[i].nonBondedInteractions[j].atomType1,
	   sysTop.molecules[i].nonBondedInteractions[j].type,
	   sysTop.molecules[i].nonBondedInteractions[j].paramArray,
	   sysTop.molecules[i].nonBondedInteractions[j].rcut);
    }
  }
  build ();
}	   
	   
void SystemNonBondedInteraction::
add (const TypeType &i,
     const TypeType &j,
     const InteractionType & type,
     const std::vector<ScalorType > & paramArray,
     const ScalorType & rcut)
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
  typeMat[ii][jj] = type;
  paramMat[ii][jj] = paramArray;

  if (rcut > maxrc){
    maxrc = rcut;
  }
}


IndexType * SystemNonBondedInteraction::
interactionTableItem (TypeType atom0, TypeType atom1)
{
  IndexType i, j;
  atom0 <= atom1 ?
      (i = IndexType(atom0), j = IndexType(atom1)) :
      (i = IndexType(atom1), j = IndexType(atom0)) ;
  return &interactionTable_[i * numAtomTypes + j - ((i*(i+1)) >> 1)];
}

IndexType * SystemNonBondedInteraction::
interactionTableItem (IndexType atom0, IndexType atom1)
{
  IndexType i, j;
  atom0 <= atom1 ?
      (i = IndexType(atom0), j = IndexType(atom1)) :
      (i = IndexType(atom1), j = IndexType(atom0)) ;
  return &interactionTable_[i * numAtomTypes + j - ((i*(i+1)) >> 1)];
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
  interactionTable_ = (IndexType *) malloc (sizeof(IndexType) *
					   interactionTableLength);
  if (interactionTable_ == NULL){
    throw MDExcptFailedMallocOnHost ("SystemNonBondedInteraction::build",
				     "interactionTable_",
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




