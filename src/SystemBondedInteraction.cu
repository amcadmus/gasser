#include "SystemBondedInteraction.h"

SystemBondedInteraction::
SystemBondedInteraction()
    : numBondedInteraction (0), memBondedInteraction(0),
      numBondedParameter (0), memBondedParameter (0),
      hasBond_(false), hasAngle_(false),
      type (NULL), bondedParameterPosition (NULL), bondedParameter (NULL)
{
}

SystemBondedInteraction::
SystemBondedInteraction(const Topology::System & sysTop)
    : numBondedInteraction (0), memBondedInteraction(0),
      numBondedParameter (0), memBondedParameter (0),
      hasBond_(false), hasAngle_(false),
      type (NULL), bondedParameterPosition (NULL), bondedParameter (NULL)
{
  reinit (sysTop);
}


SystemBondedInteraction::
~SystemBondedInteraction ()
{
  freeAPointer((void**)&type);
  freeAPointer((void**)&bondedParameterPosition);
  freeAPointer((void**)&bondedParameter);
}

void SystemBondedInteraction::
clear ()
{
  bondIndex.clear();
  bondNeighborIndex.clear();
  angleIndex.clear();
  angleNeighborIndex.clear();
  anglePosi.clear();

  hasBond_ = false;
  hasAngle_ = false;
  numBondedInteraction = 0;
  numBondedParameter = 0;
}

bool SystemBondedInteraction::
findBond (const Topology::Bond & bd, IndexType & i)
{
  for (i = 0; i < numBondedInteraction; ++i){
    if (bd.type == type[i]){
      bool allsame = true;
      for (unsigned j = 0; j < bd.paramArray.size(); ++j){
	if (bd.paramArray[j] != bondedParameter[bondedParameterPosition[i] + j]){
	  allsame = false;
	  break;
	}
      }
      if (allsame){
	return true;
      }
    }
  }
  return false;
}

bool SystemBondedInteraction::
findAngle (const Topology::Angle & ag, IndexType & i)
{
  for (i = angleInteractionShift; i < numBondedInteraction; ++i){
    if (ag.type == type[i]){
      bool allsame = true;
      for (unsigned j = 0; j < ag.paramArray.size(); ++j){
	if (ag.paramArray[j] != bondedParameter[bondedParameterPosition[i] + j]){
	  allsame = false;
	  break;
	}
      }
      if (allsame){
	return true;
      }
    }
  }
  return false;
}

void SystemBondedInteraction::
reallocInteraction ()
{
  type = (InteractionType *) realloc (
      type, sizeof(InteractionType)*memBondedInteraction);
  if (type == NULL){
    throw MDExcptFailedReallocOnHost ("SystemBondedInteraction::reallocInteraction",
				      "type",
				      sizeof(InteractionType)*memBondedInteraction);
  }
  bondedParameterPosition = (IndexType *) realloc (
      bondedParameterPosition, sizeof(InteractionType)*memBondedInteraction);
  if (bondedParameterPosition == NULL){
    throw MDExcptFailedReallocOnHost ("SystemBondedInteraction::reallocInteraction",
				      "bondedParameterPosition",
				      sizeof(IndexType)*memBondedInteraction);
  }
}

void SystemBondedInteraction::
reallocParameter ()
{
  bondedParameter = (ScalorType *) realloc (
      bondedParameter, sizeof(ScalorType)*memBondedParameter);
  if (bondedParameter == NULL){
    throw MDExcptFailedReallocOnHost ("SystemBondedInteraction::reallocParameter",
				      "bondedParameter",
				      sizeof(ScalorType)*memBondedParameter);
  }
}

void SystemBondedInteraction::
reinit (const Topology::System & sysTop)
{
  clear ();
  
  bondIndex.resize(sysTop.molecules.size());
  bondNeighborIndex.resize(sysTop.molecules.size());
  angleIndex.resize(sysTop.molecules.size());
  angleNeighborIndex.resize(sysTop.molecules.size());
  anglePosi.resize(sysTop.molecules.size());
  for (unsigned i = 0; i < sysTop.molecules.size(); ++i){
    bondIndex[i].resize(sysTop.molecules[i].atoms.size());
    bondNeighborIndex[i].resize(sysTop.molecules[i].atoms.size());
    angleIndex[i].resize(sysTop.molecules[i].atoms.size());
    angleNeighborIndex[i].resize(sysTop.molecules[i].atoms.size());
    anglePosi[i].resize(sysTop.molecules[i].atoms.size());
  }
  
  // handel the bonds
  for (unsigned i = 0; i < sysTop.molecules.size(); ++i){
    for (unsigned j = 0; j < sysTop.molecules[i].bonds.size(); ++j){
      hasBond_ = true;
      IndexType interactionIndex;
      if (!findBond (sysTop.molecules[i].bonds[j], interactionIndex)){
	if (int(numBondedInteraction) >= int(memBondedInteraction)-1){
	  memBondedInteraction += 1;
	  memBondedInteraction *= 2;
	  reallocInteraction ();
	}
	interactionIndex = numBondedInteraction ++;
	type[interactionIndex] = sysTop.molecules[i].bonds[j].type;
	bondedParameterPosition[interactionIndex] = numBondedParameter;
	if ((numBondedParameter += sysTop.molecules[i].bonds[j].paramArray.size())
	    > memBondedParameter){
	  memBondedParameter += 10;
	  memBondedParameter *= 2;
	  reallocParameter();
	}
	for (unsigned k = 0;
	     k < sysTop.molecules[i].bonds[j].paramArray.size(); ++k){
	  bondedParameter[k+bondedParameterPosition[interactionIndex]] =
	      sysTop.molecules[i].bonds[j].paramArray[k];
	}
      }
      bondNeighborIndex[i][sysTop.molecules[i].bonds[j].atom0].
	  push_back(sysTop.molecules[i].bonds[j].atom1);
      bondNeighborIndex[i][sysTop.molecules[i].bonds[j].atom1].
	  push_back(sysTop.molecules[i].bonds[j].atom0);
      bondIndex[i][sysTop.molecules[i].bonds[j].atom0].
	  push_back(interactionIndex);
      bondIndex[i][sysTop.molecules[i].bonds[j].atom1].
	  push_back(interactionIndex);
    }
  }

  angleInteractionShift = numBondedInteraction;
  if (hasBond_){
    bondedParameterPosition[numBondedInteraction] = numBondedParameter;
    sortBond();
  }

  // handel angle  
  for (unsigned i = 0; i < sysTop.molecules.size(); ++i){
    for (unsigned j = 0; j < sysTop.molecules[i].angles.size(); ++j){
      hasAngle_ = true;
      IndexType interactionIndex;
      if (!findAngle (sysTop.molecules[i].angles[j], interactionIndex)){
	if (int(numBondedInteraction) >= int(memBondedInteraction)-1){
	  memBondedInteraction += 1;
	  memBondedInteraction *= 2;
	  reallocInteraction ();
	}
	interactionIndex = numBondedInteraction ++;
	type[interactionIndex] = sysTop.molecules[i].angles[j].type;
	bondedParameterPosition[interactionIndex] = numBondedParameter;
	if ((numBondedParameter += sysTop.molecules[i].angles[j].paramArray.size())
	    > memBondedParameter){
	  memBondedParameter += 10;
	  memBondedParameter *= 2;
	  reallocParameter();
	}
	for (unsigned k = 0;
	     k < sysTop.molecules[i].angles[j].paramArray.size(); ++k){
	  bondedParameter[k+bondedParameterPosition[interactionIndex]] =
	      sysTop.molecules[i].angles[j].paramArray[k];
	}
      }
      angleNeighborIndex[i][sysTop.molecules[i].angles[j].edge0].
	  push_back(sysTop.molecules[i].angles[j].center);
      angleNeighborIndex[i][sysTop.molecules[i].angles[j].edge0].
	  push_back(sysTop.molecules[i].angles[j].edge1);
      angleNeighborIndex[i][sysTop.molecules[i].angles[j].center].
	  push_back(sysTop.molecules[i].angles[j].edge0);
      angleNeighborIndex[i][sysTop.molecules[i].angles[j].center].
	  push_back(sysTop.molecules[i].angles[j].edge1);
      angleNeighborIndex[i][sysTop.molecules[i].angles[j].edge1].
	  push_back(sysTop.molecules[i].angles[j].center);
      angleNeighborIndex[i][sysTop.molecules[i].angles[j].edge1].
	  push_back(sysTop.molecules[i].angles[j].edge0);

      anglePosi[i][sysTop.molecules[i].angles[j].edge0].push_back(0);
      anglePosi[i][sysTop.molecules[i].angles[j].center].push_back(1);
      anglePosi[i][sysTop.molecules[i].angles[j].edge1].push_back(2);

      angleIndex[i][sysTop.molecules[i].angles[j].edge0].
	  push_back(interactionIndex);
      angleIndex[i][sysTop.molecules[i].angles[j].center].
	  push_back(interactionIndex);
      angleIndex[i][sysTop.molecules[i].angles[j].edge1].
	  push_back(interactionIndex);
    }
  }
  if (hasAngle_){
    bondedParameterPosition[numBondedInteraction] = numBondedParameter;
    sortAngle();
  }
}

void SystemBondedInteraction::
printEverything ()
{
  if (hasBond_ || hasAngle_){
    printf ("#@@ here are bond parameters\n");
    printf ("#@@ type");
    for (unsigned i = 0; i < numBondedInteraction; ++i){
      printf ("  %d", type[i]);
    }
    printf ("\n");
    printf ("#@@ positions");
    for (unsigned i = 0; i < numBondedInteraction+1; ++i){
      printf ("  %d", bondedParameterPosition[i]);
    }
    printf ("\n");
    printf ("#@@ parameters\n");
    for (unsigned i = 0; i < numBondedInteraction; ++i){
      printf ("#@@@");    
      for (unsigned j = bondedParameterPosition[i];
	   j < bondedParameterPosition[i+1]; ++j){
	printf ("  %f", bondedParameter[j]);
      }
      printf ("\n");
    }
  }
  if (hasBond_){
    printf ("#@@ here are bonds in molecules (bondNeighborIndex, bondIndex)\n");
    for (unsigned i = 0; i < bondIndex.size(); ++i){
      for (unsigned j = 0; j < bondIndex[i].size(); ++j){
	printf ("#@@@ atom %d", j);    
	for (unsigned k = 0; k < bondIndex[i][j].size(); ++k){
	  printf (" (%d %d)", bondNeighborIndex[i][j][k], bondIndex[i][j][k]);
	}
	printf ("\n");
      }
      printf ("\n");
    }
  }
  if (hasAngle_){
    printf ("#@@ here are angles in molecules (angleNeighborIndex0 angleNeighborIndex1, anglePosi, angleIndex)\n");
    for (unsigned i = 0; i < angleIndex.size(); ++i){
      for (unsigned j = 0; j < angleIndex[i].size(); ++j){
	printf ("#@@@ atom %d", j);    
	for (unsigned k = 0; k < angleIndex[i][j].size(); ++k){
	  printf (" (%d %d, %d, %d)",
		  angleNeighborIndex[i][j][2*k],
		  angleNeighborIndex[i][j][2*k+1],
		  anglePosi[i][j][k],
		  angleIndex[i][j][k]);
	}
	printf ("\n");
      }
      printf ("\n");
    }
  }
}

SystemBondedInteraction::sortUnit::
sortUnit  (TypeType type_, IndexType index_)
    : mytype(type_), myindex(index_)
{
}

bool SystemBondedInteraction::sortUnit::
operator < (const sortUnit & su) const
{
  return this->mytype < su.mytype;
}

#include <algorithm>

void SystemBondedInteraction::
sortBond ()
{
  for (unsigned i = 0; i < bondIndex.size(); ++i){
    for (unsigned j = 0; j < bondIndex[i].size(); ++j){
      std::vector<IndexType > bkBondIndex (bondIndex[i][j]);
      std::vector<IndexType > bkBondNeighborIndex (bondNeighborIndex[i][j]);
      std::vector<sortUnit > sus;
      for (unsigned k = 0; k < bkBondIndex.size(); ++k){
	sus.push_back(sortUnit(type[bkBondIndex[k]], k));
      }
      std::sort (sus.begin(), sus.end());
      for (unsigned k = 0; k < bkBondIndex.size(); ++k){
	bondIndex[i][j][k] = bkBondIndex[sus[k].myindex];
	bondNeighborIndex[i][j][k] = bkBondNeighborIndex[sus[k].myindex];
      }
    }
  }
}

void SystemBondedInteraction::
sortAngle ()
{
  for (unsigned i = 0; i < angleIndex.size(); ++i){
    for (unsigned j = 0; j < angleIndex[i].size(); ++j){
      std::vector<IndexType > bkAngleIndex (angleIndex[i][j]);
      std::vector<IndexType > bkAnglePosi (anglePosi[i][j]);
      std::vector<IndexType > bkAngleNeighborIndex (angleNeighborIndex[i][j]);
      std::vector<sortUnit > sus;
      for (unsigned k = 0; k < bkAngleIndex.size(); ++k){
	sus.push_back(sortUnit(type[bkAngleIndex[k]], k));
      }
      std::sort (sus.begin(), sus.end());
      for (unsigned k = 0; k < bkAngleIndex.size(); ++k){
	angleIndex[i][j][k] = bkAngleIndex[sus[k].myindex];
	anglePosi [i][j][k] = bkAnglePosi [sus[k].myindex];
	angleNeighborIndex[i][j][k] = bkAngleNeighborIndex[sus[k].myindex];
      }
    }
  }
}


