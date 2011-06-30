#define HOST_CODE
#define CPP_FILE

#include "common.h"
#include "Interaction.h"

ScalorType NonBondedInteractionParameter::
shiftAtCut () const 
{
  return 0.f;
}

ScalorType NonBondedInteractionParameter::
energyCorrection (const ScalorType & rcut) const 
{
  return 0.f;
}

ScalorType NonBondedInteractionParameter::
pressureCorrection (const ScalorType & rcut) const 
{
  return 0.f;
}
