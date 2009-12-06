#include "NonBondedInteraction.h"

SystemNBForce::
SystemNBForce()
{
  indexTable.data = NULL;
  indexTable.dataLength = 0;
  indexTable.NAtomType = 0;
  setting.type = NULL;
  setting.NNBForce = 0;
  setting.param = NULL;
  setting.paramPosi = NULL;
  setting.paramLength = 0;
}

void SystemNBForce::
init (const IndexType NAtomType)
{
  indexTable.NAtomType = NAtomType;
  indexTable.dataLength = AtomNBForceTable::calDataLength(NAtomType);
  indexTable.data = (ForceIndexType *) malloc (
      indexTable.dataLength * sizeof(ForceIndexType));
  if (indexTable.data == NULL){
    throw MDExcptFailedMallocOnHost ("SystemNBForce::init", "indexTable.data",
				     indexTable.dataLength * sizeof(ForceIndexType));
  }
  for (IndexType i = 0; i < indexTable.dataLength; ++i){
    indexTable.data[i] = 0;
  }

  memNBForce = 1 + indexTable.dataLength;
  memNBParam = memNBForce * MaxNumberParamPerForce;
  setting.type = (mdNBInteraction_t *) malloc (
      memNBForce * sizeof(mdNBInteraction_t));
  if (setting.type == NULL){
    throw MDExcptFailedMallocOnHost ("SystemNBForce::init", "setting.type",
				     memNBForce * sizeof(mdNBInteraction_t));
  }
  setting.param = (ScalorType *) malloc (memNBParam * sizeof(ScalorType));
  if (setting.param == NULL){
    throw MDExcptFailedMallocOnHost ("SystemNBForce::init", "setting.param",
				     memNBParam * sizeof(ScalorType));
  }
  setting.paramPosi = (IndexType *) malloc (memNBForce * sizeof(IndexType));
  if (setting.paramPosi == NULL){
    throw MDExcptFailedMallocOnHost ("SystemNBForce::init", "setting.paramPosi",
				     memNBForce * sizeof(IndexType));
  }
  setting.NNBForce = 1;
  setting.paramLength = 1;
  setting.type[0] = mdForceNBNull;
  setting.param[0] = 0.f;
  setting.paramPosi[0] = 0;
}

void SystemNBForce::addNBForce (const TypeType &atomi, const TypeType &atomj, 
				const mdNBInteraction_t & forceType,
				const ScalorType * param)
{
  bool exist;
  ForceIndexType i = 0;
  IndexType Nparam = calNumNBParameter (forceType);
  for (i = 0; i < setting.NNBForce; ++i){
    exist = true;
    if (setting.type[i] == forceType ){
      if (forceType == mdForceNBNull) break;
      for (IndexType j = 0; j < Nparam; ++j){
	if (param[j] != setting.param[setting.paramPosi[i]+j]){
	  exist = false;
	  break;
	}
      }
    }
    else {
      exist = false;
    }
    if (exist){
      break;
    }
  }
  if (exist){
    AtomNBForceTable::setTableItem (indexTable.data, indexTable.NAtomType,
				    atomi, atomj, i);
  }
  else{
    AtomNBForceTable::setTableItem (indexTable.data, indexTable.NAtomType,
				    atomi, atomj, i);
    setting.type[setting.NNBForce] = forceType;
    setting.paramPosi[setting.NNBForce] = setting.paramLength;
    for (IndexType j = 0; j < Nparam; ++j){
      setting.param[j + setting.paramLength] = param[j];
    }
    setting.NNBForce ++;
    setting.paramLength += Nparam;
  }
}
  
SystemNBForce::~SystemNBForce()
{
  freeAPointer((void**)&setting.type);
  freeAPointer((void**)&setting.param);
  freeAPointer((void**)&setting.paramPosi);
  freeAPointer((void**)&indexTable.data);
}

  
  
