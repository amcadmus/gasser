#include "BondList_interface.h"

BondList::BondList ()
{
  bondType = NULL;
  paramPosi = NULL;
  NBondForce = 0;
  param = NULL;
  paramLength = 0;
  initDeviceBondList (dbdlist);
}

BondList::~BondList()
{
  freeAPointer ((void**)&bondType);
  freeAPointer ((void**)&paramPosi);
  freeAPointer ((void**)&param);
  destroyDeviceBondList (dbdlist);
}


void BondList::init (const DeviceMDData & ddata)
{
  hbdlist.init (ddata.numAtom);

  NBondForce_mem = 1024;
  paramLength_mem = NBondForce_mem * 3;
  bondType = (mdBondInteraction_t *) realloc (
      bondType, sizeof(mdBondInteraction_t) * NBondForce_mem);
  if (bondType == NULL) {
    throw MDExcptFailedReallocOnHost ("BondList::init", "bondType",
				      sizeof(mdBondInteraction_t) * NBondForce_mem);
  }
  paramPosi = (IndexType *) realloc (
      paramPosi, sizeof(mdBondInteraction_t) * NBondForce_mem);
  if (paramPosi == NULL){
    throw MDExcptFailedReallocOnHost ("BondList::init", "paramPosi",
				      sizeof(mdBondInteraction_t) * NBondForce_mem);
  }				      
  param = (ScalorType *) realloc (
      param, sizeof(ScalorType) * paramLength_mem);
  if (param == NULL){
    throw MDExcptFailedReallocOnHost ("BondList::init", "param",
				      sizeof(ScalorType) * paramLength_mem);
  }
}

  
void BondList::addBond (const IndexType & ii, const IndexType & jj,
			const mdBondInteraction_t & type,
			const ScalorType * thisparam)
{
  bool exist = false;
  ForceIndexType looking;;
  IndexType NParam = calNumBondParameter (type);
  for (looking = 0; looking < NBondForce; ++looking){
    if (type == bondType[looking]){
      bool same = true;
      for (IndexType i = 0; i < NParam; ++i){
	if (thisparam[i] != param[paramPosi[looking] + i]){
	  same = false;
	  break;
	}
      }
      if (same){
	exist = true;
	break;
      }
    }
  }
  if (!exist){
    if (NBondForce == NBondForce_mem){
      NBondForce_mem *= 2;
      bondType = (mdBondInteraction_t *) realloc (
	  bondType, sizeof(mdBondInteraction_t) * NBondForce_mem);
      if (bondType == NULL) {
	throw MDExcptFailedReallocOnHost ("BondList::init", "bondType",
					  sizeof(mdBondInteraction_t) * NBondForce_mem);
      }
      paramPosi = (IndexType *) realloc (
	  paramPosi, sizeof(mdBondInteraction_t) * NBondForce_mem);
      if (paramPosi == NULL){
	throw MDExcptFailedReallocOnHost ("BondList::init", "paramPosi",
					  sizeof(mdBondInteraction_t) * NBondForce_mem);
      }				      
    }
    if (paramLength == paramLength_mem){
      paramLength_mem *= 2;
      param = (ScalorType *) realloc (
	  param, sizeof(ScalorType) * paramLength_mem);
      if (param == NULL){
	throw MDExcptFailedReallocOnHost ("BondList::init", "param",
					  sizeof(ScalorType) * paramLength_mem);
      }
    } 
    bondType [NBondForce] = type;
    paramPosi[NBondForce] = paramLength;
    for (IndexType i = 0; i < NParam; ++i){
      param[paramPosi[NBondForce] + i] = thisparam[i];
    }
    NBondForce ++;
    paramLength += NParam;
  }
  hbdlist.addBond (ii, jj, looking);
}


// bubble sorting
// void BondList::sortBond()
// {
//   IndexType *indexMap = (IndexType *)malloc (sizeof(IndexType) * hbdlist.listLength);
//   if (indexMap == NULL){
//     MDExcptFailedMallocOnHost ("BondList::sortBond", "indexMap",
// 			       sizeof(IndexType) * hbdlist.listLength);
//   }
//   TypeType  *typeBuff = (TypeType *) malloc (sizeof(TypeType)  * hbdlist.listLength);
//   if (typeBuff == NULL){
//     MDExcptFailedMallocOnHost ("BondList::sortBond", "typeBuff",
// 			       sizeof(TypeType)  * hbdlist.listLength);
//   }
//   ForceIndexType * bkForceIndex = (ForceIndexType *) malloc (
//       sizeof (ForceIndexType) * hbdlist.listLength);
//   if (bkForceIndex == NULL){
//     MDExcptFailedMallocOnHost ("BondList::sortBond", "bkForceIndex",
// 			       sizeof (ForceIndexType) * hbdlist.listLength);
//   }			     
//   IndexType * bkData = (IndexType *) malloc (
//       sizeof (IndexType) * hbdlist.listLength);
//   if (bkData == NULL){
//     MDExcptFailedMallocOnHost ("BondList::sortBond", "bkData",
// 			       sizeof (IndexType) * hbdlist.listLength);
//   }
//   for (IndexType i = 0; i < hbdlist.stride; ++i){
//     for (IndexType j = 0; j < hbdlist.Nbond[i]; ++j){
//       indexMap[j] = j;
//       typeBuff[j] = bondType[hbdlist.bondIndex[j * hbdlist.stride + i]];
//       bkForceIndex[j] = hbdlist.bondIndex[j * hbdlist.stride + i];
//       bkData[j]       = hbdlist.data     [j * hbdlist.stride + i];
//     }
//     sortBuff (typeBuff, indexMap, hbdlist.Nbond[i]);
//     for (IndexType j = 0; j < hbdlist.Nbond[i]; ++j){
//       hbdlist.bondIndex[j * hbdlist.stride + i] = bkForceIndex[indexMap[j]];
//       hbdlist.data     [j * hbdlist.stride + i] = bkData      [indexMap[j]];
//     }
//   }
//   freeAPointer ((void**)&indexMap);
//   freeAPointer ((void**)&typeBuff);
//   freeAPointer ((void**)&bkForceIndex);
//   freeAPointer ((void**)&bkData);
// }    


void BondList::build()
{
  hbdlist.sort(bondType);
  buildDeviceBondList (hbdlist, dbdlist);
}

void initDeviceBondList (DeviceBondList & dbdlist)
{
  dbdlist.malloced = false;
  dbdlist.stride = 0;
  dbdlist.listLength = 0;
}

void destroyDeviceBondList(DeviceBondList &dbdlist )
{
  if (dbdlist.malloced) {
    cudaFree (dbdlist.data);
    cudaFree (dbdlist.bondIndex);
    cudaFree (dbdlist.Nbond);
    checkCUDAError ("destroyDeviceBondList");
  }
}

void buildDeviceBondList (HostBondList & hbdlist,
			  DeviceBondList & dbdlist)
{
  dbdlist.stride = hbdlist.stride;

  IndexType maxLength = 0;
  for (IndexType i = 0; i < hbdlist.stride; ++i){
    for (IndexType j = hbdlist.Nbond[i]; j < hbdlist.listLength_mem; ++j){
      hbdlist.data[j * hbdlist.stride + i] = MaxIndexValue;
      hbdlist.bondIndex[j * hbdlist.stride + i] = MaxForceIndexValue;
    }
    if (hbdlist.Nbond[i] > maxLength) maxLength = hbdlist.Nbond[i];
  }
  dbdlist.listLength = maxLength;
  
  cudaMalloc (&(dbdlist.data), 
	      sizeof(IndexType) * hbdlist.stride * maxLength);
  checkCUDAError ("buildDeviceBondList malloc data");
  cudaMalloc (&(dbdlist.bondIndex),
	      sizeof(TypeType) * hbdlist.stride * maxLength);
  checkCUDAError ("buildDeviceBondList malloc bondIndex");
  cudaMalloc (&(dbdlist.Nbond),
	      sizeof(IndexType) * hbdlist.stride);
  checkCUDAError ("buildDeviceBondList malloc Nbond");
  
  dbdlist.malloced = true;

  cudaMemcpy (dbdlist.data, hbdlist.data,
	      sizeof(IndexType) * hbdlist.stride * maxLength,
	      cudaMemcpyHostToDevice);
  checkCUDAError ("buildDeviceBondList cpy host data to device");
  cudaMemcpy (dbdlist.bondIndex, hbdlist.bondIndex,
	      sizeof(ForceIndexType) * hbdlist.stride * maxLength,
	      cudaMemcpyHostToDevice);
  checkCUDAError ("buildDeviceBondList cpy host bondIndex to device");
  cudaMemcpy (dbdlist.Nbond, hbdlist.Nbond,
	      sizeof(IndexType) * hbdlist.stride,
	      cudaMemcpyHostToDevice);
  checkCUDAError ("buildDeviceBondList cpy host Nbond to device");
}



HostBondList::HostBondList ()
{
  stride = 0;
  listLength_mem = 0;
  data = NULL;
  bondIndex = NULL;
  Nbond = NULL;
}

HostBondList::~HostBondList()
{
  freeAPointer ((void**)&data);
  freeAPointer ((void**)&bondIndex);
  freeAPointer ((void**)&Nbond);
}

void HostBondList::init (const IndexType & stride_)
{
  stride = stride_;
  listLength_mem = 1;
  data = (IndexType *) malloc (sizeof(IndexType) * stride * listLength_mem);
  if (data == NULL){
    throw MDExcptFailedMallocOnHost ("HostBondList::init", "data",
				     sizeof(IndexType) * stride * listLength_mem);
  }
  bondIndex = (ForceIndexType *) malloc
      (sizeof(ForceIndexType *) * stride * listLength_mem);
  if (bondIndex == NULL){
    throw MDExcptFailedMallocOnHost ("HostBondList::init", "bondIndex",
				     sizeof(ForceIndexType *) * stride * listLength_mem);
  }
  Nbond = (IndexType *) malloc (sizeof(IndexType) * stride);
  if (Nbond == NULL){
    throw MDExcptFailedMallocOnHost ("HostBondList::init", "Nbond",
				     sizeof(IndexType) * stride);
  }
  
  for (IndexType i = 0; i < stride; ++i){
    Nbond[i] = 0;
  }
  // for (IndexType i = 0; i < stride * listLength_mem; ++i){
  //   data[i] = MaxForceIndexValue;
  //   bondIndex[i] = 0;
  // }
}


void HostBondList::addBond (const IndexType & ii,
			    const IndexType & jj,
			    const ForceIndexType &looking)
{
  if (Nbond[ii] == listLength_mem || Nbond[jj] == listLength_mem){
    listLength_mem *= 2;
    data = (IndexType *) realloc (data, sizeof(IndexType) * stride * listLength_mem);
    if (data == NULL){
      throw MDExcptFailedReallocOnHost ("HostBondList::addBond", "data",
					sizeof(IndexType) * stride * listLength_mem);
    }
    bondIndex = (ForceIndexType *) realloc
	(bondIndex, sizeof(ForceIndexType *) * stride * listLength_mem);
    if (bondIndex == NULL){
      throw MDExcptFailedReallocOnHost ("HostBondList::init", "bondIndex",
					sizeof(ForceIndexType *) * stride * listLength_mem);
    }
  }
  data[Nbond[ii] * stride + ii] = jj;
  data[Nbond[jj] * stride + jj] = ii;
  bondIndex[Nbond[ii] * stride + ii] = looking;
  bondIndex[Nbond[jj] * stride + jj] = looking;
  Nbond[ii] ++;
  Nbond[jj] ++;
}


static void sortBuff (TypeType * ref, IndexType * indexMap, IndexType N)
{
  if (N == 0){
    return ;
  }
  for (IndexType i = 0; i < N - 1; ++i){
    IndexType j = i;
    while (j + 1 < N && ref[j] > ref[j+1]){
      TypeType tmptype = ref[j];
      ref[j] = ref[j+1];
      ref[j+1] = tmptype;
      IndexType tmpindex = indexMap[j];
      indexMap[j] = indexMap[j+1];
      indexMap[j+1] = tmpindex;
      j++;
    }
  }
}

void HostBondList::sort(mdBondInteraction_t * bondType)
{
  IndexType *indexMap = (IndexType *)malloc (sizeof(IndexType) * listLength_mem);
  if (indexMap == NULL){
    MDExcptFailedMallocOnHost ("BondList::sortBond", "indexMap",
			       sizeof(IndexType) * listLength_mem);
  }
  TypeType  *typeBuff = (TypeType *) malloc (sizeof(TypeType)  * listLength_mem);
  if (typeBuff == NULL){
    MDExcptFailedMallocOnHost ("BondList::sortBond", "typeBuff",
			       sizeof(TypeType)  * listLength_mem);
  }
  ForceIndexType * bkForceIndex = (ForceIndexType *) malloc (
      sizeof (ForceIndexType) * listLength_mem);
  if (bkForceIndex == NULL){
    MDExcptFailedMallocOnHost ("BondList::sortBond", "bkForceIndex",
			       sizeof (ForceIndexType) * listLength_mem);
  }			     
  IndexType * bkData = (IndexType *) malloc (
      sizeof (IndexType) * listLength_mem);
  if (bkData == NULL){
    MDExcptFailedMallocOnHost ("BondList::sortBond", "bkData",
			       sizeof (IndexType) * listLength_mem);
  }
  for (IndexType i = 0; i < stride; ++i){
    for (IndexType j = 0; j < Nbond[i]; ++j){
      indexMap[j] = j;
      typeBuff[j] = bondType[bondIndex[j * stride + i]];
      bkForceIndex[j] = bondIndex[j * stride + i];
      bkData[j]       = data     [j * stride + i];
    }
    sortBuff (typeBuff, indexMap, Nbond[i]);
    for (IndexType j = 0; j < Nbond[i]; ++j){
      bondIndex[j * stride + i] = bkForceIndex[indexMap[j]];
      data     [j * stride + i] = bkData      [indexMap[j]];
    }
  }
  freeAPointer ((void**)&indexMap);
  freeAPointer ((void**)&typeBuff);
  freeAPointer ((void**)&bkForceIndex);
  freeAPointer ((void**)&bkData);
}    
