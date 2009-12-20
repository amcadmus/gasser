#include "AngleList_interface.h"

AngleList::AngleList ()
{
  angleType = NULL;
  paramPosi = NULL;
  NAngleForce = 0;
  param = NULL;
  paramLength = 0;
  initDeviceAngleList (danglelist);
}

AngleList::~AngleList()
{
  freeAPointer ((void**)&angleType);
  freeAPointer ((void**)&paramPosi);
  freeAPointer ((void**)&param);
  destroyDeviceAngleList (danglelist);
}


void AngleList::init (const DeviceMDData & ddata)
{
  hanglelist.init (ddata.numAtom);

  NAngleForce_mem = 1024;
  paramLength_mem = NAngleForce_mem * 3;
  angleType = (mdAngleInteraction_t *) realloc (
      angleType, sizeof(mdAngleInteraction_t) * NAngleForce_mem);
  if (angleType == NULL) {
    throw MDExcptFailedReallocOnHost ("AngleList::init", "angleType",
				      sizeof(mdAngleInteraction_t) * NAngleForce_mem);
  }
  paramPosi = (IndexType *) realloc (
      paramPosi, sizeof(mdAngleInteraction_t) * NAngleForce_mem);
  if (paramPosi == NULL){
    throw MDExcptFailedReallocOnHost ("AngleList::init", "paramPosi",
				      sizeof(mdAngleInteraction_t) * NAngleForce_mem);
  }				      
  param = (ScalorType *) realloc (
      param, sizeof(ScalorType) * paramLength_mem);
  if (param == NULL){
    throw MDExcptFailedReallocOnHost ("AngleList::init", "param",
				      sizeof(ScalorType) * paramLength_mem);
  }
}

  
void AngleList::addAngle (const IndexType & ii,
			  const IndexType & jj,
			  const IndexType & kk,
			  const mdAngleInteraction_t & type,
			  const ScalorType * thisparam)
{
  bool exist = false;
  ForceIndexType looking;;
  IndexType NParam = calNumAngleParameter (type);
  for (looking = 0; looking < NAngleForce; ++looking){
    if (type == angleType[looking]){
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
    if (NAngleForce == NAngleForce_mem){
      NAngleForce_mem *= 2;
      angleType = (mdAngleInteraction_t *) realloc (
	  angleType, sizeof(mdAngleInteraction_t) * NAngleForce_mem);
      if (angleType == NULL) {
	throw MDExcptFailedReallocOnHost ("AngleList::init", "angleType",
					  sizeof(mdAngleInteraction_t) * NAngleForce_mem);
      }
      paramPosi = (IndexType *) realloc (
	  paramPosi, sizeof(mdAngleInteraction_t) * NAngleForce_mem);
      if (paramPosi == NULL){
	throw MDExcptFailedReallocOnHost ("AngleList::init", "paramPosi",
					  sizeof(mdAngleInteraction_t) * NAngleForce_mem);
      }				      
    }
    if (paramLength == paramLength_mem){
      paramLength_mem *= 2;
      param = (ScalorType *) realloc (
	  param, sizeof(ScalorType) * paramLength_mem);
      if (param == NULL){
	throw MDExcptFailedReallocOnHost ("AngleList::init", "param",
					  sizeof(ScalorType) * paramLength_mem);
      }
    } 
    angleType [NAngleForce] = type;
    paramPosi[NAngleForce] = paramLength;
    for (IndexType i = 0; i < NParam; ++i){
      param[paramPosi[NAngleForce] + i] = thisparam[i];
    }
    NAngleForce ++;
    paramLength += NParam;
  }
  hanglelist.addAngle (ii, jj, kk, looking);
}

void AngleList::build()
{
  hanglelist.sort(angleType);
  buildDeviceAngleList (hanglelist, danglelist);
}


void initDeviceAngleList (DeviceAngleList & danglelist)
{
  danglelist.malloced = false;
  danglelist.stride = 0;
  danglelist.listLength = 0;
}

void destroyDeviceAngleList(DeviceAngleList &danglelist )
{
  if (danglelist.malloced) {
    cudaFree (danglelist.angleNei);
    cudaFree (danglelist.angleIndex);
    cudaFree (danglelist.myPosi);
    cudaFree (danglelist.Nangle);
    checkCUDAError ("destroyDeviceAngleList");
  }
}

void buildDeviceAngleList (HostAngleList & hanglelist,
			  DeviceAngleList & danglelist)
{
  danglelist.stride = hanglelist.stride;

  IndexType maxLength = 0;
  for (IndexType i = 0; i < hanglelist.stride; ++i){
    for (IndexType j = hanglelist.Nangle[i]; j < hanglelist.listLength_mem; ++j){
      hanglelist.angleNei[((j<<1))   * hanglelist.stride + i] = MaxIndexValue;
      hanglelist.angleNei[((j<<1)+1) * hanglelist.stride + i] = MaxIndexValue;
      hanglelist.myPosi    [j * hanglelist.stride + i] = 0;
      hanglelist.angleIndex[j * hanglelist.stride + i] = MaxForceIndexValue;
    }
    if (hanglelist.Nangle[i] > maxLength) maxLength = hanglelist.Nangle[i];
  }
  danglelist.listLength = maxLength;
  
  cudaMalloc (&(danglelist.angleNei), 
	      sizeof(IndexType) * hanglelist.stride * 2 * maxLength);
  checkCUDAError ("buildDeviceAngleList malloc angleNei");
  cudaMalloc (&(danglelist.myPosi), 
	      sizeof(IndexType) * hanglelist.stride * maxLength);
  checkCUDAError ("buildDeviceAngleList malloc myPosi");
  cudaMalloc (&(danglelist.angleIndex),
	      sizeof(TypeType) * hanglelist.stride * maxLength);
  checkCUDAError ("buildDeviceAngleList malloc angleIndex");
  cudaMalloc (&(danglelist.Nangle),
	      sizeof(IndexType) * hanglelist.stride);
  checkCUDAError ("buildDeviceAngleList malloc Nangle");
  
  danglelist.malloced = true;

  cudaMemcpy (danglelist.angleNei, hanglelist.angleNei,
	      sizeof(IndexType) * hanglelist.stride * 2 * maxLength,
	      cudaMemcpyHostToDevice);
  checkCUDAError ("buildDeviceAngleList cpy host angleNei to device");
  cudaMemcpy (danglelist.myPosi, hanglelist.myPosi,
	      sizeof(ForceIndexType) * hanglelist.stride * maxLength,
	      cudaMemcpyHostToDevice);
  checkCUDAError ("buildDeviceAngleList cpy host myPosi to device");
  cudaMemcpy (danglelist.angleIndex, hanglelist.angleIndex,
	      sizeof(ForceIndexType) * hanglelist.stride * maxLength,
	      cudaMemcpyHostToDevice);
  checkCUDAError ("buildDeviceAngleList cpy host angleIndex to device");
  cudaMemcpy (danglelist.Nangle, hanglelist.Nangle,
	      sizeof(IndexType) * hanglelist.stride,
	      cudaMemcpyHostToDevice);
  checkCUDAError ("buildDeviceAngleList cpy host Nangle to device");
}


HostAngleList::HostAngleList ()
{
  stride = 0;
  listLength_mem = 0;
  angleNei = NULL;
  myPosi = NULL;
  angleIndex = NULL;
  Nangle = NULL;
}

HostAngleList::~HostAngleList()
{
  freeAPointer ((void**)&angleNei);
  freeAPointer ((void**)&myPosi);
  freeAPointer ((void**)&angleIndex);
  freeAPointer ((void**)&Nangle);
}

void HostAngleList::init (const IndexType & stride_)
{
  stride = stride_;
  listLength_mem = 1;
  angleNei = (IndexType *) malloc (sizeof(IndexType) * stride * 2 * listLength_mem);
  if (angleNei == NULL){
    throw MDExcptFailedMallocOnHost ("HostAngleList::init", "angleNei",
				     sizeof(IndexType) * stride * 2 * listLength_mem);
  }
  myPosi = (IndexType *) malloc (sizeof(IndexType) * stride * listLength_mem);
  if (myPosi == NULL){
    throw MDExcptFailedMallocOnHost ("HostAngleList::init", "myPosi",
				     sizeof(IndexType) * stride * listLength_mem);
  }
  angleIndex = (ForceIndexType *) malloc
      (sizeof(ForceIndexType *) * stride * listLength_mem);
  if (angleIndex == NULL){
    throw MDExcptFailedMallocOnHost ("HostAngleList::init", "angleIndex",
				     sizeof(ForceIndexType *) * stride * listLength_mem);
  }
  Nangle = (IndexType *) malloc (sizeof(IndexType) * stride);
  if (Nangle == NULL){
    throw MDExcptFailedMallocOnHost ("HostAngleList::init", "Nangle",
				     sizeof(IndexType) * stride);
  }
  
  for (IndexType i = 0; i < stride; ++i){
    Nangle[i] = 0;
  }
}


void HostAngleList::addAngle (const IndexType & ii,
			      const IndexType & jj,
			      const IndexType & kk,
			      const ForceIndexType &looking)
{
  if (Nangle[ii] == listLength_mem ||
      Nangle[jj] == listLength_mem ||
      Nangle[kk] == listLength_mem  ){
    listLength_mem *= 2;
    angleNei = (IndexType *) realloc
	(angleNei, sizeof(IndexType) * stride * 2 * listLength_mem);
    if (angleNei == NULL){
      throw MDExcptFailedReallocOnHost ("HostAngleList::addAngle", "angleNei",
					sizeof(IndexType) * stride * 2 * listLength_mem);
    }
    myPosi = (IndexType *) realloc
	(myPosi, sizeof(IndexType) * stride * listLength_mem);
    if (myPosi == NULL){
      throw MDExcptFailedReallocOnHost ("HostAngleList::addAngle", "myPosi",
					sizeof(IndexType) * stride * listLength_mem);
    }
    angleIndex = (ForceIndexType *) realloc
	(angleIndex, sizeof(ForceIndexType *) * stride * listLength_mem);
    if (angleIndex == NULL){
      throw MDExcptFailedReallocOnHost ("HostAngleList::init", "angleIndex",
					sizeof(ForceIndexType *) * stride * listLength_mem);
    }
  }

  angleNei[((Nangle[ii] << 1))     * stride + ii] = jj;
  angleNei[((Nangle[ii] << 1) + 1) * stride + ii] = kk;
  angleNei[((Nangle[jj] << 1))     * stride + jj] = ii;
  angleNei[((Nangle[jj] << 1) + 1) * stride + jj] = kk;
  angleNei[((Nangle[kk] << 1))     * stride + kk] = jj;
  angleNei[((Nangle[kk] << 1) + 1) * stride + kk] = ii;
  myPosi[Nangle[ii] * stride + ii] = 0;
  myPosi[Nangle[jj] * stride + jj] = 1;
  myPosi[Nangle[kk] * stride + kk] = 2;
  angleIndex[Nangle[ii] * stride + ii] = looking;
  angleIndex[Nangle[jj] * stride + jj] = looking;
  angleIndex[Nangle[kk] * stride + kk] = looking;
  Nangle[ii] ++;
  Nangle[jj] ++;
  Nangle[kk] ++;
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

void HostAngleList::sort(mdAngleInteraction_t * angleType)
{
  IndexType *indexMap = (IndexType *)malloc (sizeof(IndexType) * listLength_mem);
  if (indexMap == NULL){
    MDExcptFailedMallocOnHost ("AngleList::sortAngle", "indexMap",
			       sizeof(IndexType) * listLength_mem);
  }
  TypeType  *typeBuff = (TypeType *) malloc (sizeof(TypeType)  * listLength_mem);
  if (typeBuff == NULL){
    MDExcptFailedMallocOnHost ("AngleList::sortAngle", "typeBuff",
			       sizeof(TypeType)  * listLength_mem);
  }
  ForceIndexType * bkForceIndex = (ForceIndexType *) malloc (
      sizeof (ForceIndexType) * listLength_mem);
  if (bkForceIndex == NULL){
    MDExcptFailedMallocOnHost ("AngleList::sortAngle", "bkForceIndex",
			       sizeof (ForceIndexType) * listLength_mem);
  }			     
  IndexType * bkAngleNei = (IndexType *) malloc (
      sizeof (IndexType) * 2 * listLength_mem);
  if (bkAngleNei == NULL){
    MDExcptFailedMallocOnHost ("AngleList::sortAngle", "bkAngleNei",
			       sizeof (IndexType) * 2 * listLength_mem);
  }
  IndexType * bkMyPosi = (IndexType *) malloc (
      sizeof (IndexType) * listLength_mem);
  if (bkMyPosi == NULL){
    MDExcptFailedMallocOnHost ("AngleList::sortAngle", "bkMyPosi",
			       sizeof (IndexType) * listLength_mem);
  }
  
  for (IndexType i = 0; i < stride; ++i){
    for (IndexType j = 0; j < Nangle[i]; ++j){
      indexMap[j] = j;
      typeBuff[j] = angleType[angleIndex[j * stride + i]];
      bkForceIndex[j] = angleIndex[j * stride + i];
      bkMyPosi[j] = myPosi[j * stride + i];
      bkAngleNei[(j<<1)]       = angleNei[(j<<1)       * stride + i];
      bkAngleNei[(j<<1) + 1]   = angleNei[((j<<1) + 1) * stride + i];
    }
    sortBuff (typeBuff, indexMap, Nangle[i]);
    for (IndexType j = 0; j < Nangle[i]; ++j){
      angleIndex[j * stride + i] = bkForceIndex[indexMap[j]];
      myPosi[j * stride + i] = bkMyPosi[indexMap[j]];
      angleNei[(j<<1)     * stride + i] = bkAngleNei[(indexMap[j]<<1)];
      angleNei[((j<<1)+1) * stride + i] = bkAngleNei[(indexMap[j]<<1)+1];
    }
  }
  freeAPointer ((void**)&indexMap);
  freeAPointer ((void**)&typeBuff);
  freeAPointer ((void**)&bkForceIndex);
  freeAPointer ((void**)&bkAngleNei);
  freeAPointer ((void**)&bkMyPosi);
}    



// AngleList::AngleList ()
// {
//   angleType = NULL;
//   paramPosi = NULL;
//   NAngleForce = 0;
//   param = NULL;
//   paramLength = 0;
//   initDeviceAngleList (danglelist);
// }

// AngleList::~AngleList()
// {
//   freeAPointer ((void**)&angleType);
//   freeAPointer ((void**)&paramPosi);
//   freeAPointer ((void**)&param);
//   destroyDeviceAngleList (danglelist);
// }

// void AngleList::init (const DeviceMDData & ddata,
// 		      const IndexType & listLength)
// {
//   hanglelist.init (ddata.numAtom, listLength);

//   NAngleForce_mem = 1024;
//   paramLength_mem = NAngleForce_mem * 3;
//   angleType = (mdAngleInteraction_t *) realloc (
//       angleType, sizeof(mdAngleInteraction_t) * NAngleForce_mem);
//   if (angleType == NULL) {
//     throw MDExcptFailedReallocOnHost ("AngleList::init", "angleType",
// 				      sizeof(mdAngleInteraction_t) * NAngleForce_mem);
//   }
//   paramPosi = (IndexType *) realloc (
//       paramPosi, sizeof(mdAngleInteraction_t) * NAngleForce_mem);
//   if (paramPosi == NULL){
//     throw MDExcptFailedReallocOnHost ("AngleList::init", "paramPosi",
// 				      sizeof(mdAngleInteraction_t) * NAngleForce_mem);
//   }				      
//   param = (ScalorType *) realloc (
//       param, sizeof(ScalorType) * paramLength_mem);
//   if (param == NULL){
//     throw MDExcptFailedReallocOnHost ("AngleList::init", "param",
// 				      sizeof(ScalorType) * paramLength_mem);
//   }
// }

// void AngleList::addAngle (const IndexType & ii,
// 			  const IndexType & jj,
// 			  const IndexType & kk,
// 			  const mdAngleInteraction_t & type,
// 			  const ScalorType * thisparam)
// {
//   bool exist = false;
//   ForceIndexType looking;;
//   IndexType NParam = calNumAngleParameter (type);
//   for (looking = 0; looking < NAngleForce; ++looking){
//     if (type == angleType[looking]){
//       bool same = true;
//       for (IndexType i = 0; i < NParam; ++i){
// 	if (thisparam[i] != param[paramPosi[looking] + i]){
// 	  same = false;
// 	  break;
// 	}
//       }
//       if (same){
// 	exist = true;
// 	break;
//       }
//     }
//   }
//   if (!exist){
//     if (NAngleForce == NAngleForce_mem){
//       NAngleForce_mem *= 2;
//       angleType = (mdAngleInteraction_t *) realloc (
// 	  angleType, sizeof(mdAngleInteraction_t) * NAngleForce_mem);
//       if (angleType == NULL) {
// 	throw MDExcptFailedReallocOnHost ("AngleList::init", "angleType",
// 					  sizeof(mdAngleInteraction_t) * NAngleForce_mem);
//       }
//       paramPosi = (IndexType *) realloc (
// 	  paramPosi, sizeof(mdAngleInteraction_t) * NAngleForce_mem);
//       if (paramPosi == NULL){
// 	throw MDExcptFailedReallocOnHost ("AngleList::init", "paramPosi",
// 					  sizeof(mdAngleInteraction_t) * NAngleForce_mem);
//       }				      
//     }
//     if (paramLength == paramLength_mem){
//       paramLength_mem *= 2;
//       param = (ScalorType *) realloc (
// 	  param, sizeof(ScalorType) * paramLength_mem);
//       if (param == NULL){
// 	throw MDExcptFailedReallocOnHost ("AngleList::init", "param",
// 					  sizeof(ScalorType) * paramLength_mem);
//       }
//     } 
//     angleType [NAngleForce] = type;
//     paramPosi[NAngleForce] = paramLength;
//     for (IndexType i = 0; i < NParam; ++i){
//       param[paramPosi[NAngleForce] + i] = thisparam[i];
//     }
//     NAngleForce ++;
//     paramLength += NParam;
//   }
//   hanglelist.addAngle (ii, jj, kk, looking);
// }




