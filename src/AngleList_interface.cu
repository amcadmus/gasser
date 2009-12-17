// #include "AngleList_interface.h"

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




