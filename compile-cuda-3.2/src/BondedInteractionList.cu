#define DEVICE_CODE

#include "BondedInteractionList.h"

BondedInteractionList::
BondedInteractionList ()
{
  initDeviceBondList (dbondlist);
  initDeviceAngleList (danglelist);
  initDeviceBondList (bkdbondlist);
  initDeviceAngleList (bkdanglelist);
}

BondedInteractionList::
~BondedInteractionList ()
{
  destroyDeviceBondList (dbondlist);
  destroyDeviceAngleList (danglelist);
}


void BondedInteractionList::
reinit (const MDSystem & sysData,
	const Topology::System & sysTop,
	const SystemBondedInteraction & sysBdInter)
{
  if (sysData.hdata.numAtom != sysTop.indexShift.back()){
    throw MDExcptWrongNumberAtomDataTopology ();
  }

  IndexType maxNumBond = 0;
  for (unsigned i = 0; i < sysBdInter.bondIndex.size(); ++i){
    for (unsigned j = 0; j < sysBdInter.bondIndex[i].size(); ++j){
      IndexType c ;
      if ((c=sysBdInter.bondIndex[i][j].size()) > maxNumBond){
	maxNumBond = c;
      }
    }
  }
  IndexType maxNumAngle = 0;
  for (unsigned i = 0; i < sysBdInter.angleIndex.size(); ++i){
    for (unsigned j = 0; j < sysBdInter.angleIndex[i].size(); ++j){
      IndexType c ;
      if ((c=sysBdInter.angleIndex[i][j].size()) > maxNumAngle){
	maxNumAngle = c;
      }
    }
  }
  
  hbondlist .reinit (sysData.hdata.numAtom, maxNumBond);
  hanglelist.reinit (sysData.hdata.numAtom, maxNumAngle);
  
  {
    IndexType shift0 = 0;
    for (unsigned i = 0; i < sysTop.molecules.size(); ++i){
      shift0 = sysTop.indexShift[i];
      IndexType molSize = sysTop.molecules[i].size();
      for (unsigned j = 0; j < sysTop.numbers[i]; ++j){
	IndexType shift1 = j * molSize;
	IndexType indexSift = shift0 + shift1;
	for (unsigned k = 0; k < molSize; ++k){
	  for (unsigned l = 0; l < sysBdInter.bondIndex[i][k].size(); ++l){
	    hbondlist.addBond (indexSift + k,
			       indexSift + sysBdInter.bondNeighborIndex[i][k][l],
			       sysBdInter.bondIndex[i][k][l]);
	  }
	}
      }
    }
  }
  {
    IndexType shift0 = 0;
    for (unsigned i = 0; i < sysTop.molecules.size(); ++i){
      shift0 = sysTop.indexShift[i];
      IndexType molSize = sysTop.molecules[i].size();
      for (unsigned j = 0; j < sysTop.numbers[i]; ++j){
	IndexType shift1 = j * molSize;
	IndexType indexSift = shift0 + shift1;
	for (unsigned k = 0; k < molSize; ++k){
	  for (unsigned l = 0; l < sysBdInter.angleIndex[i][k].size(); ++l){
	    hanglelist.addAngle (
		indexSift + k,
		indexSift + sysBdInter.angleNeighborIndex[i][k][2*l],
		indexSift + sysBdInter.angleNeighborIndex[i][k][2*l+1],
		sysBdInter.angleIndex[i][k][l],
		sysBdInter.anglePosi [i][k][l]);
	  }
	}
      }
    }
  }
  
  destroyDeviceBondList (dbondlist);
  destroyDeviceBondList (bkdbondlist);
  initDeviceBondList (dbondlist);
  initDeviceBondList (bkdbondlist);
  mallocDeviceBondList (hbondlist, dbondlist);
  mallocDeviceBondList (hbondlist, bkdbondlist);
  copyDeviceBondList (hbondlist, dbondlist);

  destroyDeviceAngleList (danglelist);
  destroyDeviceAngleList (bkdanglelist);
  initDeviceAngleList (danglelist);
  initDeviceAngleList (bkdanglelist);
  mallocDeviceAngleList (hanglelist, danglelist);
  mallocDeviceAngleList (hanglelist, bkdanglelist);
  copyDeviceAngleList (hanglelist, danglelist);
}


static __global__ void
Reshuffle_reshuffleBondList (const IndexType numAtom,
			     const IndexType * bdlistData2,
			     const IndexType * bdlistBondIndex2,
			     const IndexType * bdlistNumB2,
			     const IndexType stride,
			     const IndexType listLength,
			     const IndexType * idxTable,
			     IndexType * bdlistData,
			     IndexType * bdlistBondIndex,
			     IndexType * bdlistNumB)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;

  if (ii < numAtom){
    IndexType toid = idxTable[ii];
    bdlistNumB[toid] = bdlistNumB2[ii];
    for (IndexType jj = 0; jj < bdlistNumB2[ii]; ++jj){
      bdlistData    [jj * stride + toid] = 
	  idxTable[bdlistData2[jj * stride + ii]];
      bdlistBondIndex[jj * stride + toid] = 
	  bdlistBondIndex2[jj * stride + ii];
    }
  }
}

static __global__ void
Reshuffle_reshuffleAngleList (const IndexType numAtom,
			      const IndexType * bkAngleListNei,
			      const IndexType * bkAngleListPosi,
			      const IndexType * bkAngleListAngleIndex,
			      const IndexType * bkAngleListNangle,
			      const IndexType stride,
			      const IndexType listLength,
			      const IndexType * idxTable,
			      IndexType * angleListNei,
			      IndexType * angleListPosi,
			      IndexType * angleListAngleIndex,
			      IndexType * angleListNangle)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;

  if (ii < numAtom){
    IndexType toid = idxTable[ii];
    angleListNangle[toid] = bkAngleListNangle[ii];
    for (IndexType jj = 0; jj < bkAngleListNangle[ii]; ++ jj){
      angleListNei[((jj<<1)+0)*stride + toid] =
	  idxTable[bkAngleListNei[((jj<<1)+0)*stride + ii]];
      angleListNei[((jj<<1)+1)*stride + toid] =
	  idxTable[bkAngleListNei[((jj<<1)+1)*stride + ii]];
      angleListPosi[jj * stride + toid] = bkAngleListPosi[jj * stride + ii];
      angleListAngleIndex[jj * stride + toid] =
	  bkAngleListAngleIndex[jj * stride + ii];
    }
  }
}

void BondedInteractionList::
reshuffle (const IndexType * indexTable,
	   const IndexType & numAtom,
	   MDTimer * timer) 
{
  if (timer != NULL) timer->tic(mdTimeReshuffleSystem);

  copyDeviceBondList (dbondlist, bkdbondlist);
  copyDeviceAngleList (danglelist, bkdanglelist);

  dim3 myBlockDim, atomGridDim;
  myBlockDim.y = 1;
  myBlockDim.z = 1;
  myBlockDim.x = DefaultNThreadPerBlock;
  IndexType nob;
  if (numAtom % myBlockDim.x == 0){
    nob = numAtom / myBlockDim.x;
  } else {
    nob = numAtom / myBlockDim.x + 1;
  }
  atomGridDim = toGridDim (nob);
  
  Reshuffle_reshuffleBondList
      <<<atomGridDim,myBlockDim>>> (
  	  numAtom,
  	  bkdbondlist.bondNeighborIndex, 
	  bkdbondlist.bondIndex, 
	  bkdbondlist.numBond,
  	  dbondlist.stride,
	  dbondlist.maxNumBond, 
  	  indexTable,
  	  dbondlist.bondNeighborIndex, 
	  dbondlist.bondIndex, 
	  dbondlist.numBond);
  checkCUDAError ("BondedInteractionList::reshuffle, reshuffle bond list");
  Reshuffle_reshuffleAngleList
      <<<atomGridDim,myBlockDim>>> (
	  numAtom,
	  bkdanglelist.angleNeighborIndex,
	  bkdanglelist.anglePosi,
	  bkdanglelist.angleIndex,
	  bkdanglelist.numAngle,
  	  danglelist.stride,
	  danglelist.maxNumAngle, 
  	  indexTable,
	  danglelist.angleNeighborIndex,
	  danglelist.anglePosi,
	  danglelist.angleIndex,
	  danglelist.numAngle);
  checkCUDAError ("BondedInteractionList::reshuffle, reshuffle angle list");

  if (timer != NULL) timer->toc(mdTimeReshuffleSystem);
}

