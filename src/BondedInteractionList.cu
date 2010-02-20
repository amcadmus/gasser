#include "BondedInteractionList.h"

BondedInteractionList::
BondedInteractionList ()
{
  initDeviceBondList (dbondlist);
  initDeviceAngleList (danglelist);
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
      IndexType molSize = sysTop.molecules[i].size();
      for (unsigned j = 0; j < sysTop.numbers[i]; ++j){
	IndexType shift1 = j * molSize;
	for (unsigned k = 0; k < molSize; ++k){
	  IndexType indexSift = shift0 + shift1;
	  for (unsigned l = 0; l < sysBdInter.bondIndex[i][k].size(); ++l){
	    hbondlist.addBond (indexSift + k,
			       indexSift + sysBdInter.bondNeighborIndex[i][k][l],
			       sysBdInter.bondIndex[i][k][l]);
	  }
	}
      }
      shift0 += sysTop.indexShift[i];
    }
  }
  {
    IndexType shift0 = 0;
    for (unsigned i = 0; i < sysTop.molecules.size(); ++i){
      IndexType molSize = sysTop.molecules[i].size();
      for (unsigned j = 0; j < sysTop.numbers[i]; ++j){
	IndexType shift1 = j * molSize;
	for (unsigned k = 0; k < molSize; ++k){
	  IndexType indexSift = shift0 + shift1;
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
      shift0 += sysTop.indexShift[i];
    }
  }
  
  destroyDeviceBondList (dbondlist);
  initDeviceBondList (dbondlist);
  mallocDeviceBondList (hbondlist, dbondlist);
  copyDeviceBondList (hbondlist, dbondlist);
  destroyDeviceAngleList (danglelist);
  initDeviceAngleList (danglelist);
  mallocDeviceAngleList (hanglelist, danglelist);
  copyDeviceAngleList (hanglelist, danglelist);
}

	
