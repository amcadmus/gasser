#ifndef __BondList_interface_h_wanghan__
#define __BondList_interface_h_wanghan__

#include "common.h"
#include "BondList.h"
#include "MDSystem.h"
#include "BondedInteraction.h"

class BondList
{
  HostBondList hbdlist;
  void sortBond();
  IndexType NBondForce_mem;
  IndexType paramLength_mem;
public:
  DeviceBondList dbdlist;
  mdBondInteraction_t * bondType;
  IndexType * paramPosi;
  IndexType NBondForce;
  ScalorType * param;
  IndexType paramLength;
public:
  BondList () ;
  ~BondList ();
  BondList(const DeviceMDData & ddata, const IndexType & listLength) 
      {init (ddata, listLength);}
  void init (const DeviceMDData & ddata, const IndexType & listLength);
  void addBond (const IndexType & ii, const IndexType & jj,
		const mdBondInteraction_t & type,
		const ScalorType * param);
  void build ();
};

	     

#endif
