#ifndef __BondedInteractionList_h_wanghan__
#define __BondedInteractionList_h_wanghan__

#include "common.h"
#include "BondList.h"
#include "AngleList.h"
#include "SystemBondedInteraction.h"
#include "Topology.h"
#include "MDSystem_interface.h"
#include "MDException.h"
#include "Reshuffle_interface.h"

class BondedInteractionList : public Reshufflable
{
  HostBondList    hbondlist;
  HostAngleList   hanglelist;
public:
  DeviceBondList  dbondlist;
  DeviceAngleList danglelist;
private: // reshuffle backup things
  DeviceBondList  bkdbondlist;
  DeviceAngleList bkdanglelist;
public:
  BondedInteractionList (); 
  ~BondedInteractionList ();
  void reinit (const MDSystem & sysData,
	       const Topology::System & sysTop,
	       const SystemBondedInteraction & sysBdInter);
public:
  // const DeviceBondList & deviceBondList () const
  //     {return dbondlist;}
  // const DeviceAngleList & deviceAngleList () const
  //     {return danglelist;}
  virtual void reshuffle (const IndexType * indexTable,
			  const IndexType & numAtom,
			  MDTimer * timer = NULL) ;
};


#endif
