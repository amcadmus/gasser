#ifndef __BondedInteractionList_h_wanghan__
#define __BondedInteractionList_h_wanghan__

#include "common.h"
#include "BondList.h"
#include "AngleList.h"
#include "SystemBondedInteraction.h"
#include "Topology.h"
#include "MDSystem_interface.h"
#include "MDException.h"


class BondedInteractionList
{
  HostBondList   hbondlist;
  DeviceBondList dbondlist;
  HostAngleList   hanglelist;
  DeviceAngleList danglelist;
public:
  BondedInteractionList (); 
  ~BondedInteractionList ();
  void reinit (const MDSystem & sysData,
	       const Topology::System & sysTop,
	       const SystemBondedInteraction & sysBdInter);
public:
  const DeviceBondList & deviceBondList () const
      {return dbondlist;}
  const DeviceAngleList & deviceAngleList () const
      {return danglelist;}
};


#endif
