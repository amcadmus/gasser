#ifndef __ExclusionList_h_wanghan__
#define __ExclusionList_h_wanghan__

#include "common.h"
#include "SystemNonBondedInteraction.h"
#include "Topology.h"
#include "MDSystem_interface.h"
#include "MDException.h"
#include "Reshuffle_interface.h"

struct HostExclusionList 
{
  IndexType stride;
  IndexType maxNumExclusion;
  IndexType * exclusionNeighborIndex;
  IndexType * numExclusion;
private:
  void clearMem ();
public:
  HostExclusionList ();
  ~HostExclusionList ();
  void reinit (const IndexType & stride,
	       const IndexType & maxNumExclusion);
  void clearExclusion ();
  void addExclusion (const IndexType & i,
		     const IndexType & j);
};

struct DeviceExclusionList 
{
  bool malloced;
  IndexType stride;
  IndexType maxNumExclusion;
  IndexType * exclusionNeighborIndex;
  IndexType * numExclusion;
};

void initDeviceExclusionList   (DeviceExclusionList & dexcllist);
void mallocDeviceExclusionList (const HostExclusionList & hexcllist,
				DeviceExclusionList & dexcllist);
void copyDeviceExclusionList   (const HostExclusionList & hexcllist,
				DeviceExclusionList & dexcllist);
void copyDeviceExclusionList   (const DeviceExclusionList & dexcllist1,
				DeviceExclusionList & dexcllist);
void destroyDeviceBondList     (DeviceExclusionList & dexcllist);


class ExclusionList : public Reshufflable
{
private:
  DeviceExclusionList bkdexcllist;
  HostExclusionList   hexcllist;
public:
  DeviceExclusionList dexcllist;
public:
  ExclusionList ();
  ~ExclusionList ();
  void reinit (const MDSystem & sysData,
	       const Topology::System & sysTop,
	       const SystemNonBondedInteraction & sysNbInter);
  virtual void reshuffle (const IndexType * indexTable,
			  const IndexType & numAtom,
			  MDTimer * timer = NULL) ;
};


#endif

