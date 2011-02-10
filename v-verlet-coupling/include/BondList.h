#ifndef __BondList_h_wanghan__
#define __BondList_h_wanghan__

#include "common.h"

struct HostBondList 
{
  IndexType stride; // is the expected larger than or equal to the number of Atoms
  // IndexType listLength;
  IndexType maxNumBond;
  IndexType * bondNeighborIndex;
  IndexType * bondIndex;
  IndexType * numBond;
private:
  void clearMem ();
public:
  HostBondList ();
  ~HostBondList ();
  void reinit (const IndexType & stride,
	       const IndexType & maxNumBond);
  void clearBond ();
  void addBond (const IndexType &i,
		const IndexType &j,
		const IndexType &bondIndex);
};


struct DeviceBondList 
{
  bool malloced ;
  IndexType stride; // is the expected larger than or equal to the number of Atoms
  IndexType maxNumBond;
  IndexType * bondNeighborIndex;
  IndexType * bondIndex;
  IndexType * numBond;
};

void initDeviceBondList (DeviceBondList & dbdlist) ;
void mallocDeviceBondList (const HostBondList & hbdlist,
			   DeviceBondList & dbdlist);
void copyDeviceBondList (const HostBondList & hbdlist,
			 DeviceBondList & dbdlist);
void copyDeviceBondList (const DeviceBondList & dbdlist1,
			 DeviceBondList & dbdlist);
void destroyDeviceBondList (DeviceBondList & dbdlist) ;

#endif
