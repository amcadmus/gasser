#ifndef __BondList_h_wanghan__
#define __BondList_h_wanghan__

#include "common.h"
#include "BondedInteraction.h"

struct HostBondList 
{
  IndexType stride; // is the expected larger than or equal to the number of Atoms
  // IndexType listLength;
  IndexType listLength_mem;
  IndexType * data;
  ForceIndexType * bondIndex;
  IndexType * Nbond;
public:
  HostBondList ();
  ~HostBondList ();
  void init (const IndexType & stride);
  void addBond (const IndexType &i, const IndexType &j,
		const ForceIndexType &fidx);
  void sort (InteractionType * bondType);
};


struct DeviceBondList 
{
  bool malloced ;
  IndexType stride; // is the expected larger than or equal to the number of Atoms
  IndexType listLength;
  IndexType * data;
  ForceIndexType * bondIndex;
  IndexType * Nbond;
};

void initDeviceBondList (DeviceBondList & dbdlist) ;
void buildDeviceBondList (HostBondList & hbdlist,
			  DeviceBondList & dbdlist);
void destroyDeviceBondList (DeviceBondList & dbdlist) ;

// get the jj th bond neighbor of atom ii
__device__ IndexType getBondAtomIndex (DeviceBondList dbdlist,
				       IndexType jj,
				       IndexType ii) 
{ return dbdlist.data[dbdlist.stride * jj + ii];}
__device__ ForceIndexType getBondForceIndex (DeviceBondList dbdlist,
					     IndexType jj,
					     IndexType ii)
{ return dbdlist.bondIndex[dbdlist.stride * jj + ii];}
// number of bonds of ii th atom
__device__ IndexType getNBond (DeviceBondList dbdlist,
			       IndexType ii)
{ return dbdlist.Nbond[ii]; }

#endif
