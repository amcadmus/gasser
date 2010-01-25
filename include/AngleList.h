#ifndef __AngleList_h_wanghan__
#define __AngleList_h_wanghan__

#include "common.h"
#include "AngleInteraction.h"

struct HostAngleList 
{
  IndexType stride; // is the expected larger than or equal to the number of Atoms
  // IndexType listLength;
  IndexType listLength_mem;
  IndexType * angleNei;
  IndexType * myPosi;
  ForceIndexType * angleIndex;
  IndexType * Nangle;
public:
  HostAngleList ();
  ~HostAngleList ();
  void init (const IndexType & stride);
  void addAngle (const IndexType &i,
		 const IndexType &j,
		 const IndexType &k,
		 const ForceIndexType &fidx);
  void sort (InteractionType * angleType);
};


struct DeviceAngleList 
{
  bool malloced ;
  IndexType stride; // is the expected larger than or equal to the number of Atoms
  IndexType listLength;
  IndexType * angleNei;
  IndexType * myPosi;
  ForceIndexType * angleIndex;
  IndexType * Nangle;
};

void initDeviceAngleList (DeviceAngleList & danglelist) ;
void buildDeviceAngleList (HostAngleList & hanglelist,
			  DeviceAngleList & danglelist);
void destroyDeviceAngleList (DeviceAngleList & danglelist) ;

// // get the jj th angle neighbor of atom ii
// __device__ IndexType getAngleAtomIndex (DeviceAngleList danglelist,
// 					IndexType jj,
// 					IndexType ii) 
// { return danglelist.data[danglelist.stride * jj + ii];}
// __device__ ForceIndexType getAngleForceIndex (DeviceAngleList danglelist,
// 					     IndexType jj,
// 					     IndexType ii)
// { return danglelist.angleIndex[danglelist.stride * jj + ii];}
// // number of angles of ii th atom
// __device__ IndexType getNAngle (DeviceAngleList danglelist,
// 			       IndexType ii)
// { return danglelist.Nangle[ii]; }

#endif

