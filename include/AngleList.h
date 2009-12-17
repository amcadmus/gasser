#ifndef __AngleList_h_wanghan__
#define __AngleList_h_wanghan__

#include "common.h"

struct HostAngleList 
{
  IndexType stride; // is the expected larger than or equal to the number of Atoms
  IndexType listLength;
  IndexType * angleNei;
  IndexType * myPosi;
  ForceIndexType * angleIndex;
  IndexType * Nangle;
public:
  HostAngleList ();
  ~HostAngleList ();
  void init (const IndexType & stride,
	     const IndexType & listLength);
  void addAngle (const IndexType &i,
		 const IndexType &j,
		 const IndexType &k,
		 const ForceIndexType &fidx);
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
void buildDeviceAngleList (const HostAngleList & hanglelist,
			   DeviceAngleList & danglelist);
void destroyDeviceAngleList (DeviceAngleList & danglelist) ;




#endif

