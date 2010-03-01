#ifndef __AngleList_h_wanghan__
#define __AngleList_h_wanghan__

#include "common.h"
#include "AngleInteraction.h"


struct HostAngleList 
{
  IndexType stride; // is the expected larger than or equal to the number of Atoms
  // IndexType listLength;
  IndexType maxNumAngle;
  IndexType * angleNeighborIndex;
  IndexType * angleIndex;
  IndexType * anglePosi;
  IndexType * numAngle;
private:
  void clearMem ();
public:
  HostAngleList ();
  ~HostAngleList ();
  void clearAngle ();
  void reinit (const IndexType & stride,
	       const IndexType & maxNumAngle);
  void addAngle (const IndexType &i,
		 const IndexType &j,
		 const IndexType &k,
		 const IndexType &angleIndex,
		 const IndexType &anglePosi);
};


struct DeviceAngleList 
{
  bool malloced ;
  IndexType stride; // is the expected larger than or equal to the number of Atoms
  IndexType maxNumAngle;
  IndexType * angleNeighborIndex;
  IndexType * angleIndex;
  IndexType * anglePosi;
  IndexType * numAngle;
};

void initDeviceAngleList (DeviceAngleList & dbdlist) ;
void mallocDeviceAngleList (const HostAngleList & hbdlist,
			    DeviceAngleList & dbdlist);
void copyDeviceAngleList (const HostAngleList & hbdlist,
			  DeviceAngleList & dbdlist);
void copyDeviceAngleList (const DeviceAngleList & dbdlist1,
			  DeviceAngleList & dbdlist);
void destroyDeviceAngleList (DeviceAngleList & dbdlist) ;


#endif

