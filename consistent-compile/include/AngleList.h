#ifndef __AngleList_h_wanghan__
#define __AngleList_h_wanghan__

#include "common.h"
#include "AngleInteraction.h"

/// Angle list on host.

struct HostAngleList 
{
  IndexType stride; /**< The stride of the list, which is the expected
			 * larger than or equal to the number of
			 * Atoms. */
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

/// Angle list on device.

struct DeviceAngleList 
{
  bool malloced ;		/**< Tell us if it is allocated. */
  IndexType stride; 	/**< The stride of the list, which is the
			 * expected larger than or equal to the number
			 * of Atoms. */
  IndexType maxNumAngle;	/**< Maximum number of angle
				 * interaction per atom. */
  IndexType * angleNeighborIndex; /**< Neighbor index of the angle interaction. */
  IndexType * angleIndex;	/**< Index of the angle interaction. */
  IndexType * anglePosi;	/**< Position of atom in the angle interaction. */
  IndexType * numAngle;		/**< The number of angle interaction of atoms. */
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

