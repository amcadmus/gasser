#ifndef __AngleList_interface_h_wanghan__
#define __AngleList_interface_h_wanghan__

#include "common.h"
#include "AngleList.h"
#include "MDSystem.h"
#include "AngleInteraction.h"

class AngleList 
{
  HostAngleList hanglelist;
  void sortAngle ();
  IndexType NAngleForce_mem;
  IndexType paramLength_mem;
public:
  DeviceAngleList danglelist;
public:
  IndexType NAngleForce;
  InteractionType * angleType;
  IndexType * paramPosi;
  ScalorType * param;
  IndexType paramLength;
public:
  AngleList () ;
  ~AngleList ();
  AngleList(const DeviceMDData & ddata)
      {init (ddata);}
  void init (const DeviceMDData & ddata);
  void addAngle (const IndexType & ii,
		 const IndexType & jj,
		 const IndexType & kk,
		 const AngleInteractionParameter & param_);
  void build ();
};



#endif
