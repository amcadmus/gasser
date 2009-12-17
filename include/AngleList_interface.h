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
  mdAngleInteraction_t * angleType;
  IndexType * paramPosi;
  ScalorType * param;
  IndexType paramLength;
public:
  AngleList () ;
  ~AngleList ();
  AngleList(const DeviceMDData & ddata, const IndexType & listLength) 
      {init (ddata, listLength);}
  void init (const DeviceMDData & ddata, const IndexType & listLength);
  void addAngle (const IndexType & ii,
		 const IndexType & jj,
		 const IndexType & kk,
		 const mdAngleInteraction_t & type,
		 const ScalorType * param);
  void build ();
};



#endif
