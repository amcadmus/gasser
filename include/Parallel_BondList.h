#ifndef __Paralle_BondeList_h_wanghan__
#define __Paralle_BondeList_h_wanghan__

#include "common.h"
#include "Parallel_CellList.h"
#include "SystemBondedInteraction.h"

namespace Parallel{
  class DeviceBondList;

  class HostBondList
  {
    friend class DeviceBondList;
    
    IndexType stride;
    IndexType maxNumBond;
    IndexType maxNumAngle;
    IndexType maxNumDihedral;
    IndexType * bondNeighbor_localIndex;
    IndexType * angleNeighbor_localIndex;
    IndexType * dihedralNeighbor_localIndex;
private:
    void clear ();
    void easyMalloc (const IndexType & stride,
		     const IndexType & maxNumBond,
		     const IndexType & maxNumAngle,
		     const IndexType & maxNumDihedral);
    void fillzero ();
public:
    HostBondList ();
    ~HostBondList ();
public:
    IndexType getStride () const {return stride;}
    IndexType getMaxNumBond () const {return maxNumBond;}
    IndexType getMaxNumAngle () const {return maxNumAngle;}
    IndexType getMaxNumDihedral () const {return maxNumDihedral;}
    IndexType * cptr_bondNeighbor_localIndex ()
	{return bondNeighbor_localIndex;}
    IndexType * cptr_angleNeighbor_localIndex ()
	{return angleNeighbor_localIndex;}
    IndexType * cptr_dihedralNeighbor_localIndex ()
	{return dihedralNeighbor_localIndex;}
    const IndexType * cptr_bondNeighbor_localIndex () const
	{return bondNeighbor_localIndex;}
    const IndexType * cptr_angleNeighbor_localIndex () const
	{return angleNeighbor_localIndex;}
    const IndexType * cptr_dihedralNeighbor_localIndex () const
	{return dihedralNeighbor_localIndex;}
  };					  

#ifdef DEVICE_CODE
  class DeviceBondList
  {
public:
    bool malloced;
    
    IndexType stride;
    IndexType maxNumBond;
    IndexType maxNumAngle;
    IndexType maxNumDihedral;
    IndexType * bondNeighbor_localIndex;
    IndexType * angleNeighbor_localIndex;
    IndexType * dihedralNeighbor_localIndex;
private:
    void clear ();
    void easyMalloc (const IndexType & stride,
		     const IndexType & maxNumBond,
		     const IndexType & maxNumAngle,
		     const IndexType & maxNumDihedral);
public:
    DeviceBondList ();
    DeviceBondList (const DeviceCellListedMDData & data);
    ~DeviceBondList ();
public:
    void reinit (const DeviceCellListedMDData & data);
    IndexType getStride () const {return stride;}
    IndexType getMaxNumBond () const {return maxNumBond;}
    IndexType getMaxNumAngle () const {return maxNumAngle;}
    IndexType getMaxNumDihedral () const {return maxNumDihedral;}
    IndexType * dptr_bondNeighbor_localIndex () {return bondNeighbor_localIndex;}
    IndexType * dptr_angleNeighbor_localIndex () {return angleNeighbor_localIndex;}
    IndexType * dptr_dihedralNeighbor_localIndex () {return dihedralNeighbor_localIndex;}
    const IndexType * dptr_bondNeighbor_localIndex () const
	{return bondNeighbor_localIndex;}
    const IndexType * dptr_angleNeighbor_localIndex () const
	{return angleNeighbor_localIndex;}
    const IndexType * dptr_dihedralNeighbor_localIndex () const
	{return dihedralNeighbor_localIndex;}
public:
    void copyFromHost (const HostBondList & hbdlist);
    void copyToHost (HostBondList & hbdlist) const;
  };

  void buildDeviceBondList (const DeviceCellListedMDData & ddata,
			    const DeviceCellRelation & relation,
			    DeviceBondList & dbdlist);
  
  namespace CudaGlobal{
    __global__ void 
    buildDeviceBondList (const IndexType * numAtomInCell,
			 const IndexType * globalIndex,
			 const IndexType * numNeighborCell,
			 const IndexType * neighborCellIndex,
			 const IndexType   cellRelationStride,
			 const IndexType   maxNumBond,
			 const IndexType * numBond,
			 const IndexType * bondNeighbor_globalIndex,
			 IndexType * bondNeighbor_localIndex,
			 const IndexType   maxNumAngle,
			 const IndexType * numAngle,
			 const IndexType * angleNeighbor_globalIndex,
			 IndexType * angleNeighbor_localIndex,
			 const IndexType   maxNumDihedral,
			 const IndexType * numDihedral,
			 const IndexType * dihedralNeighbor_globalIndex,
			 IndexType * dihedralNeighbor_localIndex,
			 const IndexType   bondTopStride,
			 const IndexType   listTopStride);
  }

#endif
}


#endif
