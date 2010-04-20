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
    
    IndexType myStride;
    IndexType maxLength;
    IndexType * global_numBond;
    IndexType * global_neighborIndex;
    IndexType * global_bondIndex;

    IndexType * neighborIndex;
    // IndexType * bondActive;
private:
    void clear ();
    void easyMalloc (const IndexType & stride,
		     const IndexType & length);
    void initZero ();
    IndexType convertIndex (const IndexType & localIndex,
			    const IndexType & ithBond) const;
public:
    HostBondList ();
    ~HostBondList ();
public:
    void reinit (HostCellListedMDData & hcellData,
		 Topology::System & sysTop,
		 SystemBondedInteraction & sysBdInter);
public:
    IndexType stride () const {return myStride;}
    IndexType & item_global_numBond      (const IndexType & localIndex); 
    IndexType & item_global_neighborIndex(const IndexType & localIndex,
					  const IndexType & ithBond);
    IndexType & item_global_bondIndex    (const IndexType & localIndex,
					  const IndexType & ithBond);
    IndexType & item_neighborIndex       (const IndexType & localIndex,
					  const IndexType & ithBond);
    // IndexType & item_bondActive          (const IndexType & localIndex,
    // 					  const IndexType & ithBond);

    const IndexType & item_global_numBond      (const IndexType & localIndex) const; 
    const IndexType & item_global_neighborIndex(const IndexType & localIndex,
						const IndexType & ithBond) const;
    const IndexType & item_global_bondIndex    (const IndexType & localIndex,
						const IndexType & ithBond) const;
    const IndexType & item_neighborIndex       (const IndexType & localIndex,
						const IndexType & ithBond) const;
    // const IndexType & item_bondActive          (const IndexType & localIndex,
    // 						const IndexType & ithBond) const;
  };					  

#ifdef DEVICE_CODE
  class DeviceBondList
  {
public:
    bool malloced;
    
    IndexType myStride;
    IndexType maxLength;
    IndexType * global_numBond;
    IndexType * global_neighborIndex;
    IndexType * global_bondIndex;

    IndexType * neighborIndex;
    // IndexType * bondActive;
private:
    void clear ();
    void easyMalloc (const IndexType & stride,
		     const IndexType & length);
public:
    DeviceBondList ();
    DeviceBondList (const HostBondList & hbdlist);
    ~DeviceBondList ();
public:
    void reinit (const HostBondList & hbdlist);
    IndexType stride () const {return myStride;}
    IndexType convertIndex (const IndexType & localIndex,
			    const IndexType & ithBond) const;
    IndexType * dptr_global_numBond       () {return global_numBond;}
    IndexType * dptr_global_neighborIndex () {return global_neighborIndex;}
    IndexType * dptr_global_bondIndex     () {return global_bondIndex;}
    IndexType * dptr_neighborIndex        () {return neighborIndex;}
    const IndexType * dptr_global_numBond       () const {return global_numBond;}
    const IndexType * dptr_global_neighborIndex () const {return global_neighborIndex;}
    const IndexType * dptr_global_bondIndex     () const {return global_bondIndex;}
    const IndexType * dptr_neighborIndex        () const {return neighborIndex;}
public:
    void copyFromHost (const HostBondList & hbdlist);
    void copyToHost (HostBondList & hbdlist) const;
  };

  void buildDeviceBondList (const DeviceCellListedMDData & ddata,
			    const DeviceCellRelation & relation,
			    DeviceBondList & dbdlist);

  namespace DeviceBondList_cudaDevice {
    __device__ IndexType
    indexConvert (const IndexType & stride,
		  const IndexType & localIndex,
		  const IndexType & ithBond)
    { return ithBond * stride + localIndex; }
  }
  
  namespace CudaGlobal{
    __global__ void 
    buildDeviceBondList (const IndexType * numAtomInCell,
			 const IndexType * globalIndex,
			 const IndexType * numNeighborCell,
			 const IndexType * neighborCellIndex,
			 const IndexType   cellRelationStride,
			 const IndexType * global_neighborIndex,
			 const IndexType * global_numBond,
			 const IndexType   bondListStride,
			 IndexType * neighborIndex);
  }

#endif
}


#endif
