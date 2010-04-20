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
    
    IndexType totalNumCell;
    IndexType maxNumBond;
    IndexType * global_numBond;
    IndexType * global_neighborIndex;
    IndexType * global_bondIndex;

    IndexType * neighborIndex;
    // IndexType * bondActive;
private:
    void clear ();
    void easyMalloc (const IndexType & totalNumCell,
		     const IndexType & maxNumBond);
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
    IndexType stride () const {return maxNumBond;}
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
    
    IndexType totalNumCell;
    IndexType maxNumBond;
    IndexType * global_numBond;
    IndexType * global_neighborIndex;
    IndexType * global_bondIndex;

    IndexType * neighborIndex;
    // IndexType * bondActive;
private:
    void clear ();
    void easyMalloc (const IndexType & totalNumCell,
		     const IndexType & maxNumBond);
public:
    DeviceBondList ();
    ~DeviceBondList ();
public:
    IndexType stride () const {return maxNumBond;}
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

  namespace DeviceBondList_cudaDevice {
    __device__ IndexType
    indexConvert (const IndexType & stride,
		  const IndexType & localIndex,
		  const IndexType & ithBond)
    { return localIndex * stride + ithBond; }
  }
  
#endif
}


#endif
