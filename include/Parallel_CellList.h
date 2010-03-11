#ifndef __Parallel_CellList_h__
#define __Parallel_CellList_h__

#include "common.h"
#include "Parallel_MDData.h"
#include "MDError_interface.h"

namespace Parallel{
  
  // class MDExcptTooFewCells	: public MDException {};
  // class MDExcptWrongAtomInProc	: public MDException {};
  // class MDExcptWrongCellIndex	: public MDException {};
  
  class MDExcptCellList : public MDException {
    char message[MaxExceptionMsgLength];
public:
    MDExcptCellList (const char * description) 
	{strncpy (message, description, MaxExceptionMsgLength);}
    virtual const char* what() const throw()
	{return message;}
  };

  class SubCellList;
  
  class HostCellListedMDData : public HostMDData
  {
    friend class HostTransferPackage;
    
    ScalorType rlist;
    IndexType devideLevel;
    HostVectorType frameLow;
    HostVectorType frameUp;
    HostIntVectorType numCell;

    IndexType * numAtomInCell;

    // IndexType * numNeighborCell;
    // IndexType * neighborCellIndex;
    // IndexType maxNumNeighborCell;
public:
    void reallocAll (const IndexType & totalNumCell,
		     const IndexType & maxNumNeighborCell);
    IndexType * getNumAtomInCell () {return numAtomInCell;}
public:
    HostCellListedMDData ();
    ~HostCellListedMDData ();
    // HostCellList ();
    // void formCellStructure (const ScalorType & rlist,
    // 			    const IndexType & devideLevel = 1,
    // 			    const BoxDirection_t & bdir = 7);
    const HostVectorType & getFrameUp  () const {return frameUp;}
    const HostVectorType & getFrameLow () const {return frameLow;}
    const IndexType & getDevideLevel () const {return devideLevel;}
    const HostIntVectorType & getNumCell () const {return numCell;}
    const ScalorType & getRlist () const {return rlist;}
    IndexType D3toD1 (const IndexType & ix,
		      const IndexType & iy,
		      const IndexType & iz) const
	{return iz + numCell.z * (iy + numCell.y * ix);}
    void D1toD3 (const IndexType & i,
		 IndexType & x,
		 IndexType & y,
		 IndexType & z)
	{IndexType tmp = i;  z = tmp % (numCell.z); tmp = (tmp - z) / numCell.z;
	  y = tmp % (numCell.y); x = (tmp - y) / numCell.y;}
public:
    void buildSubList (const IndexType & xIdLo,
		       const IndexType & xIdUp,
		       const IndexType & yIdLo,
		       const IndexType & yIdUp,
		       const IndexType & zIdLo,
		       const IndexType & zIdUp,
		       SubCellList & subList);
    inline void buildSubListAllCell   (SubCellList & subList);
    inline void buildSubListRealCell  (SubCellList & subList);
    inline void buildSubListGhostCell (SubCellList & subList);
  };
  
  class SubCellList : public std::vector<IndexType > 
  {
public:
    SubCellList ();
    ~SubCellList () {}
public:
    void build ();
    bool isBuilt ();
    void add (const SubCellList & a);
    void sub (const SubCellList & a);    
  };



  class HostTransferPackage : public HostMDData
  {
    IndexType numCell;
    IndexType memSize;
    IndexType * cellIndex;
    IndexType * cellStartIndex;
    MDDataItemMask_t myMask;
public:
    inline void clearMe ();
    void easyMallocMe (IndexType memSize);
public:
    HostTransferPackage ();
    ~HostTransferPackage ();
public:
    const IndexType & getTotalNumCell   () const {return numCell;}
    const IndexType & getMemSize        () const {return memSize;}
    const MDDataItemMask_t & getMask	() const {return myMask;}
    const IndexType * getCellIndex      () const {return cellIndex;}
    const IndexType * getCellStartIndex () const {return cellStartIndex;}
    IndexType & getTotalNumCell   () {return numCell;}
    IndexType & getMemSize        () {return memSize;}
    MDDataItemMask_t & getMask	  () {return myMask;}
    IndexType * getCellIndex      () {return cellIndex;}
    IndexType * getCellStartIndex () {return cellStartIndex;}
public:
    void reinit (const SubCellList & SubCellList);
    void pack (const HostCellListedMDData & hdata,
	       const MDDataItemMask_t mask = MDDataItemMask_All);
  };  

  
  
#ifdef DEVICE_CODE
  class SubCellList;
  class DeviceTransferPackage ;
  
  class DeviceCellListedMDData : public DeviceMDData 
  {
    friend class DeviceTransferPackage ;
    
    bool malloced ;
    
    ScalorType rlist;
    IndexType devideLevel;
    IntVectorType numCell;
    VectorType frameLow;
    VectorType frameUp;

    IndexType * numAtomInCell;

    // IndexType * numNeighborCell;
    // IndexType * neighborCellIndex;
    // IndexType maxNumNeighborCell;

    MDError err;
private:
    void mallocCell (const IndexType & totalNumCell,
		     const IndexType & maxNumNeighborCell);
    void initZeroCell ();
    void clearCell ();
public:
    DeviceCellListedMDData () ;
    ~DeviceCellListedMDData () ;
    const VectorType & getFrameUp  () const {return frameUp;}
    const VectorType & getFrameLow () const {return frameLow;}
    const IndexType  & getDevideLevel () const {return devideLevel;}
    const IntVectorType & getNumCell () const {return numCell;}
    const ScalorType & getRlist () const {return rlist;}
    IndexType D3toD1 (const IndexType & ix,
		      const IndexType & iy,
		      const IndexType & iz) const
	{return iz + numCell.z * (iy + numCell.y * ix);}
    void D1toD3 (const IndexType & i,
		 IndexType & x,
		 IndexType & y,
		 IndexType & z)
	{IndexType tmp = i;  z = tmp % (numCell.z); tmp = (tmp - z) / numCell.z;
	  y = tmp % (numCell.y); x = (tmp - y) / numCell.y;}

public:
    void initCellStructure (const ScalorType & rlist,
			    const IndexType & devideLevel = 1,
			    const BoxDirection_t & bdir = 7);
    void rebuild ();
public:
    void copyToHost   (HostCellListedMDData & hdata) const;
    void copyFromHost (const HostCellListedMDData & hdata);
public:
    void buildSubList (const IndexType & xIdLo,
		       const IndexType & xIdUp,
		       const IndexType & yIdLo,
		       const IndexType & yIdUp,
		       const IndexType & zIdLo,
		       const IndexType & zIdUp,
		       SubCellList & subList);
    inline void buildSubListAllCell   (SubCellList & subList);
    inline void buildSubListRealCell  (SubCellList & subList);
    inline void buildSubListGhostCell (SubCellList & subList);
  };


  class DeviceTransferPackage : public DeviceMDData
  {
    IndexType numCell;
    IndexType memSize;
    IndexType * hcellIndex;
    IndexType * cellIndex;
    IndexType * hcellStartIndex;
    IndexType * cellStartIndex;
    MDDataItemMask_t myMask;
    mutable MDError err;
    inline void clearMe ();
    inline void easyMallocMe (IndexType memSize);
public:
    DeviceTransferPackage ();
    ~DeviceTransferPackage ();
public:
    const MDDataItemMask_t & getMask () const {return myMask;}
    void reinit (const SubCellList & subCellList);
    void pack (const DeviceCellListedMDData & ddata,
	       const MDDataItemMask_t mask = MDDataItemMask_All);
    void unpack_replace (DeviceCellListedMDData & ddata) const;
    void unpack_add     (DeviceCellListedMDData & ddata) const;
public:
    void copyToHost   (HostTransferPackage & hpkg) const;
    void copyFromHost (const HostTransferPackage & hpkg);
  };
  
  
#endif // DEVICE_CODE
    
}


void Parallel::HostCellListedMDData::
buildSubListAllCell (SubCellList & subList)
{
  buildSubList (0, getNumCell().x,
		0, getNumCell().y,
		0, getNumCell().z,
		subList);
}

  
void Parallel::HostCellListedMDData::
buildSubListRealCell  (SubCellList & subList)
{
  IndexType devideLevel = getDevideLevel();
  buildSubList (devideLevel, getNumCell().x - devideLevel,
		devideLevel, getNumCell().y - devideLevel,
		devideLevel, getNumCell().z - devideLevel,
		subList);
}

      
void Parallel::HostCellListedMDData::
buildSubListGhostCell (SubCellList & subList)
{
  buildSubListAllCell (subList);
  SubCellList temp;
  buildSubListRealCell (temp);
  subList.sub(temp);
}


#ifdef DEVICE_CODE
void Parallel::DeviceCellListedMDData::
buildSubListAllCell (SubCellList & subList)
{
  buildSubList (0, getNumCell().x,
		0, getNumCell().y,
		0, getNumCell().z,
		subList);
}

  
void Parallel::DeviceCellListedMDData::
buildSubListRealCell  (SubCellList & subList)
{
  IndexType devideLevel = getDevideLevel();
  buildSubList (devideLevel, getNumCell().x - devideLevel,
		devideLevel, getNumCell().y - devideLevel,
		devideLevel, getNumCell().z - devideLevel,
		subList);
}

      
void Parallel::DeviceCellListedMDData::
buildSubListGhostCell (SubCellList & subList)
{
  buildSubListAllCell (subList);
  SubCellList temp;
  buildSubListRealCell (temp);
  subList.sub(temp);

  // int count = 0;
  // for (IndexType i = 0; i < subList.size(); ++i){
  //   IndexType ix, iy, iz;
  //   D1toD3 (subList[i], ix, iy, iz);
  //   printf ("%d %d %d %d\n", ++count, ix, iy, iz);
  // }  
}
#endif // DEVICE_CODE




#endif

