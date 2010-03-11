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

  class SubHostCellList ;
  
  class HostCellListedMDData : public HostMDData
  {
    ScalorType rlist;
    IndexType devideLevel;
    HostVectorType frameLow;
    HostVectorType frameUp;
    HostIntVectorType numCell;

    IndexType * numAtomInCell;

    IndexType * numNeighborCell;
    IndexType * neighborCellIndex;
    IndexType maxNumNeighborCell;
private:
    void reallocAll (const IndexType & totalNumCell,
		     const IndexType & maxNumNeighborCell);
public:
    HostCellListedMDData ();
    ~HostCellListedMDData ();
    // HostCellList ();
    // void formCellStructure (const ScalorType & rlist,
    // 			    const IndexType & devideLevel = 1,
    // 			    const BoxDirection_t & bdir = 7);
    void buildSubList (const IndexType & xIdLo,
		       const IndexType & xIdUp,
		       const IndexType & yIdLo,
		       const IndexType & yIdUp,
		       const IndexType & zIdLo,
		       const IndexType & zIdUp,
		       SubHostCellList & subList);
    void buildSubListLevel1 (const IndexType & xIdLo,
			     const IndexType & xIdUp,
			     const IndexType & yIdLo,
			     const IndexType & yIdUp,
			     const IndexType & zIdLo,
			     const IndexType & zIdUp,
			     SubHostCellList & subList);
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



  class HostTransferapackage : public HostMDData
  {
    IndexType numCell;
    IndexType memSize;
    IndexType * cellIndex;
    IndexType * cellStartIndex;
    inline void clearMe ();
    inline void easyMallocMe (IndexType memSize);
public:
    HostTransferapackage ();
    ~HostTransferapackage ();
public:
    void reinit (const SubCellList & SubCellList);
    void packAtom (HostCellListedMDData & hdata);
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

    IndexType * numNeighborCell;
    IndexType * neighborCellIndex;
    IndexType maxNumNeighborCell;

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
    // void buildSubListLevel1 (const IndexType & xIdLo,
    // 			     const IndexType & xIdUp,
    // 			     const IndexType & yIdLo,
    // 			     const IndexType & yIdUp,
    // 			     const IndexType & zIdLo,
    // 			     const IndexType & zIdUp,
    // 			     SubCellList & subList);    
  };


  class DeviceTransferPackage : public DeviceMDData
  {
    IndexType numCell;
    IndexType memSize;
    IndexType * hcellIndex;
    IndexType * cellIndex;
    IndexType * hcellStartIndex;
    IndexType * cellStartIndex;
    inline void clearMe ();
    inline void easyMallocMe (IndexType memSize);
public:
    DeviceTransferPackage ();
    ~DeviceTransferPackage ();
public:
    void reinit (const SubCellList & subCellList);
    void pack (DeviceCellListedMDData & ddata,
	       const MDDataItemMask_t mask = MDDataItemMask_All);
  };
  
  
#endif // DEVICE_CODE
    
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

