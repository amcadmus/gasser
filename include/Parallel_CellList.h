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

  
  class HostSubCellList : std::vector<IndexType > 
  {
    HostCellListedMDData * ptr_data;
    
    // IndexType * begin;
    // IndexType * end;
  //   class Iterator : public std::vector<IndexType >::iterator
  //   {
  // public:
  //     // IndexType cellStartIndex ()
  //     // 	  {return Parallel::Environment::numAtomInCell() * **this;}
  //     // IndexType numAtomInCell ();
  //   };
  //   friend class Iterator;
public:
    HostSubCellList () :ptr_data(NULL) {};
    ~HostSubCellList () {};
public:
    void build ();
    bool isBuilt ();
    void bond (HostCellListedMDData & ddata);
    void add (const HostSubCellList & a);
    void sub (const HostSubCellList & a);    
  };


  
#ifdef DEVICE_CODE
  class DeviceSubCellList;
  
  class DeviceCellListedMDData : public DeviceMDData 
  {
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
		       DeviceSubCellList & subList);
    void buildSubListLevel1 (const IndexType & xIdLo,
			     const IndexType & xIdUp,
			     const IndexType & yIdLo,
			     const IndexType & yIdUp,
			     const IndexType & zIdLo,
			     const IndexType & zIdUp,
			     DeviceSubCellList & subList);    
  };

  
  class DeviceSubCellList : std::vector<IndexType > 
  {
    friend  class DeviceCellListedMDData;
    DeviceCellListedMDData * ptr_data;
    
    // IndexType * begin;
    // IndexType * end;
  //   class Iterator : public std::vector<IndexType >::iterator
  //   {
  // public:
  //     // IndexType cellStartIndex ()
  //     // 	  {return Parallel::Environment::numAtomInCell() * **this;}
  //     // IndexType numAtomInCell ();
  //   };
  //   friend class Iterator;
public:
    DeviceSubCellList ();
    ~DeviceSubCellList () {}
public:
    void bond (DeviceCellListedMDData & ddata);
    void build ();
    bool isBuilt ();
    void add (const DeviceSubCellList & a);
    void sub (const DeviceSubCellList & a);    
  };

  
  // namespace GPUGlobalCall{
#endif // DEVICE_CODE
    
}

    


#endif

