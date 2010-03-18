#ifndef __Parallel_CellList_h__
#define __Parallel_CellList_h__

#include "common.h"
#include "Parallel_MDData.h"
#include "Parallel_TransferEngineCompatible.h"
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
  class DeviceCellListedMDData;
  class DeviceTransferPackage;
  class HostSubCellList;
  class DeviceSubCellList;
  
  class HostCellListedMDData : public HostMDData
  {
    friend class HostTransferPackage;
    friend class DeviceCellListedMDData;
    
    ScalorType rlist;
    IndexType devideLevel;
    HostVectorType frameLow;
    HostVectorType frameUp;
    HostIntVectorType numCell;

    IndexType * numAtomInCell;
    IndexType memSize;

    // IndexType * numNeighborCell;
    // IndexType * neighborCellIndex;
    // IndexType maxNumNeighborCell;
private:
    void easyReallocCell (const IndexType & totalNumCell);
public:
    HostCellListedMDData ();
    HostCellListedMDData (const HostCellListedMDData & hdata);
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
    IndexType * cptr_numAtomInCell () {return numAtomInCell;}
    const IndexType * cptr_numAtomInCell () const {return numAtomInCell;}
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
    void clearData ();
    void clearData (const SubCellList & subList);
    void copy (const HostCellListedMDData & hdata,
	       const MDDataItemMask_t mask = MDDataItemMask_All);
    void add  (const HostCellListedMDData & hdata,
	       const MDDataItemMask_t mask = MDDataItemMask_All);
public:
    void buildSubList (const IndexType & xIdLo,
		       const IndexType & xIdUp,
		       const IndexType & yIdLo,
		       const IndexType & yIdUp,
		       const IndexType & zIdLo,
		       const IndexType & zIdUp,
		       SubCellList & subList);
    void buildSubList (const IndexType & xIdLo,
		       const IndexType & xIdUp,
		       const IndexType & yIdLo,
		       const IndexType & yIdUp,
		       const IndexType & zIdLo,
		       const IndexType & zIdUp,
		       HostSubCellList & subList);
    inline void buildSubListAllCell   (SubCellList & subList);
    inline void buildSubListRealCell  (SubCellList & subList);
    inline void buildSubListGhostCell (SubCellList & subList);
    inline void buildSubListAllCell   (HostSubCellList & subList);
    inline void buildSubListRealCell  (HostSubCellList & subList);
    inline void buildSubListGhostCell (HostSubCellList & subList);
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

  class TransNumAtomInSubList : public TransferEngineCompatible 
  {
    void ** buffs;
    size_t * sizes;
    IndexType num;
    void clear();
public:
    TransNumAtomInSubList ();
    ~TransNumAtomInSubList ();
    void reinit (HostSubCellList & list);
    virtual void getTransBuffs (IndexType * num,
				void *** buffs,
				size_t ** sizes);
  };

  class TransSubListData : public TransferEngineCompatible
  {
    void ** buffs;
    size_t * sizes;
    IndexType num;
    MDDataItemMask_t myMask;
    HostSubCellList * ptr_list;
    void clear ();
public:
    TransSubListData ();
    ~TransSubListData ();
    void reinit (HostSubCellList & list,
		 const MDDataItemMask_t mask = MDDataItemMask_All);
    void build ();
    virtual void getTransBuffs (IndexType * num,
				void *** buffs,
				size_t ** sizes);
  };

  
  class HostSubCellList : public SubCellList
  {
    HostCellListedMDData * ptr_hdata;
public:
    HostSubCellList ()
	: ptr_hdata(NULL) {}
    HostSubCellList (HostCellListedMDData & hdata)
	: ptr_hdata(&hdata) {}
    void setHostData (HostCellListedMDData & hdata)
	{ ptr_hdata = &hdata; }
    HostCellListedMDData * host_ptr ()
	{ return ptr_hdata; }
    const HostCellListedMDData * host_ptr () const
	{ return ptr_hdata; }
    bool isBuilt ()
	{return ptr_hdata != NULL && SubCellList::isBuilt();}
public:
    void clearData ()
	{ ptr_hdata->clearData(*this);}
    void add (const HostSubCellList & clist,
	      const MDDataItemMask_t mask = MDDataItemMask_All);
    
// IndexType mallocSendRecvBuffRecord (const MDDataItemMask_t mask,
    // 					void *** buffs,
    // 					size_t ** sizes);
    // IndexType prepareSendData (const MDDataItemMask_t mask,
    // 			       IndexType * num,
    // 			       void *** buffs,
    // 			       size_t ** sizes);
    // IndexType prepareRecvData_replace (const MDDataItemMask_t mask,
    // 				       IndexType * numRecv,
    // 				       IndexType * num);
    
    // void mallocNumAtomTrans (IndexType ** numSent, size_t * size);
    // IndexType calculateNumAtomSend (IndexType * numSent);

  };

  
  
  class HostTransferPackage : public HostMDData
  {
    friend class DeviceTransferPackage;
    
    IndexType numCell;
    IndexType memSize;
    IndexType * cellIndex;
    IndexType * cellStartIndex;
    MDDataItemMask_t myMask;
private:
    inline void clearMe ();
    void easyMallocMe (IndexType memSize);
private:
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
    HostTransferPackage ();
    ~HostTransferPackage ();
public:
    void reinit (const SubCellList & SubCellList);
    void pack (const HostCellListedMDData & hdata,
	       const MDDataItemMask_t mask = MDDataItemMask_All);
    void unpack_replace (HostCellListedMDData & ddata) const;
    void unpack_add     (HostCellListedMDData & ddata) const;
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
    IndexType memSize;
    
    // IndexType * numNeighborCell;
    // IndexType * neighborCellIndex;
    // IndexType maxNumNeighborCell;

    mutable MDError err;
private:
    void easyMallocCell (const IndexType & totalNumCell);
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
    void clearData (const SubCellList & subList);
public:
    void copyToHost   (HostCellListedMDData & hdata,
		       const MDDataItemMask_t mask = MDDataItemMask_All) const;
    void copyFromHost (const HostCellListedMDData & hdata,
		       const MDDataItemMask_t mask = MDDataItemMask_All);
    void copyFromDevice (const DeviceCellListedMDData & ddata,
			 const MDDataItemMask_t = MDDataItemMask_All);
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

  
  class DeviceSubCellList : public SubCellList
  {
    DeviceCellListedMDData * ptr_ddata;
public:
    DeviceSubCellList ()
	: ptr_ddata(NULL) {}
    DeviceSubCellList (DeviceCellListedMDData & ddata)
	: ptr_ddata(&ddata) {}
    void setDeviceData (DeviceCellListedMDData & ddata)
	{ ptr_ddata = &ddata;}
    bool isBuilt ()
	{return ptr_ddata != NULL && SubCellList::isBuilt();}
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



void Parallel::HostCellListedMDData::
buildSubListAllCell (HostSubCellList & subList)
{
  subList.setHostData (*this);
  buildSubList (0, getNumCell().x,
		0, getNumCell().y,
		0, getNumCell().z,
		(SubCellList &)subList);
}

  
void Parallel::HostCellListedMDData::
buildSubListRealCell  (HostSubCellList & subList)
{
  subList.setHostData (*this);
  IndexType devideLevel = getDevideLevel();
  buildSubList (devideLevel, getNumCell().x - devideLevel,
		devideLevel, getNumCell().y - devideLevel,
		devideLevel, getNumCell().z - devideLevel,
		(SubCellList &)subList);
}

      
void Parallel::HostCellListedMDData::
buildSubListGhostCell (HostSubCellList & subList)
{
  subList.setHostData (*this);
  buildSubListAllCell ((SubCellList &)subList);
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

