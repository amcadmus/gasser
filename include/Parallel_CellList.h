#ifndef __Parallel_CellList_h__
#define __Parallel_CellList_h__

#include "common.h"
#include "Parallel_MDData.h"
#include "Parallel_TransferEngineCompatible.h"
#include "MDError_interface.h"
#include "SystemBondedInteraction.h"

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
  class DeviceCellRelation;
  
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
public:
    const HostVectorType & getFrameUp  () const {return frameUp;}
    const HostVectorType & getFrameLow () const {return frameLow;}
    const IndexType & getDevideLevel () const {return devideLevel;}
    const HostIntVectorType & getNumCell () const {return numCell;}
    const ScalorType & getRlist () const {return rlist;}
    IndexType * cptr_numAtomInCell () {return numAtomInCell;}
    const IndexType * cptr_numAtomInCell () const {return numAtomInCell;}
public:
    void initTopology (const Topology::System & sysTop,
		       const SystemBondedInteraction & sysBdInter);
public:
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
    void clearData ();
    void clearData (const SubCellList & subList);
    void mallocFromHost (const HostCellListedMDData & hdata);
    void copy (const HostCellListedMDData & hdata,
	       const MDDataItemMask_t mask = MDDataItemMask_All);
    void add  (const HostCellListedMDData & hdata,
	       const MDDataItemMask_t mask = MDDataItemMask_All);
    void writeData_SimpleFile (const char * filename);
public:
    void buildSubList (const IndexType & xIdLo,
		       const IndexType & xIdUp,
		       const IndexType & yIdLo,
		       const IndexType & yIdUp,
		       const IndexType & zIdLo,
		       const IndexType & zIdUp,
		       SubCellList & subList) const;
    void buildSubList (const IndexType & xIdLo,
		       const IndexType & xIdUp,
		       const IndexType & yIdLo,
		       const IndexType & yIdUp,
		       const IndexType & zIdLo,
		       const IndexType & zIdUp,
		       HostSubCellList & subList) ;
    inline void buildSubListAllCell    (SubCellList & subList) const;
    inline void buildSubListRealCell   (SubCellList & subList) const;
    inline void buildSubListGhostCell  (SubCellList & subList) const;
    inline void buildSubListInnerShell (SubCellList & subList) const;
    inline void buildSubListAllCell    (HostSubCellList & subList) ;
    inline void buildSubListRealCell   (HostSubCellList & subList) ;
    inline void buildSubListGhostCell  (HostSubCellList & subList) ;
    inline void buildSubListInnerShell (HostSubCellList & subList) ;
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

  class HostCellRelation
  {
    friend class DeviceCellRelation;
    
    IndexType totalNumCell;
    IndexType MaxNeiPerCell;
    IndexType * numNeighbor;
    IndexType * neighborCellIndex;
    HostCoordNoiType * neighborShiftNoi;
    void clear ();
    void easyMalloc (const IndexType & totalNumCell,
		     const IndexType & MaxNeiPerCell);
public:
    HostCellRelation ();
    ~HostCellRelation ();
public:
    IndexType * cptr_numNeighborCell   () {return numNeighbor;}
    IndexType * cptr_neighborCellIndex () {return neighborCellIndex;}
    HostCoordNoiType * cptr_neighborShiftNoi () {return neighborShiftNoi;}
    const IndexType * cptr_numNeighborCell   () const {return numNeighbor;}
    const IndexType * cptr_neighborCellIndex () const {return neighborCellIndex;}
    const HostCoordNoiType * cptr_neighborShiftNoi () const {return neighborShiftNoi;}
    IndexType stride_neighborCellIndex () const {return MaxNeiPerCell;}
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
    IndexType getMaxNum () const ;
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
    void build (const IndexType & maxNumTrans);
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
    HostTransferPackage (const Parallel::HostAllocator::HostMallocType_t & t =
			 Parallel::HostAllocator::hostMallocCDefault);
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
  class DeviceBondList;
  
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
    IndexType * forwardMap_step1;
    IndexType * forwardMap_step2;
    IndexType memSize;
    
    // IndexType * numNeighborCell;
    // IndexType * neighborCellIndex;
    // IndexType maxNumNeighborCell;

    mutable MDError err;
private:
    void easyMallocCell (const IndexType & totalNumCell);
    void initZeroCell ();
    void clearCell ();
    void calNumCell (const ScalorType & rlist,
		     const IndexType & devideLevel,
		     const BoxDirection_t & bdir,
		     IntVectorType & numCell,
		     VectorType & result_frameLow,
		     VectorType & result_frameUp);
public:
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
    DeviceCellListedMDData () ;
    ~DeviceCellListedMDData () ;
    const VectorType & getFrameUp  () const {return frameUp;}
    const VectorType & getFrameLow () const {return frameLow;}
    const IndexType  & getDevideLevel () const {return devideLevel;}
    const IntVectorType & getNumCell () const {return numCell;}
    const ScalorType & getRlist () const {return rlist;}
    IndexType * dptr_numAtomInCell () {return numAtomInCell;}
    const IndexType * dptr_numAtomInCell () const {return numAtomInCell;}
public:
    void initCellStructure (const ScalorType & rlist,
			    const IndexType & devideLevel = 1,
			    const BoxDirection_t & bdir = 7);
    bool reinitCellStructure (const ScalorType & rlist,
			      const IndexType & devideLevel = 1,
			      const BoxDirection_t & bidr = 7);
    void rebuild ();
    void rebuild (DeviceBondList & dbdlist);
    void applyPeriodicBondaryCondition ();
    void applyPeriodicBondaryConditionOnGhostCells ();
    void rescaleCoordinate (const HostVectorType & scale);
    void rescaleVelocity   (const HostVectorType & scale);
public:
    void clearData (const SubCellList & subList);
public:
    void mallocFromHost   (const HostCellListedMDData & hdata);
    void mallocToHost     (HostCellListedMDData & hdata) const;
    void mallocFromDevice (const DeviceCellListedMDData & ddata);
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
		       SubCellList & subList) const;
    inline void buildSubListAllCell    (SubCellList & subList) const;
    inline void buildSubListRealCell   (SubCellList & subList) const;
    inline void buildSubListGhostCell  (SubCellList & subList) const;
    inline void buildSubListInnerShell (SubCellList & subList) const;
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
  

  class DeviceCellRelation 
  {
    const DeviceCellListedMDData * ptr_list;
    
    bool		 malloced;
    IndexType		 * numNeighbor;
    IndexType		 * neighborCellIndex;
    CoordNoiType	 * neighborShiftNoi;
    IndexType		 totalNumCell;
    IndexType		 MaxNeiPerCell;
    void clear ();
    void easyMalloc (const IndexType & totalNumCell,
		     const IndexType & MaxNeiPerCell);
public:
    DeviceCellRelation  () : ptr_list(NULL), malloced(false) {};
    ~DeviceCellRelation ();
    void rebuild (const DeviceCellListedMDData & list);
    void rebuild (const DeviceCellListedMDData & list,
		  const SubCellList & sub0,
		  const SubCellList & sub1);
public:
    IndexType * dptr_numNeighborCell () {return numNeighbor;}
    IndexType * dptr_neighborCellIndex () {return neighborCellIndex;}
    CoordNoiType * dptr_neighborShiftNoi () {return neighborShiftNoi;}
    const IndexType * dptr_numNeighborCell () const {return numNeighbor;}
    const IndexType * dptr_neighborCellIndex () const {return neighborCellIndex;}
    const CoordNoiType * dptr_neighborShiftNoi () const {return neighborShiftNoi;}
    IndexType stride_neighborCellIndex () const {return MaxNeiPerCell;}
public:
    void copyToHost (HostCellRelation & hrelation);
  };
  
#endif // DEVICE_CODE
    
}


void Parallel::HostCellListedMDData::
buildSubListAllCell (SubCellList & subList) const
{
  buildSubList (0, getNumCell().x,
		0, getNumCell().y,
		0, getNumCell().z,
		subList);
}

  
void Parallel::HostCellListedMDData::
buildSubListRealCell  (SubCellList & subList) const
{
  IndexType devideLevel = getDevideLevel();
  buildSubList (devideLevel, getNumCell().x - devideLevel,
		devideLevel, getNumCell().y - devideLevel,
		devideLevel, getNumCell().z - devideLevel,
		subList);
}

      
void Parallel::HostCellListedMDData::
buildSubListGhostCell (SubCellList & subList) const
{
  buildSubListAllCell (subList);
  SubCellList temp;
  buildSubListRealCell (temp);
  subList.sub(temp);
}

void Parallel::HostCellListedMDData::
buildSubListInnerShell (SubCellList & subList) const
{
  IndexType devideLevel = getDevideLevel();
  buildSubList (devideLevel, getNumCell().x - devideLevel,
		devideLevel, getNumCell().y - devideLevel,
		devideLevel, getNumCell().z - devideLevel,
		subList);
  SubCellList temp;
  devideLevel *= 2;
  buildSubList (devideLevel, getNumCell().x - devideLevel,
		devideLevel, getNumCell().y - devideLevel,
		devideLevel, getNumCell().z - devideLevel,
		temp);
  subList.sub (temp);
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

void Parallel::HostCellListedMDData::
buildSubListInnerShell (HostSubCellList & subList)
{
  subList.setHostData (*this);
  IndexType devideLevel = getDevideLevel();
  buildSubList (devideLevel, getNumCell().x - devideLevel,
		devideLevel, getNumCell().y - devideLevel,
		devideLevel, getNumCell().z - devideLevel,
		(SubCellList &)subList);
  SubCellList temp;
  devideLevel *= 2;
  buildSubList (devideLevel, getNumCell().x - devideLevel,
		devideLevel, getNumCell().y - devideLevel,
		devideLevel, getNumCell().z - devideLevel,
		temp);
  subList.sub (temp);
}


#ifdef DEVICE_CODE
void Parallel::DeviceCellListedMDData::
buildSubListAllCell (SubCellList & subList) const
{
  buildSubList (0, getNumCell().x,
		0, getNumCell().y,
		0, getNumCell().z,
		subList);
}

  
void Parallel::DeviceCellListedMDData::
buildSubListRealCell  (SubCellList & subList) const
{
  IndexType devideLevel = getDevideLevel();
  buildSubList (devideLevel, getNumCell().x - devideLevel,
		devideLevel, getNumCell().y - devideLevel,
		devideLevel, getNumCell().z - devideLevel,
		subList);
}

      
void Parallel::DeviceCellListedMDData::
buildSubListGhostCell (SubCellList & subList) const
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

void Parallel::DeviceCellListedMDData::
buildSubListInnerShell (SubCellList & subList) const
{
  IndexType devideLevel = getDevideLevel();
  buildSubList (devideLevel, getNumCell().x - devideLevel,
		devideLevel, getNumCell().y - devideLevel,
		devideLevel, getNumCell().z - devideLevel,
		subList);
  SubCellList temp;
  devideLevel *= 2;
  buildSubList (devideLevel, getNumCell().x - devideLevel,
		devideLevel, getNumCell().y - devideLevel,
		devideLevel, getNumCell().z - devideLevel,
		temp);
  subList.sub (temp);
}


#endif // DEVICE_CODE




#endif

