#ifndef __Parallel_MDSystem_h_wanghan__
#define __Parallel_MDSystem_h_wanghan__

#include "Parallel_MDData.h"
#include "Parallel_CellList.h"

// #define DEVICE_CODE

namespace Parallel {
  class SystemRedistributeTransferUtil
  {
    HostCellListedMDData * ptr_hdata;
    HostCellListedMDData * ptr_buff ;

    MDDataItemMask_t mask;

    HostSubCellList xsend0;
    HostSubCellList xrecv0;
    HostSubCellList xrecv0h;
    HostSubCellList xsend1;
    HostSubCellList xrecv1;
    HostSubCellList xrecv1h;
    HostSubCellList ysend0;
    HostSubCellList yrecv0;
    HostSubCellList yrecv0h;
    HostSubCellList ysend1;
    HostSubCellList yrecv1;
    HostSubCellList yrecv1h;
    HostSubCellList zsend0;
    HostSubCellList zrecv0;
    HostSubCellList zrecv0h;
    HostSubCellList zsend1;
    HostSubCellList zrecv1;
    HostSubCellList zrecv1h;
    
    int xdest0;
    int xsrc0;
    int xdest1;
    int xsrc1;
    int ydest0;
    int ysrc0;
    int ydest1;
    int ysrc1;
    int zdest0;
    int zsrc0;
    int zdest1;
    int zsrc1;

    TransNumAtomInSubList xsendNum0;
    TransNumAtomInSubList xrecvNum0;
    TransNumAtomInSubList xsendNum1;
    TransNumAtomInSubList xrecvNum1;
    TransNumAtomInSubList ysendNum0;
    TransNumAtomInSubList yrecvNum0;
    TransNumAtomInSubList ysendNum1;
    TransNumAtomInSubList yrecvNum1;
    TransNumAtomInSubList zsendNum0;
    TransNumAtomInSubList zrecvNum0;
    TransNumAtomInSubList zsendNum1;
    TransNumAtomInSubList zrecvNum1;

    TransSubListData xsendData0;
    TransSubListData xrecvData0;
    TransSubListData xsendData1;
    TransSubListData xrecvData1;    
    TransSubListData ysendData0;
    TransSubListData yrecvData0;
    TransSubListData ysendData1;
    TransSubListData yrecvData1;
    TransSubListData zsendData0;
    TransSubListData zrecvData0;
    TransSubListData zsendData1;
    TransSubListData zrecvData1;

private:
public:
    SystemRedistributeTransferUtil ();
    void setHostData (HostCellListedMDData & hdata,
		      HostCellListedMDData & buffdata);
    void redistributeHost ();
  };


  class SystemCollectDataUtil
  {
    HostCellListedMDData * ptr_hdata;
    HostCellListedMDData * ptr_buff ;
    HostMDData * ptr_gdata;
    
    MDDataItemMask_t mask;

    HostSubCellList sendlist;
    HostSubCellList recvlist;
    TransNumAtomInSubList sendNum;
    TransNumAtomInSubList recvNum;
    TransSubListData sendData;
    TransSubListData recvData;
private:
    void addBuffToGlobal ();
public:
    SystemCollectDataUtil ();
    void setHostData (HostCellListedMDData & hdata,
		      HostCellListedMDData & buffdata,
		      HostMDData & globalData);
    void collect ();
  };
  
  

#ifdef DEVICE_CODE

  class SystemRedistributeCopyUtil
  {
    SubCellList hostSubInner;
    SubCellList hostSubOuter;
    SubCellList deviceSubInner;
    SubCellList deviceSubOuter;
    
    HostTransferPackage		hpkgInner;
    DeviceTransferPackage	dpkgInner;
    HostTransferPackage		hpkgOuter;
    DeviceTransferPackage	dpkgOuter;
    HostCellListedMDData *	ptr_hdata;
    DeviceCellListedMDData *	ptr_ddata;
    MDDataItemMask_t		mask;
public:
    SystemRedistributeCopyUtil ();
    void setData (HostCellListedMDData & hdata,
		  DeviceCellListedMDData & ddata);
    inline void clearDeviceSent ();
    inline void copyToHost ();
    inline void copyFromHost ();
  };
  
  
  class MDSystem 
  {
public:
    // gro file related
    char * atomName;
    IndexType * atomIndex;
    char * resdName;
    IndexType * resdIndex;
    void reallocGroProperty (const IndexType & memSize);

    IndexType				globalNumAtom;
    GlobalHostMDData			globalHostData;
    HostCellListedMDData		localHostData;
    HostCellListedMDData		hostBuff;
    DeviceCellListedMDData		deviceData;
private:
    SystemRedistributeTransferUtil	redistribtransUtil;
    SystemRedistributeCopyUtil		redistribcopyUtil;
    SystemCollectDataUtil		collectUtil;
public:
    MDSystem ();
    ~MDSystem();
    // void initConf_GroFile (const char * filename);
    IndexType numAtomInGroFile (const char * filename) 
	{ return globalHostData.numAtomInGroFile(filename); }
    void init (const char * confFileName,
	       const Topology::System & sysTop);
    void writeLocalData_SimpleFile (const char * filename)
	{ localHostData.writeData_SimpleFile(filename); }
    void writeGlobalData_GroFile (const char * filename);
public:
    void updateHost ()
	{ deviceData.copyToHost (localHostData); }
    void collectLocalData ();
    void redistribute ();
  };
#endif // DEVICE_CODE

#ifdef DEVICE_CODE
  void Parallel::SystemRedistributeCopyUtil::
  copyToHost ()
  {
    dpkgOuter.pack (*ptr_ddata, mask);
    dpkgOuter.copyToHost (hpkgOuter);
    hpkgOuter.unpack_replace (*ptr_hdata);
  }

  void Parallel::SystemRedistributeCopyUtil::
  copyFromHost ()
  {
    hpkgInner.pack (*ptr_hdata, mask);
    dpkgInner.copyFromHost (hpkgInner);
    dpkgInner.unpack_add (*ptr_ddata);
  }

  void Parallel::SystemRedistributeCopyUtil::
  clearDeviceSent ()
  {
    ptr_ddata->clearData (deviceSubOuter);
  }
#endif // DEVICE_CODE
  
}




#endif
