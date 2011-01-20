#ifndef __Parallel_MDSystem_h_wanghan__
#define __Parallel_MDSystem_h_wanghan__

#include "Parallel_MDData.h"
#include "Parallel_CellList.h"
#include "Parallel_Timer.h"

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


    IndexType maxNum_xsend0;
    IndexType maxNum_xrecv0;
    IndexType maxNum_xsend1;
    IndexType maxNum_xrecv1;
    IndexType maxNum_ysend0;
    IndexType maxNum_yrecv0;
    IndexType maxNum_ysend1;
    IndexType maxNum_yrecv1;
    IndexType maxNum_zsend0;
    IndexType maxNum_zrecv0;
    IndexType maxNum_zsend1;
    IndexType maxNum_zrecv1;

    IndexType thisNum_xsend0;
    IndexType thisNum_xrecv0;
    IndexType thisNum_xsend1;
    IndexType thisNum_xrecv1;
    IndexType thisNum_ysend0;
    IndexType thisNum_yrecv0;
    IndexType thisNum_ysend1;
    IndexType thisNum_yrecv1;
    IndexType thisNum_zsend0;
    IndexType thisNum_zrecv0;
    IndexType thisNum_zsend1;
    IndexType thisNum_zrecv1;

private:
    void calTransNumX ();
    void calTransNumY ();
    void calTransNumZ ();
public:
    SystemRedistributeTransferUtil ();
    ~SystemRedistributeTransferUtil ();
    void clear ();
    void setHostData (HostCellListedMDData & hdata,
		      HostCellListedMDData & buffdata);
    void redistributeHost ();
  };


  class SystemTransCoordsTransferUtil
  {
    HostCellListedMDData * ptr_hdata;
    HostCellListedMDData * ptr_buff ;

    MDDataItemMask_t mask;

    HostSubCellList xsend0;
    HostSubCellList xrecv0;
    HostSubCellList xsend1;
    HostSubCellList xrecv1;
    HostSubCellList ysend0;
    HostSubCellList yrecv0;
    HostSubCellList ysend1;
    HostSubCellList yrecv1;
    HostSubCellList zsend0;
    HostSubCellList zrecv0;
    HostSubCellList zsend1;
    HostSubCellList zrecv1;
    
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

    IndexType maxNum_xsend0;
    IndexType maxNum_xrecv0;
    IndexType maxNum_xsend1;
    IndexType maxNum_xrecv1;
    IndexType maxNum_ysend0;
    IndexType maxNum_yrecv0;
    IndexType maxNum_ysend1;
    IndexType maxNum_yrecv1;
    IndexType maxNum_zsend0;
    IndexType maxNum_zrecv0;
    IndexType maxNum_zsend1;
    IndexType maxNum_zrecv1;

    IndexType thisNum_xsend0;
    IndexType thisNum_xrecv0;
    IndexType thisNum_xsend1;
    IndexType thisNum_xrecv1;
    IndexType thisNum_ysend0;
    IndexType thisNum_yrecv0;
    IndexType thisNum_ysend1;
    IndexType thisNum_yrecv1;
    IndexType thisNum_zsend0;
    IndexType thisNum_zrecv0;
    IndexType thisNum_zsend1;
    IndexType thisNum_zrecv1;

private:
    void calTransNumX ();
    void calTransNumY ();
    void calTransNumZ ();
public:
    SystemTransCoordsTransferUtil ();
    ~SystemTransCoordsTransferUtil ();
    void clear ();
    void setHostData (HostCellListedMDData & hdata,
		      HostCellListedMDData & buffdata);
    void transCoords ();
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
    SubCellList			hostSubInner;
    SubCellList			hostSubOuter;
    SubCellList			deviceSubInner;
    SubCellList			deviceSubOuter;
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

  class SystemTransCoordsCopyUtil
  {
    SubCellList			hostSubInner;
    SubCellList			hostSubOuter;
    SubCellList			deviceSubInner;
    SubCellList			deviceSubOuter;
    HostTransferPackage		hpkgInner;
    DeviceTransferPackage	dpkgInner;
    HostTransferPackage		hpkgOuter;
    DeviceTransferPackage	dpkgOuter;
    HostCellListedMDData *	ptr_hdata;
    DeviceCellListedMDData *	ptr_ddata;
    MDDataItemMask_t		mask;
public:
    SystemTransCoordsCopyUtil ();
    void setData (HostCellListedMDData & hdata,
		  DeviceCellListedMDData & ddata);
    inline void copyToHost ();
    inline void copyFromHost ();
    inline void clearGhost ();
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
    IndexType				numFreedom;
    GlobalHostMDData			globalHostData;
    HostCellListedMDData		localHostData;
    HostCellListedMDData		hostBuff;
    DeviceCellListedMDData		deviceData;
    DeviceCellRelation			cellRelation;
private:
    SystemRedistributeTransferUtil	redistribtransUtil;
    SystemRedistributeCopyUtil		redistribcopyUtil;
    SystemTransCoordsTransferUtil	transCoordstransUtil;
    SystemTransCoordsCopyUtil		transCoordscopyUtil;
    SystemCollectDataUtil		collectUtil;
public:
    MDSystem ();
    ~MDSystem();
    /** 
     * Initialize the MD system
     * 
     * @param confFileName	file that provides the configuration
     * @param sysTop		system topology			
     * @param cellSize		the expected cell size (the actual cell size may
     *				be larger
     * @param divideLevel	the divide level of the cell
     */
    void init (const char * confFileName,
	       const Topology::System & sysTop,
	       const ScalorType & cellSize,
	       const IndexType  & divideLevel = 1);
    /** 
     * When the size of the system box, expected cell size or the divide level
     * changes, reinit cell structure accordingly. 
     *
     * This function first calculates the number of cells on each
     * direction. If any number is different than the old system, it
     * renints the cell structure, otherwise it will not do
     * anything. This function only applies to the device MD data. A
     * additional updateHost is needed to update the host data.
     * 
     * @param cellSize		the expected cell size
     * @param divideLevel	the divide level of the cell
     * 
     * @return			true if the cell structure changes,
     *				false otherwise.
     */
    bool reinitCellStructure (const ScalorType & cellSize,
			      const IndexType  & divideLevel = 1);
    /** 
     * Finalize the system.
     * 
     */
    void finalize ();
public:
    /** 
     * Update local host MD data from device MD data
     * 
     */
    void updateHost ()
	{ deviceData.copyToHost (localHostData); }
    /** 
     * Collect local host MD data to global MD data.
     * 
     */
    void collectLocalData ();
    void redistribute ();
    void transferGhost ();
    void clearGhost ();
public:
    const IndexType &  getNumFreedom () const {return numFreedom;}
public:
    IndexType numAtomInGroFile (const char * filename) 
	{ return globalHostData.numAtomInGroFile(filename); }
    void writeLocalData_SimpleFile (const char * filename)
	{ localHostData.writeData_SimpleFile(filename); }
    void writeGlobalData_GroFile (const char * filename);
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


#ifdef DEVICE_CODE  
  
  void Parallel::SystemTransCoordsCopyUtil::
  copyToHost ()
  {
    Parallel::Timer::HostTimer::tic (Parallel::Timer::item_TransferGhost_DHCopy_Pack);
    dpkgInner.pack (*ptr_ddata, mask);
    Parallel::Timer::HostTimer::toc (Parallel::Timer::item_TransferGhost_DHCopy_Pack);
    Parallel::Timer::HostTimer::tic (Parallel::Timer::item_TransferGhost_DHCopy_Copy);
    dpkgInner.copyToHost (hpkgInner);
    Parallel::Timer::HostTimer::toc (Parallel::Timer::item_TransferGhost_DHCopy_Copy);
    Parallel::Timer::HostTimer::tic (Parallel::Timer::item_TransferGhost_DHCopy_Unpack);
    hpkgInner.unpack_replace (*ptr_hdata);
    Parallel::Timer::HostTimer::toc (Parallel::Timer::item_TransferGhost_DHCopy_Unpack);
  }

  void Parallel::SystemTransCoordsCopyUtil::
  copyFromHost ()
  {
    Parallel::Timer::HostTimer::tic (Parallel::Timer::item_TransferGhost_DHCopy_Pack);
    hpkgOuter.pack (*ptr_hdata, mask);
    Parallel::Timer::HostTimer::toc (Parallel::Timer::item_TransferGhost_DHCopy_Pack);
    Parallel::Timer::HostTimer::tic (Parallel::Timer::item_TransferGhost_DHCopy_Copy);
    dpkgOuter.copyFromHost (hpkgOuter);
    Parallel::Timer::HostTimer::toc (Parallel::Timer::item_TransferGhost_DHCopy_Copy);
    Parallel::Timer::HostTimer::tic (Parallel::Timer::item_TransferGhost_DHCopy_Unpack);
    dpkgOuter.unpack_replace (*ptr_ddata);
    Parallel::Timer::HostTimer::toc (Parallel::Timer::item_TransferGhost_DHCopy_Unpack);
  }

  void Parallel::SystemTransCoordsCopyUtil::
  clearGhost()
  {
    ptr_ddata->clearData (deviceSubOuter);
  }
  
#endif // DEVICE_CODE

}




#endif
