#ifndef __Parallel_MDSystem_h_wanghan__
#define __Parallel_MDSystem_h_wanghan__

#include "Parallel_MDData.h"
#include "Parallel_CellList.h"

// #define DEVICE_CODE

namespace Parallel {
  class SystemTranferUtils
  {
    HostCellListedMDData * ptr_hdata;
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
private:
public:
    SystemTranferUtils ();
    void setHostData (HostCellListedMDData & hdata);
    void redistribute ();
  };
  

#ifdef DEVICE_CODE
  class MDSystem 
  {
    // SubCellList xsend0;
    // SubCellList xrecv0;
    // IndexType * recv0Num;
    // SubCellList xsend1;
    // SubCellList xrecv1;
    // SubCellList ysend0;
    // SubCellList yrecv0;
    // SubCellList ysend1;
    // SubCellList yrecv1;
    // SubCellList zsend0;
    // SubCellList zrecv0;
    // SubCellList zsend1;
    // SubCellList zrecv1;
private:
    // gro file related
    char * atomName;
    IndexType * atomIndex;
    char * resdName;
    IndexType * resdIndex;
    void reallocGroProperty (const IndexType & memSize);

    IndexType			globalNumAtom;
    GlobalHostMDData		globalHostData;
    HostCellListedMDData	localHostData;
    DeviceCellListedMDData	deviceData;
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
public:
    // void redistribute ();
    void tryHostSend ();
  };
#endif
  
}




#endif
