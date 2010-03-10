#ifndef __Parallel_MDSystem_h_wanghan__
#define __Parallel_MDSystem_h_wanghan__

#include "Parallel_MDData.h"
#include "Parallel_CellList.h"

#define DEVICE_CODE

namespace Parallel {
    
  class MDSystem 
  {
    HostSubCellList xsend0;
    HostSubCellList xrecv0;
    IndexType * recv0Num;
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
    
  };
}




#endif
