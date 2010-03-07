#ifndef __Parallel_MDSystem_h_wanghan__
#define __Parallel_MDSystem_h_wanghan__

#include "Parallel_MDData.h"

#define DEVICE_CODE

namespace Parallel {
    
  class MDSystem 
  {
    HostCellList & myCellList;
    SubHostCellList xsend0;
    SubHostCellList xrecv0;
    SubHostCellList xsend1;
    SubHostCellList xrecv1;
    SubHostCellList ysend0;
    SubHostCellList yrecv0;
    SubHostCellList ysend1;
    SubHostCellList yrecv1;
    SubHostCellList zsend0;
    SubHostCellList zrecv0;
    SubHostCellList zsend1;
    SubHostCellList zrecv1;
private:
    // gro file related
    char * atomName;
    IndexType * atomIndex;
    char * resdName;
    IndexType * resdIndex;
    void reallocGroProperty (const IndexType & memSize);

    IndexType  globalNumAtom;
    GlobalHostMDData globalHostData;
    HostMDData localHostData;
    DeviceMDData deviceData;
    HostMDData sendHostData;
public:
    MDSystem (HostCellList & myCellList);
    ~MDSystem();
    // void initConf_GroFile (const char * filename);
    IndexType numAtomInGroFile (const char * filename) 
	{ return globalHostData.numAtomInGroFile(filename); }
    void init (const char * confFileName,
	       const Topology::System & sysTop);
    void writeLocalData_SimpleFile (const char * filename)
	{ localHostData.writeData_SimpleFile(filename); }
public:
    void redistribute ();
    
  };
}




#endif
