#ifndef __Parallel_MDSystem_h_wanghan__
#define __Parallel_MDSystem_h_wanghan__

#include "Parallel_MDData.h"
#include "Parallel_CellList.h"

#define DEVICE_CODE

namespace Parallel {
    
  class MDSystem 
  {
    SubCellList xsend0;
    SubCellList xrecv0;
    IndexType * recv0Num;
    SubCellList xsend1;
    SubCellList xrecv1;
    SubCellList ysend0;
    SubCellList yrecv0;
    SubCellList ysend1;
    SubCellList yrecv1;
    SubCellList zsend0;
    SubCellList zrecv0;
    SubCellList zsend1;
    SubCellList zrecv1;
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
