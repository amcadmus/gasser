#ifndef __Parallel_MDSystem_h_wanghan__
#define __Parallel_MDSystem_h_wanghan__

#include "Parallel_Environment.h"
#include "Parallel_MDData.h"

namespace Parallel {
    
  class MDSystem 
  {
private:
    // gro file related
    char * atomName;
    IndexType * atomIndex;
    char * resdName;
    IndexType * resdIndex;
    void reallocGroProperty (const IndexType & memSize);

    const Environment * env;
    IndexType  globalNumAtom;
    GlobalHostMDData globalHostData;
    HostMDData sendHostData;
    HostMDData localHostData;
    DeviceMDData deviceData;
public:
    MDSystem (const Environment & env_);
    ~MDSystem();
    // void initConf_GroFile (const char * filename);
    IndexType numAtomInGroFile (const char * filename) 
	{ return globalHostData.numAtomInGroFile(filename); }
    void init (const char * confFileName,
	       const Topology::System & sysTop);
    void writeLocalData_SimpleFile (const char * filename)
	{ localHostData.writeData_SimpleFile(filename); }
  };
}




#endif
