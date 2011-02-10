#ifndef __GPU_Environment_h__
#define __GPU_Environment_h__

#define MaxLengthDeviceName			256
#define MaxNumDevicePerNode			32

namespace GPU{
  class Environment
  {
    static unsigned numThreadsInCell;
    static int my_deviceId;
    static int numActiveDevice;
    static int memActiveDevice;
    static int * activeDeviceId;
    static char activeDeviceName[MaxLengthDeviceName];
public:
    static void initialize (const unsigned & numThreadsInCell = 64,
			    const char * activeDeviceName = "Tesla C1060");
    static void finalize ();
    static int getNumActiveDevice () 
	{return GPU::Environment::numActiveDevice;}
    static const int * cptr_activeDeviceId () 
	{return GPU::Environment::activeDeviceId;}
    static void setDeviceId (int id);
    static int  getDeviceId () {return GPU::Environment:: my_deviceId;}
    static int  getNumThreadsInCell () {return numThreadsInCell;}
  };
}

#endif
