#ifndef __GPU_Environment_h__
#define __GPU_Environment_h__

namespace GPU{
  class Environment
  {
    int deviceId;
public:
    Environment();
    void setDeviceId (int id);
  };
}

#endif
