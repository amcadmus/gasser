#ifndef __Parallel_HostAllocator_h_wanghan__
#define __Parallel_HostAllocator_h_wanghan__


namespace Parallel{  
  namespace HostAllocator{
    enum HostMallocType {
      hostMallocCDefault			= 1,
      hostMallocPageLocked			= 2
    };
    typedef enum HostMallocType HostMallocType_t;

    void hostMalloc (void** ptr,
		     size_t size,
		     HostMallocType_t type = hostMallocCDefault);
    void hostFree   (void** ptr,
		     HostMallocType_t type = hostMallocCDefault);
  }
}


#endif
