#ifndef __DataTransferBlock_h_wanghan__
#define __DataTransferBlock_h_wanghan__

#define NumDataItemsInMDData 14

namespace Parallel{
  struct DataTransferBlock 
  {
    void * pointer[NumDataItemsInMDData];
    size_t size [NumDataItemsInMDData];
  };
}



#endif


