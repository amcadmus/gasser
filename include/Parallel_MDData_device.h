#ifndef __Parallel_MDData_device_h__
#define __Parallel_MDData_device_h__

namespace Parallel {
  namespace CudaGlobal {
    __global__ void
    initZeroDeviceData (const IndexType num,
			CoordType  * coord,
			IntScalorType * coordNoix,
			IntScalorType * coordNoiy,
			IntScalorType * coordNoiz,
			ScalorType * velox,
			ScalorType * veloy,
			ScalorType * veloz,
			IndexType  * globalIndex,
			TypeType   * type,
			ScalorType * mass,
			ScalorType * charge);
  }
}




#endif
