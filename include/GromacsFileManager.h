#ifndef __GromacsFileManager_h_wanghan__
#define __GromacsFileManager_h_wanghan__

#include "common.h"

namespace GromacsFileManager{
    void readGroFile (const char * filename,
		      IndexType * resdindex,
		      char * resdname,
		      char * atomname,
		      IndexType * atomindex,
		      ScalorType * posix,
		      ScalorType * posiy,
		      ScalorType * posiz,
		      ScalorType * velox,
		      ScalorType * veloy,
		      ScalorType * veloz,
		      ScalorType * boxx,  
		      ScalorType * boxy,  
		      ScalorType * boxz);
};




#endif
