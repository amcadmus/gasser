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
    void writeGroFile (FILE * fp,
		       const IndexType num,
		       const IndexType * resdindex,
		       const char * resdname,
		       const char * atomname,
		       const IndexType * atomindex,
		       const ScalorType * posix,
		       const ScalorType * posiy,
		       const ScalorType * posiz,
		       const ScalorType * velox,
		       const ScalorType * veloy,
		       const ScalorType * veloz,
		       const ScalorType boxx,  
		       const ScalorType boxy,  
		       const ScalorType boxz);
};




#endif
