#ifndef __FileManager_h_wanghan__
#define __FileManager_h_wanghan__

#include "common.h"
#include <stdio.h>

IndexType  readAtomNameMapFile (const char * filename,
				const IndexType numAtom,
				const char * atomName,
				TypeType * type,
				ScalorType * mass,
				ScalorType * charge);

			 

#endif 
