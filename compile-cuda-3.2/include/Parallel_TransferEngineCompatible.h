#ifndef __Parallel_TransferEngineCompatible_h__
#define __Parallel_TransferEngineCompatible_h__

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>


namespace Parallel{
  
  class TransferEngineCompatible 
  {
public:
    virtual void getTransBuffs (IndexType * num,
				void *** buffs,
				size_t ** sizes) = 0;
  };
  
}

#endif
