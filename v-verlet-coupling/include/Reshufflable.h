#ifndef __Reshufflable_h_wanghan__
#define __Reshufflable_h_wanghan__

#include "systemDefines.h"
#include "MDTimer_interface.h"

class Reshufflable
{
  virtual void reshuffle (const IndexType * indexTable,
			  const IndexType & numAtom,
			  MDTimer * timer) = 0;
};

#endif
