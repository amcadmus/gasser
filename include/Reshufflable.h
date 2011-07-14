#ifndef __Reshufflable_h_wanghan__
#define __Reshufflable_h_wanghan__

#include "systemDefines.h"
#include "MDTimer_interface.h"

/// The class that can be reshuffled.

/**
 * Reshuffling can improve the efficient a lot.
 */

class Reshufflable
{
 /** 
  * Reshuffle according to the index table that tells how the atom
  * indexes are mapped.
  * 
  * @param indexTable The index table.
  * @param numAtom Number of atom in the index table.
  * @param timer Timer measuring the performance.
  */
  virtual void reshuffle (const IndexType * indexTable,
			  const IndexType & numAtom,
			  MDTimer * timer) = 0;
};

#endif
