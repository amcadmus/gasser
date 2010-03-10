#ifndef __Parallel_Algorithm_h_wanghan__
#define __Parallel_Algorithm_h_wanghan__

#include <vector>
#include "common.h"

namespace Parallel{
  namespace Interface {
    void sort (std::vector<IndexType>::iterator a,
	       std::vector<IndexType>::iterator b);
    void unique (std::vector<IndexType>::iterator a,
		 std::vector<IndexType>::iterator b);
    bool is_sorted (std::vector<IndexType>::iterator a,
		    std::vector<IndexType>::iterator b);
    std::vector<IndexType >::iterator
    set_difference (std::vector<IndexType>::const_iterator a0,
		    std::vector<IndexType>::const_iterator a1,
		    std::vector<IndexType>::const_iterator b0,
		    std::vector<IndexType>::const_iterator b1,
		    std::vector<IndexType>::iterator c0);
    std::vector<IndexType >::iterator 
    copy (std::vector<IndexType>::const_iterator b0,
	  std::vector<IndexType>::const_iterator b1,
	  std::vector<IndexType>::iterator c0);    
  }
}



#endif

