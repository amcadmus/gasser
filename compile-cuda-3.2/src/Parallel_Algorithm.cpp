#include <algorithm>
#include "Parallel_Algorithm.h"
#include "common.h"
#include "compile_error_mixcode.h"


void Parallel::Interface::
sort (std::vector<IndexType>::iterator a, std::vector<IndexType>::iterator b)
{
  std::sort (a, b);
}


bool Parallel::Interface::
is_sorted (std::vector<IndexType>::iterator a, std::vector<IndexType>::iterator b)
{
  bool state = true;
  if (a == b) return state;
  std::vector<IndexType>::iterator i (a), j(a);
  ++j;
  for (; j != b; ++i, ++j){
    if (*j > *i) {
      state = false;
      break;
    }
  }
  return state;
}

void Parallel::Interface::
unique (std::vector<IndexType>::iterator a, std::vector<IndexType>::iterator b)
{
  std::unique (a, b);
}

std::vector<IndexType >::iterator Parallel::Interface::
set_difference (std::vector<IndexType>::const_iterator a0,
		std::vector<IndexType>::const_iterator a1,
		std::vector<IndexType>::const_iterator b0,
		std::vector<IndexType>::const_iterator b1,
		std::vector<IndexType>::iterator c0)
{
  return std::set_difference (a0, a1, b0, b1, c0);
}

std::vector<IndexType >::iterator Parallel::Interface::
copy (std::vector<IndexType>::const_iterator b0,
      std::vector<IndexType>::const_iterator b1,
      std::vector<IndexType>::iterator c0)
{
  return std::copy (b0, b1, c0);
}

