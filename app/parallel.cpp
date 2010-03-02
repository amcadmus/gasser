#define CPP_FILE

#include "common.h"
#include "Parallel_Environment.h"


int main(int argc, char * argv[])
{
  Parallel::Environment env (&argc, &argv);

  int dim[3];
  dim[0] = 2;
  dim[1] = dim[2] = 1;

  env.init (dim);

  
  return 0;
}

