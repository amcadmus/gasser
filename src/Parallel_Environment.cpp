#include "Parallel_Environment.h"

Parallel::Environment::
Environment (int * argc, char *** argv)
{
  int flag;
  MPI_Initialized (&flag);
  if (flag) return;

  MPI_Init (argc, argv);
}

Parallel::Environment::
~Environment ()
{
  MPI_Finalize();
}


void Parallel::Environment::
init (const int division[3])
{
  MPI_Comm_size (MPI_COMM_WORLD, &numProc);
  
  dims[0] = division[0];
  dims[1] = division[1];
  dims[2] = division[2];
  if (dims[0] * dims[1] * dims[2] != numProc){
    throw MDExcptDimsNotConsistentWithNProc ();
  }

  int periodic [3];
  periodic[0] = periodic[1] = periodic[2] = 1;
  int reorder = 1;
  int ndims = 3;
  MPI_Cart_create (MPI_COMM_WORLD, ndims, dims, periodic, reorder, &commCart);

  MPI_Comm_rank (commCart, &myRank);
  printf ("# myrank is %d in %d\n", myRank, numProc);
}

