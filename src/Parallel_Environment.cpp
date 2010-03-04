#include "Parallel_Environment.h"




Parallel::Environment::
Environment (int * argc, char *** argv)
{
  int flag;
  MPI_Initialized (&flag);
  if (flag) return;

  MPI_Init (argc, argv);

  MPI_Comm_size (MPI_COMM_WORLD, &numProc_);
  MPI_Comm_rank (MPI_COMM_WORLD, &myRank_);  
}

Parallel::Environment::
~Environment ()
{
  MPI_Finalize();
}


void Parallel::Environment::
init (const int division[3])
{
  dims[0] = division[0];
  dims[1] = division[1];
  dims[2] = division[2];
  if (dims[0] * dims[1] * dims[2] != numProc_){
    throw MDExcptDimsNotConsistentWithNProc ();
  }

  int periodic [3];
  periodic[0] = periodic[1] = periodic[2] = 1;
  int reorder = 1;
  int ndims = 3;
  MPI_Cart_create (MPI_COMM_WORLD, ndims, dims, periodic, reorder, &commCart);

  MPI_Comm_rank (commCart, &myRank_);
  // printf ("# myrank is %d in %d\n", myRank_, numProc_);
}

void Parallel::Environment::
cartCoordToRank(const int & ix,
		const int & iy,
		const int & iz,
		int & rank) const
{
  int coord[3];
  coord[CoordXIndex] = ix;
  coord[CoordYIndex] = iy;
  coord[CoordZIndex] = iz;
  MPI_Cart_rank (commCart, coord, &rank);
}

void Parallel::Environment::
randToCartCoord (const int & rank,
		 int & ix,
		 int & iy,
		 int & iz) const
{
  int coord[3];
  MPI_Cart_coords (commCart, rank, 3, coord);
  ix = coord[CoordXIndex];
  iy = coord[CoordYIndex];
  iz = coord[CoordZIndex];
}

void Parallel::Environment::
numProcDim (int & nx, int & ny, int & nz) const
{
  nx = dims[CoordXIndex];
  ny = dims[CoordYIndex];
  nz = dims[CoordZIndex];
}

