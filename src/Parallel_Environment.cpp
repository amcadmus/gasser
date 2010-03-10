#define MPI_CODE
#include "Parallel_Environment.h"
#include "compile_error_mixcode.h"

// int Parallel::Environment::
// myRank  ()  {return myRank_;}
// int Parallel::Environment::
// numProc ()  {return numProc_;}


MPI_Comm Parallel::Environment::commCart;
int Parallel::Environment::dims[3] = {0, 0, 0};
int Parallel::Environment::myRank_ = 0;
int Parallel::Environment::numProc_ = 0;


void Parallel::Environment::
init (int * argc, char *** argv)
{
  int flag;
  MPI_Initialized (&flag);
  if (flag) return;

  MPI_Init (argc, argv);

  MPI_Comm_size (MPI_COMM_WORLD, &numProc_);
  MPI_Comm_rank (MPI_COMM_WORLD, &myRank_);  
}

void Parallel::Environment::
finalize ()
{
  MPI_Finalize();
}

void Parallel::Environment::
initCart  (const int & nx,
	   const int & ny,
	   const int & nz)
{
  dims[0] = nx;
  dims[1] = ny;
  dims[2] = nz;
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
		int & rank) 
{
  int coord[3];
  coord[CoordXIndex] = ix;
  coord[CoordYIndex] = iy;
  coord[CoordZIndex] = iz;
  MPI_Cart_rank (commCart, coord, &rank);
}

void Parallel::Environment::
rankToCartCoord (const int & rank,
		 int & ix,
		 int & iy,
		 int & iz) 
{
  int coord[3];
  MPI_Cart_coords (commCart, rank, 3, coord);
  ix = coord[CoordXIndex];
  iy = coord[CoordYIndex];
  iz = coord[CoordZIndex];
}

void Parallel::Environment::
numProcDim (int & nx, int & ny, int & nz) 
{
  nx = dims[CoordXIndex];
  ny = dims[CoordYIndex];
  nz = dims[CoordZIndex];
}

void Parallel::Environment::
barrier ()
{
  MPI_Barrier (commCart);
}


void Parallel::Environment::
neighborProcIndex (int direction,
		   int displacement,
		   int & src,
		   int & dest)
{
  MPI_Cart_shift (Parallel::Environment::commCart,
		  direction,
		  displacement,
		  & src, & dest);
}

