#define MPI_CODE

#include "Parallel_Interface.h"
#include "Parallel_Environment.h"
#include "compile_error_mixcode.h"

#define NUMTHREADSINCELL 20

unsigned Parallel::Interface::
numThreadsInCell ()
{
  return NUMTHREADSINCELL;
}

void Parallel::Interface::
initEnvironment (int * argc, char ***argv)
{
  Parallel::Environment::init (argc, argv);
}

void Parallel::Interface::
finalizeEnvironment ()
{
  Parallel::Environment::finalize ();
}

int Parallel::Interface::
myRank()
{
  return Parallel::Environment::myRank();
}

int Parallel::Interface::
numProc()
{
  return Parallel::Environment::numProc();
}

void Parallel::Interface::
initCart (const int & nx,
	  const int & ny,
	  const int & nz)
{
  Parallel::Environment::initCart(nx, ny, nz);
}

void Parallel::Interface::
cartCoordToRank (const int & ix,
		 const int & iy,
		 const int & iz,
		 int & rank )
{
  Parallel::Environment::cartCoordToRank (ix, iy, iz, rank);
}

void Parallel::Interface::
rankToCartCoord (const int & rank,
		 int & ix,
		 int & iy,
		 int & iz) 
{
  Parallel::Environment::rankToCartCoord (rank, ix, iy, iz);
}

void Parallel::Interface::
numProcDim (int & nx,
	    int & ny,
	    int & nz)
{
  Parallel::Environment::numProcDim (nx, ny, nz);
}

void Parallel::Interface::
barrier ()
{
  Parallel::Environment::barrier();
}

void Parallel::Interface::
shiftNeighbor (int direction,
	       int displacement,
	       int & src,
	       int & dest)
{
  Parallel::Environment::neighborProcIndex (direction, displacement, src, dest);
}

