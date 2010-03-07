#define MPI_CODE

#include "Parallel_Interface.h"
#include "Parallel_Environment.h"
#include "compile_error_mixcode.h"

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
  Environment env;
  return env.myRank();
}

int Parallel::Interface::
numProc()
{
  Environment env;  
  return env.numProc();
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
