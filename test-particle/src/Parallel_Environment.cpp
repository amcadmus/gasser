#define MPI_CODE
#include "Parallel_Environment.h"
#include "GPU_Environment.h"
#include "compile_error_mixcode.h"

// int Parallel::Environment::
// myRank  ()  {return myRank_;}
// int Parallel::Environment::
// numProc ()  {return numProc_;}


MPI_Comm Parallel::Environment::commCart;
MPI_Comm Parallel::Environment::commActive;
int Parallel::Environment::dims[3] = {0, 0, 0};
int Parallel::Environment::myRank_ = 0;
int Parallel::Environment::numProc_ = 0;
int Parallel::Environment::active = 0;
// int Parallel::Environment::rank_worldRoot = 0;
int Parallel::Environment::inited = 0;
unsigned Parallel::Environment::cellCapacity = 0;

struct StrCmp 
{
  bool operator () (const char * a, const char * b)
      { return strcmp (a, b);}
};  

void Parallel::Environment::
initCart  (const int & nx,
	   const int & ny,
	   const int & nz)
{
  int ndims = 3;
  dims[0] = nx;
  dims[1] = ny;
  dims[2] = nz;
  // printf ("%d %d %d %d\n",
  // 	  Parallel::Environment::numProc(), nx, ny, nz);

  MPI_Dims_create (numProc(), ndims, dims);
  if (dims[0] * dims[1] * dims[2] != numProc()){
    throw Parallel::MDExcptDimsNotConsistentWithNProc ();
  }

  if (myRank () == 0){
    printf ("# ndims are %d %d %d\n", dims[0], dims[1], dims[2]);
  }
  
  int periodic [3];
  periodic[0] = periodic[1] = periodic[2] = 1;
  int reorder = 1;
  MPI_Cart_create (commActive, ndims, dims, periodic, reorder, &commCart);

  // printf ("# myrank is %d in %d\n", myRank_, numProc_);
}

void Parallel::Environment::
init_mpi (int * argc, char *** argv)
{
  int flag;
  MPI_Initialized (&flag);
  if (! flag) {
    MPI_Init (argc, argv);
  }
  else{
    std::cout << "duplicate call of Parallel::Environment::init_mpi may cause problems"
	      << std::endl;
    return;
  }
}

void Parallel::Environment::
init_env (const unsigned & cellCapacity_,
	  const char * deviceName,
	  const int & nx,
	  const int & ny,
	  const int & nz)
{
  if (inited){
    std::cout << "duplicate call of Parallel::Environment::init_env may cause problems"
	      << std::endl;
    return;
  }    
  cellCapacity = cellCapacity_;
  
  char my_name [MPI_MAX_PROCESSOR_NAME];
  int numProc_world;
  int myRank_world;
  MPI_Comm_size (MPI_COMM_WORLD, &numProc_world);
  MPI_Comm_rank (MPI_COMM_WORLD, &myRank_world);
  char * all_name;
  all_name = (char *) malloc (
      MPI_MAX_PROCESSOR_NAME * numProc_world * sizeof(char));
  if (all_name == NULL) {
    std::cerr << "not enough space\n";
    abort (1);
  }  
  int * all_numActiveDevice;
  all_numActiveDevice = (int *) malloc (numProc_world * sizeof(int));
  if (all_numActiveDevice == NULL){
    std::cerr << "not enough space\n";
    abort (1);
  }
  int * all_activity;
  all_activity = (int *) malloc (numProc_world * sizeof(int));
  if (all_activity == NULL){
    std::cerr << "not enough space\n";
    abort (1);
  }
  
  {
    int length;
    MPI_Get_processor_name (my_name, &length);
    MPI_Datatype dataType_procName;
    MPI_Type_contiguous (MPI_MAX_PROCESSOR_NAME, MPI_CHAR, &dataType_procName);
    MPI_Type_commit (&dataType_procName);
    MPI_Allgather (my_name,  1, dataType_procName,
		   all_name, 1, dataType_procName,
		   MPI_COMM_WORLD);
    MPI_Type_free(&dataType_procName);

    // char deviceName[] ="Tesla C1060";
    // char deviceName[] = "Device Emulation (CPU)";
    GPU::Environment::initialize (cellCapacity, deviceName);
    int my_numActiveDevice = GPU::Environment::getNumActiveDevice();
    MPI_Allgather (&my_numActiveDevice, 1, MPI_INT,
		   all_numActiveDevice, 1, MPI_INT,
		   MPI_COMM_WORLD);

    int countProcOnNode = 0;
    for (int i = 0; i < numProc_world; ++i){
      if (strcmp (my_name, &all_name[i*MPI_MAX_PROCESSOR_NAME]) == 0 &&
	  countProcOnNode < all_numActiveDevice[i]){
	all_activity[i] = GPU::Environment::cptr_activeDeviceId()[countProcOnNode++];
      }
      else {
	all_activity[i] = -1;
      }
    }

    active = (all_activity[myRank_world] != -1);
    if (isActive()){
      GPU::Environment::setDeviceId (all_activity[myRank_world]);
      // GPU::Environment::setDeviceId (1);
    }
  }
  
  {
    int my_activity = isActive();
    MPI_Allgather (&my_activity, 1, MPI_INT,
  		   all_activity, 1, MPI_INT,
  		   MPI_COMM_WORLD);
    int numActiveProc = 0;
    for (int i = 0; i < numProc_world; ++i){
      if (all_activity[i]) numActiveProc ++;
    }
    if (numActiveProc == 0){
      throw MDExcptNoActiveProcess ();
    }
    size_t tmpSize = numActiveProc * sizeof(int);
    int * rank_activeNodes = (int *) malloc (tmpSize);
    if (rank_activeNodes == NULL){
      std::cerr << "not enough space\n";
      abort (1);
    }
    int count = 0;
    for (int i = 0; i < numProc_world; ++i){
      if (all_activity[i]){
  	rank_activeNodes[count++] = i;
      }
    }
    // rank_worldRoot = rank_activeNodes[0];
    MPI_Group group_world, group_active;
    MPI_Comm_group (MPI_COMM_WORLD, &group_world);
    MPI_Group_incl (group_world, numActiveProc, rank_activeNodes, &group_active);
    MPI_Comm_create (MPI_COMM_WORLD, group_active, &commActive);
  
    free(rank_activeNodes);
  }

  if (isActive()){
    printf ("# myrank: %d, active, deviceId %d\n",
  	    myRank_world, GPU::Environment::getDeviceId());
  }
  else {
    printf ("# myrank: %d, not active\n",
  	    myRank_world);
  }
  
  free (all_activity);
  free (all_numActiveDevice);
  free (all_name);

  // int myRank_world1 = 0;
  // MPI_Comm_rank (MPI_COMM_WORLD, &myRank_world1);
  // printf ("##### myRank_world: %d is here \n", myRank_world1);
  // fflush (stdout);
  // MPI_Barrier (MPI_COMM_WORLD);
  
  if (isActive()){
    MPI_Comm_size (commActive, &numProc_);
    MPI_Comm_rank (commActive, &myRank_);
    initCart (nx, ny, nz);
    MPI_Comm_size (commCart,   &numProc_);
    MPI_Comm_rank (commCart,   &myRank_);
  }

  // printf ("@@@ myrank: %d of %d\n", myRank(), numProc());
  inited = 1;
}

void Parallel::Environment::
finalize ()
{
  // int myRank_world;
  // MPI_Comm_rank (MPI_COMM_WORLD, &myRank_world);
  // int flag;
  // if (myRank_world == rank_worldRoot){
  //   flag = 0;
  // }
  // MPI_Bcast (&flag, 1, MPI_INT, rank_worldRoot, MPI_COMM_WORLD);
  MPI_Finalize();
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
abort (int errorCode)
{
  MPI_Abort (commCart, errorCode);
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

