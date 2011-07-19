#define DEVICE_CODE

#include "MDError_interface.h"

#define NErrorIndex	2
#define NErrorScalor	13

MDError::MDError ()
{
  he = mdSuccess;
  hindex = 0;
  hscalor = 0;
  hindex = (IndexType *) malloc (sizeof(IndexType) * NErrorIndex);
  hscalor = (ScalorType *) malloc (sizeof(ScalorType) * NErrorScalor);
  for (IndexType i = 0; i < NErrorIndex; ++i){
    hindex[i] = 0;
  }
  for (IndexType i = 0; i < NErrorScalor; ++i){
    hscalor[i] = 0.f;
  }
  cudaMalloc ((void **)&ptr_de, sizeof(mdError_t));
  cudaMemcpy (ptr_de, &he, sizeof(mdError_t), cudaMemcpyHostToDevice);
  cudaMalloc ((void **)&ptr_dindex, NErrorIndex * sizeof(IndexType));
  cudaMalloc ((void **)&ptr_dscalor, NErrorScalor * sizeof(ScalorType));
  cudaMemcpy (ptr_dindex, &hindex, NErrorIndex * sizeof(IndexType),
	      cudaMemcpyHostToDevice);
  cudaMemcpy (ptr_dscalor, &hscalor, NErrorScalor * sizeof(ScalorType),
	      cudaMemcpyHostToDevice);
  // cudaMemset (ptr_dindex, 0, NErrorIndex * sizeof(IndexType));
  // cudaMemset (ptr_dscalor, 0, NErrorScalor * sizeof(ScalorType));
  checkCUDAError ("MDError::MDError");
}

MDError::~MDError()
{
  freeAPointer((void**)&hindex);
  freeAPointer((void**)&hscalor);
  cudaFree (ptr_de);
  cudaFree (ptr_dindex);
  cudaFree (ptr_dscalor);
  checkCUDAError ("MDError::~MDError");
}

inline char * MDError::getErrorString (mdError_t err) {
  switch (err){
  case mdSuccess:
      return "Success";
  case mdErrorShortCellList:
      return "The cell list is too short, increase the number of thread per block";
  case mdErrorShortNeighborList:
      return "The neighbor list is too shor, increase the DeviceNeighborListExpansion";
  case mdErrorOverFlowCellIdx:
      return "Detect an over flown cell index";
  case mdErrorBreakFENEBond:
      return "Detect a broken FENE bond";
  default:
      return "Unknow error status";
  }
}

void MDError::updateHost ()
{
  cudaMemcpy (hindex, ptr_dindex, NErrorIndex * sizeof(IndexType),
	      cudaMemcpyDeviceToHost);
  cudaMemcpy (hscalor, ptr_dscalor, NErrorScalor * sizeof(ScalorType),
	      cudaMemcpyDeviceToHost);
  checkCUDAError ("MDError::updateHost");
}


// #include "Parallel_Interface.h"
void MDError::check (const char * msg)
{
  cudaMemcpy (&he, ptr_de, sizeof(mdError_t), cudaMemcpyDeviceToHost);
  updateHost();
  if (mdSuccess != he){
//    fprintf (stderr, "myrank: %d, Md error: %s: %s.\n", Parallel::Interface::myRank(), msg, getErrorString(he));
    fprintf (stderr, "recorded indexes are");
    for (IndexType i = 0; i < NErrorIndex; ++i){
      printf ("%d  ", hindex[i]);
    }
    printf ("\n");
    for (IndexType i = 0; i < NErrorScalor; ++i){
      printf ("%f  ", hscalor[i]);
    }
    printf ("\n");
    // exit (EXIT_FAILURE);
    throw MDExcptCuda ();
  }
}

  
