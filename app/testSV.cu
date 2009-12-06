#include "SumVector.h"

// #define N 640000
#define N 640

int main(int argc, char * argv[])
{
  
  ScalorType hdata [N];
  for (IndexType i = 0; i < N; ++i){
    hdata[i] = i;
  }

  IndexType Nt = 32;
  if (argc != 1){
    Nt = atoi(argv[1]);
    printf("Nthread  is %d\n", Nt);
  }
  
  SumVector<ScalorType> sv;
  sv.init (N, Nt);
  cudaMemcpy (sv.getBuff(), hdata, sizeof(ScalorType) * N, cudaMemcpyHostToDevice);
  
  ScalorType * dresult;
  cudaMalloc ((void**)&dresult, sizeof(ScalorType));
  ScalorType hresule;


  cudaEvent_t start, stop;
  float tmptime;
  cudaEventCreate (&start);
  cudaEventCreate (&stop );
  cudaEventRecord(start, 0);
  // for (unsigned i = 0; i < 10000; ++i){
  for (unsigned i = 0; i < 1; ++i){
    sv.sumBuff (dresult, 0);
  }
  cudaEventRecord(stop, 0);
  cudaEventSynchronize (stop);
  cudaEventElapsedTime (&tmptime, start, stop);
  
  cudaMemcpy (&hresule, dresult, sizeof(ScalorType), cudaMemcpyDeviceToHost);
  printf ("result: %f, time %.3e\n", hresule, tmptime);

  
  cudaFree(dresult);
  cudaEventDestroy (start);
  cudaEventDestroy (stop);
  
  return 0;
}
