// #define CPP_FILE

// #include "common.h"
// #include "Parallel_Environment.h"
// #include "Parallel_TransferEngine.h"
#include <mpi.h>

#define NUM 5

int main(int argc, char * argv[])
{

	
	
// Parallel::Environment env (&argc, &argv);

  // int dim[3];
  // dim[0] = 2;
  // dim[1] = dim[2] = 1;

  // env.init (dim);

  // int a;
  // double b[NUM];

  // int c;
  // double d[NUM];
  
  // Parallel::TransferEngine transS (env);
  // Parallel::TransferEngine transR (env);

  // for (unsigned i = 0; i < env.numProc(); ++i){
  //   if (env.myRank() == 0){
  //     a = i+1;
  //     for (unsigned j = 0; j < NUM; ++j){
  // 	b[j] = i+5;
  //     }
  //     transS.clearRegistered ();
  //     transS.registerBuff (&a, sizeof(int));
  //     transS.registerBuff (b,  sizeof(double) * NUM);
  //     transS.build ();
  //     transS.Isend (i, i);
  //   }
  //   if (env.myRank() == i){
  //     transR.clearRegistered ();
  //     transR.registerBuff (&c, sizeof(int));
  //     transR.registerBuff (d,  sizeof(double) * NUM);
  //     transR.build ();
  //     transR.Irecv (0, i);
  //     transR.wait ();
  //   }
  //   if (env.myRank() == 0){
  //     transS.wait ();
  //   }
  // }

  // printf ("myrank is %d, c is %d, d is %f", env.myRank(), c, d[0]);
  // fflush (stdout);
  // for (unsigned i = 1; i < NUM; i ++){
  //   printf (" %f", d[i]);
  //   fflush (stdout);
  // }
  // printf ("\n");
  // fflush (stdout);
  



      
  return 0;
}

