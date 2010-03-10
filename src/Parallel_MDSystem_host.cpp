#define MPI_CODE

// #include "Parallel_MDSystem.h"
// #include "Parallel_TransferEngine.h"

// #include "compile_error_mixcode.h"

// void Parallel::MDSystem::
// redistribute()
// {
//   Parallel::TransferEngine tr_send;
//   Parallel::TransferEngine tr_recv;

//   int xsend0_nei, xrecv0_nei;
//   int xsend1_nei, xrecv1_nei;
//   int ysend0_nei, yrecv0_nei;
//   int ysend1_nei, yrecv1_nei;
//   int zsend0_nei, zrecv0_nei;
//   int zsend1_nei, zrecv1_nei;
  
//   Parallel::Environment env;
//   env.neighborProcIndex (Parallel::CoordXIndex,  1, xsend0_nei, xrecv0_nei);
//   env.neighborProcIndex (Parallel::CoordXIndex, -1, xsend1_nei, xrecv1_nei);
//   env.neighborProcIndex (Parallel::CoordYIndex,  1, ysend0_nei, yrecv0_nei);
//   env.neighborProcIndex (Parallel::CoordYIndex, -1, ysend1_nei, yrecv1_nei);
//   env.neighborProcIndex (Parallel::CoordZIndex,  1, zsend0_nei, zrecv0_nei);
//   env.neighborProcIndex (Parallel::CoordZIndex, -1, zsend1_nei, zrecv1_nei);

//   IndexType i = 0;

  
//   for (Parallel::SubHostCellList::iterator it = xsend0.begin();
//        it != xsend0.end(); ++it){
//     tr_send.registerBuff (&xsend0.numAtomInCell(), sizeof(IndexType));
//   }
//   tr_send.build ();
//   tr_recv.registerBuff (recv0Num, sizeof(IndexType)*xrecv0.numCell());
//   tr_recv.build ();
//   tr_send.Isend (xsend0_nei, 0);
//   tr_recv.Irecv (xrecv0_nei, 0);
//   tr_send.wait();
//   tr_send.clearRegistered();
//   tr_recv.wait();
//   tr_recv.clearRegistered();
//   for (Parallel::SubHostCellList::iterator it = xsend0.begin();
//        it != xsend0.end(); ++it){
//     DataTransferBlock block;
//     localHostData.formDataTransferBlock (it.cellStartIndex(),
//   					 it.numAtomInCell(),
//   					 block);
//     tr_send.registerBuff (block);
//   }
//   tr_send.build();
//   for (Parallel::SubHostCellList::iterator it = xrecv0.begin(), i = 0;
//        it != xrecv0.end(); ++it, ++i){
//     DataTransferBlock block;
//     localHostData.formDataTransferBlock (it.cellStartIndex() + it.numAtomInCell(),
//   					 recv0Num[i],
// 					 block);
//     it.cellStartIndex() += recv0Num[i];
//     if (it.cellStartIndex() > Parallel::Environment::cellSize()){
//       throw MDExcptTooSmallCell ("Parallel::MDSystem::redistribute");
//     }
//     tr_recv.registerBuff (block);
//   }
//   tr_recv.build ();
//   tr_send.Isend(xsend0_nei, 1);
//   tr_recv.Irecv (xrecv0_nei, 1);
//   tr_send.wait();
//   tr_send.clearRegistered();
//   tr_recv.wait();
//   tr_recv.clearRegistered();



//   for (Parallel::SubHostCellList::iterator it = xsend1.begin();
//        it != xsend1.end(); ++it){
//     tr_send.registerBuff (&xsend1.numAtomInCell(), sizeof(IndexType));
//   }
//   tr_send.build ();
//   tr_recv.registerBuff (recv1Num, sizeof(IndexType)*xrecv1.numCell());
//   tr_recv.build ();
//   tr_send.Isend (xsend1_nei, 2);
//   tr_recv.Irecv (xrecv1_nei, 2);
//   tr_send.wait();
//   tr_send.clearRegistered();
//   tr_recv.wait();
//   tr_recv.clearRegistered();
//   for (Parallel::SubHostCellList::iterator it = xsend1.begin();
//        it != xsend1.end(); ++it){
//     DataTransferBlock block;
//     localHostData.formDataTransferBlock (it.cellStartIndex(),
//   					 it.numAtomInCell(),
//   					 block);
//     tr_send.registerBuff (block);
//   }
//   tr_send.build();
//   for (Parallel::SubHostCellList::iterator it = xrecv1.begin(), i = 0;
//        it != xrecv1.end(); ++it, ++i){
//     DataTransferBlock block;
//     localHostData.formDataTransferBlock (it.cellStartIndex() + it.numAtomInCell(),
//   					 recv1Num[i],
// 					 block);
//     it.cellStartIndex() += recv1Num[i];
//     if (it.cellStartIndex() > Parallel::Environment::cellSize()){
//       throw MDExcptTooSmallCell ("Parallel::MDSystem::redistribute");
//     }
//     tr_recv.registerBuff (block);
//   }
//   tr_recv.build ();
//   tr_send.Isend (xsend1_nei, 3);
//   tr_recv.Irecv (xrecv1_nei, 3);
//   tr_send.wait();
//   tr_send.clearRegistered();
//   tr_recv.wait();
//   tr_recv.clearRegistered();

// }

	
