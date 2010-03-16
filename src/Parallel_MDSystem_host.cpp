#define MPI_CODE

#include "Parallel_MDSystem.h"
#include "Parallel_TransferEngine.h"
#include "Parallel_CellList.h"
#include "Parallel_Interface.h"

#include "compile_error_mixcode.h"

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

	

// void Parallel::MDSystem::
// tryHostSend ()
// {
//   HostSubCellList sublist0;
//   HostSubCellList sublist1;

//   localHostData.buildSubListAllCell (sublist0);
//   localHostData1.buildSubListAllCell (sublist1);
  
//   Parallel::TransferEngine tran0;
//   Parallel::TransferEngine tran1;
  
//   tran0.registerBuff (sublist0, MDDataItemMask_All);
//   tran1.registerBuff (sublist1, MDDataItemMask_All);
  
//   tran0.Isend (0, 0);
//   tran1.Irecv (0, 0);

//   tran0.wait();
//   tran1.wait();

//   return;
// }

Parallel::SystemTranferUtils::
SystemTranferUtils ()
    : ptr_hdata (NULL)
{
  Parallel::Interface::shiftNeighbor (CoordXIndex,  1, xsrc0, xdest0);
  Parallel::Interface::shiftNeighbor (CoordXIndex, -1, xsrc1, xdest1);
  Parallel::Interface::shiftNeighbor (CoordYIndex,  1, ysrc0, ydest0);
  Parallel::Interface::shiftNeighbor (CoordYIndex, -1, ysrc1, ydest1);
  Parallel::Interface::shiftNeighbor (CoordZIndex,  1, zsrc0, zdest0);
  Parallel::Interface::shiftNeighbor (CoordZIndex, -1, zsrc1, zdest1);
}

void Parallel::SystemTranferUtils::
setHostData (HostCellListedMDData & hdata)
{
  ptr_hdata = &hdata;

  HostIntVectorType numCell = ptr_hdata->getNumCell ();
  IndexType devideLevel = ptr_hdata->getDevideLevel ();

  ptr_hdata->buildSubList (numCell.x - devideLevel, numCell.x,
			   0, numCell.y,
			   0, numCell.z,
			   xsend0);
  ptr_hdata->buildSubList (0, devideLevel,
			   0, numCell.y,
			   0, numCell.z,
			   xsend1);
  ptr_hdata->buildSubList (devideLevel, devideLevel + devideLevel,
			   0, numCell.y,
			   0, numCell.z,
			   xrecv0);
  ptr_hdata->buildSubList (numCell.x - 2 * devideLevel, numCell.x - devideLevel,
			   0, numCell.y,
			   0, numCell.z,
			   xrecv1);
  ptr_hdata->buildSubList (devideLevel, numCell.x - devideLevel,
			   numCell.y - devideLevel, numCell.y,
			   0, numCell.z,
			   ysend0);
  ptr_hdata->buildSubList (devideLevel, numCell.x - devideLevel,
			   0, devideLevel,
			   0, numCell.z,
			   ysend1);
  ptr_hdata->buildSubList (devideLevel, numCell.x - devideLevel,
			   devideLevel, 2 * devideLevel,
			   0, numCell.z,
			   yrecv0);
  ptr_hdata->buildSubList (devideLevel, numCell.x - devideLevel,
			   numCell.y - 2 * devideLevel, numCell.y - devideLevel,
			   0, numCell.z,
			   yrecv1);
  ptr_hdata->buildSubList (devideLevel, numCell.x - devideLevel,
			   devideLevel, numCell.y - devideLevel,
			   numCell.z - devideLevel, numCell.z,
			   zsend0);
  ptr_hdata->buildSubList (devideLevel, numCell.x - devideLevel,
			   devideLevel, numCell.y - devideLevel,
			   0, devideLevel,
			   zsend1);
  ptr_hdata->buildSubList (devideLevel, numCell.x - devideLevel,
			   devideLevel, numCell.y - devideLevel,
			   devideLevel, 2 * devideLevel,
			   zrecv0);
  ptr_hdata->buildSubList (devideLevel, numCell.x - devideLevel,
			   devideLevel, numCell.y - devideLevel,
			   numCell.z - 2 * devideLevel, numCell.z - devideLevel,
			   zrecv1);
}


void Parallel::SystemTranferUtils::
redistribute ()
{
  MDDataItemMask_t mask = MDDataItemMask_AllExceptForce;

  Parallel::TransferEngine sender;
  Parallel::TransferEngine recver;

  sender.clearRegistered();
  recver.clearRegistered();
  sender.registerBuff (xsend0, mask);
  recver.registerBuff (xrecv0, mask);
  sender.Isend (xdest0, 0);
  recver.Irecv (xsrc0, 0);
  sender.wait();
  recver.wait();  
  sender.clearRegistered();
  recver.clearRegistered();
  sender.registerBuff (xsend1, mask);
  recver.registerBuff (xrecv1, mask);
  sender.Isend (xdest1, 1);
  recver.Irecv (xsrc1, 1);
  sender.wait();
  recver.wait();


  sender.clearRegistered();
  recver.clearRegistered();
  sender.registerBuff (ysend0, mask);
  recver.registerBuff (yrecv0, mask);
  sender.Isend (ydest0, 0);
  recver.Irecv (ysrc0, 0);
  sender.wait();
  recver.wait();  
  sender.clearRegistered();
  recver.clearRegistered();
  sender.registerBuff (ysend1, mask);
  recver.registerBuff (yrecv1, mask);
  sender.Isend (ydest1, 1);
  recver.Irecv (ysrc1, 1);
  sender.wait();
  recver.wait();

  sender.clearRegistered();
  recver.clearRegistered();
  sender.registerBuff (zsend0, mask);
  recver.registerBuff (zrecv0, mask);
  sender.Isend (zdest0, 0);
  recver.Irecv (zsrc0, 0);
  sender.wait();
  recver.wait();  
  sender.clearRegistered();
  recver.clearRegistered();
  sender.registerBuff (zsend1, mask);
  recver.registerBuff (zrecv1, mask);
  sender.Isend (zdest1, 1);
  recver.Irecv (zsrc1, 1);
  sender.wait();
  recver.wait();
}

