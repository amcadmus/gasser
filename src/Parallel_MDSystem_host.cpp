#define MPI_CODE

#include "Parallel_MDSystem.h"
#include "Parallel_TransferEngine.h"
#include "Parallel_CellList.h"
#include "Parallel_Interface.h"
#include "Parallel_Timer.h"

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

Parallel::SystemRedistributeTransferUtil::
SystemRedistributeTransferUtil ()
    : ptr_hdata (NULL), ptr_buff(NULL), mask (MDDataItemMask_All)
{
  Parallel::Interface::shiftNeighbor (CoordXIndex,  1, xsrc0, xdest0);
  Parallel::Interface::shiftNeighbor (CoordXIndex, -1, xsrc1, xdest1);
  Parallel::Interface::shiftNeighbor (CoordYIndex,  1, ysrc0, ydest0);
  Parallel::Interface::shiftNeighbor (CoordYIndex, -1, ysrc1, ydest1);
  Parallel::Interface::shiftNeighbor (CoordZIndex,  1, zsrc0, zdest0);
  Parallel::Interface::shiftNeighbor (CoordZIndex, -1, zsrc1, zdest1);
}

void Parallel::SystemRedistributeTransferUtil::
setHostData (HostCellListedMDData & hdata,
	     HostCellListedMDData & buffdata)
{  
  ptr_hdata = &hdata;
  ptr_buff  = &buffdata;
  
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
  ptr_buff ->buildSubList (devideLevel, devideLevel + devideLevel,
			   0, numCell.y,
			   0, numCell.z,
			   xrecv0);
  ptr_buff ->buildSubList (numCell.x - 2 * devideLevel, numCell.x - devideLevel,
			   0, numCell.y,
			   0, numCell.z,
			   xrecv1);
  ptr_hdata->buildSubList (devideLevel, devideLevel + devideLevel,
			   0, numCell.y,
			   0, numCell.z,
			   xrecv0h);
  ptr_hdata->buildSubList (numCell.x - 2 * devideLevel, numCell.x - devideLevel,
			   0, numCell.y,
			   0, numCell.z,
			   xrecv1h);
  ptr_hdata->buildSubList (devideLevel, numCell.x - devideLevel,
			   numCell.y - devideLevel, numCell.y,
			   0, numCell.z,
			   ysend0);
  ptr_hdata->buildSubList (devideLevel, numCell.x - devideLevel,
			   0, devideLevel,
			   0, numCell.z,
			   ysend1);
  ptr_buff ->buildSubList (devideLevel, numCell.x - devideLevel,
			   devideLevel, 2 * devideLevel,
			   0, numCell.z,
			   yrecv0);
  ptr_buff ->buildSubList (devideLevel, numCell.x - devideLevel,
			   numCell.y - 2 * devideLevel, numCell.y - devideLevel,
			   0, numCell.z,
			   yrecv1);
  ptr_hdata->buildSubList (devideLevel, numCell.x - devideLevel,
			   devideLevel, 2 * devideLevel,
			   0, numCell.z,
			   yrecv0h);
  ptr_hdata->buildSubList (devideLevel, numCell.x - devideLevel,
			   numCell.y - 2 * devideLevel, numCell.y - devideLevel,
			   0, numCell.z,
			   yrecv1h);
  ptr_hdata->buildSubList (devideLevel, numCell.x - devideLevel,
			   devideLevel, numCell.y - devideLevel,
			   numCell.z - devideLevel, numCell.z,
			   zsend0);
  ptr_hdata->buildSubList (devideLevel, numCell.x - devideLevel,
			   devideLevel, numCell.y - devideLevel,
			   0, devideLevel,
			   zsend1);
  ptr_buff ->buildSubList (devideLevel, numCell.x - devideLevel,
			   devideLevel, numCell.y - devideLevel,
			   devideLevel, 2 * devideLevel,
			   zrecv0);
  ptr_buff ->buildSubList (devideLevel, numCell.x - devideLevel,
			   devideLevel, numCell.y - devideLevel,
			   numCell.z - 2 * devideLevel, numCell.z - devideLevel,
			   zrecv1);
  ptr_hdata->buildSubList (devideLevel, numCell.x - devideLevel,
			   devideLevel, numCell.y - devideLevel,
			   devideLevel, 2 * devideLevel,
			   zrecv0h);
  ptr_hdata->buildSubList (devideLevel, numCell.x - devideLevel,
			   devideLevel, numCell.y - devideLevel,
			   numCell.z - 2 * devideLevel, numCell.z - devideLevel,
			   zrecv1h);

  // xsend0.mallocNumAtomSend (&xNumSend0, &xNumSend0size);
  // xrecv0.mallocNumAtomSend (&xNumRecv0, &xNumRecv0size);
  // xsend1.mallocNumAtomSend (&xNumSend1, &xNumSend1size);
  // xrecv1.mallocNumAtomSend (&xNumRecv1, &xNumRecv1size);
  // ysend0.mallocNumAtomSend (&yNumSend0, &yNumSend0size);
  // yrecv0.mallocNumAtomSend (&yNumRecv0, &yNumRecv0size);
  // ysend1.mallocNumAtomSend (&yNumSend1, &yNumSend1size);
  // yrecv1.mallocNumAtomSend (&yNumRecv1, &yNumRecv1size);
  // zsend0.mallocNumAtomSend (&zNumSend0, &zNumSend0size);
  // zrecv0.mallocNumAtomSend (&zNumRecv0, &zNumRecv0size);
  // zsend1.mallocNumAtomSend (&zNumSend1, &zNumSend1size);
  // zrecv1.mallocNumAtomSend (&zNumRecv1, &zNumRecv1size);  

  xsendNum0.reinit (xsend0);
  xrecvNum0.reinit (xrecv0);
  xsendNum1.reinit (xsend1);
  xrecvNum1.reinit (xrecv1);
  ysendNum0.reinit (ysend0);
  yrecvNum0.reinit (yrecv0);
  ysendNum1.reinit (ysend1);
  yrecvNum1.reinit (yrecv1);
  zsendNum0.reinit (zsend0);
  zrecvNum0.reinit (zrecv0);
  zsendNum1.reinit (zsend1);
  zrecvNum1.reinit (zrecv1);

  xsendData0.reinit (xsend0, mask);
  xrecvData0.reinit (xrecv0, mask);
  xsendData1.reinit (xsend1, mask);
  xrecvData1.reinit (xrecv1, mask);
  ysendData0.reinit (ysend0, mask);
  yrecvData0.reinit (yrecv0, mask);
  ysendData1.reinit (ysend1, mask);
  yrecvData1.reinit (yrecv1, mask);
  zsendData0.reinit (zsend0, mask);
  zrecvData0.reinit (zrecv0, mask);
  zsendData1.reinit (zsend1, mask);
  zrecvData1.reinit (zrecv1, mask);
  
}

using namespace Parallel::Timer;

void Parallel::SystemRedistributeTransferUtil::
redistributeHost ()
{
  Parallel::TransferEngine sender;
  Parallel::TransferEngine recver;
  
  // ptr_buff->clearData();

  HostTimer::tic (item_Redistribute_Data0);
  sender.clearRegistered();
  recver.clearRegistered();
  sender.registerBuff (xsendNum0);
  recver.registerBuff (xrecvNum0);
  sender.build ();
  recver.build ();
  HostTimer::toc (item_Redistribute_Data0);
  HostTimer::tic (item_Redistribute_Transfer);
  sender.Isend (xdest0, 0);
  recver.Irecv (xsrc0,  0);
  sender.wait();
  recver.wait();
  HostTimer::toc (item_Redistribute_Transfer);

  xsendData0.build ();
  xrecvData0.build ();
  HostTimer::tic (item_Redistribute_Data);
  sender.clearRegistered();
  recver.clearRegistered();
  sender.registerBuff (xsendData0);
  recver.registerBuff (xrecvData0);
  sender.build ();
  recver.build ();
  HostTimer::toc (item_Redistribute_Data);
  HostTimer::tic (item_Redistribute_Transfer);
  sender.Isend (xdest0, 1);
  recver.Irecv (xsrc0,  1);
  sender.wait();
  recver.wait();
  xsend0.clearData();
  xrecv0h.add (xrecv0, mask);
  HostTimer::toc (item_Redistribute_Transfer);

  HostTimer::tic (item_Redistribute_Data0);
  sender.clearRegistered();
  recver.clearRegistered();
  sender.registerBuff (xsendNum1);
  recver.registerBuff (xrecvNum1);
  sender.build ();
  recver.build ();
  HostTimer::toc (item_Redistribute_Data0);
  HostTimer::tic (item_Redistribute_Transfer);
  sender.Isend (xdest1, 2);
  recver.Irecv (xsrc1,  2);
  sender.wait();
  recver.wait();  
  HostTimer::toc (item_Redistribute_Transfer);

  xsendData1.build ();
  xrecvData1.build ();
  HostTimer::tic (item_Redistribute_Data);
  sender.clearRegistered();
  recver.clearRegistered();
  sender.registerBuff (xsendData1);
  recver.registerBuff (xrecvData1);
  sender.build ();
  recver.build ();
  HostTimer::toc (item_Redistribute_Data);
  HostTimer::tic (item_Redistribute_Transfer);
  sender.Isend (xdest1, 3);
  recver.Irecv (xsrc1,  3);
  sender.wait();
  recver.wait();
  xsend1.clearData();  
  xrecv1h.add (xrecv1, mask);
  HostTimer::toc (item_Redistribute_Transfer);

  Parallel::Interface::barrier();

  HostTimer::tic (item_Redistribute_Data0);
  sender.clearRegistered();
  recver.clearRegistered();
  sender.registerBuff (ysendNum0);
  recver.registerBuff (yrecvNum0);
  sender.build ();
  recver.build ();
  HostTimer::toc (item_Redistribute_Data0);
  HostTimer::tic (item_Redistribute_Transfer);
  sender.Isend (ydest0, 0);
  recver.Irecv (ysrc0,  0);
  sender.wait();
  recver.wait();  
  HostTimer::toc (item_Redistribute_Transfer);

  ysendData0.build ();
  yrecvData0.build ();
  HostTimer::tic (item_Redistribute_Data);
  sender.clearRegistered();
  recver.clearRegistered();
  sender.registerBuff (ysendData0);
  recver.registerBuff (yrecvData0);
  sender.build ();
  recver.build ();
  HostTimer::toc (item_Redistribute_Data);
  HostTimer::tic (item_Redistribute_Transfer);
  sender.Isend (ydest0, 1);
  recver.Irecv (ysrc0,  1);
  sender.wait();
  recver.wait();
  ysend0.clearData();
  yrecv0h.add (yrecv0, mask);
  HostTimer::toc (item_Redistribute_Transfer);

  HostTimer::tic (item_Redistribute_Data0);
  sender.clearRegistered();
  recver.clearRegistered();
  sender.registerBuff (ysendNum1);
  recver.registerBuff (yrecvNum1);
  sender.build ();
  recver.build ();
  HostTimer::toc (item_Redistribute_Data0);
  HostTimer::tic (item_Redistribute_Transfer);
  sender.Isend (ydest1, 2);
  recver.Irecv (ysrc1,  2);
  sender.wait();
  recver.wait();  
  HostTimer::toc (item_Redistribute_Transfer);

  ysendData1.build ();
  yrecvData1.build ();
  HostTimer::tic (item_Redistribute_Data);
  sender.clearRegistered();
  recver.clearRegistered();
  sender.registerBuff (ysendData1);
  recver.registerBuff (yrecvData1);
  sender.build ();
  recver.build ();
  HostTimer::toc (item_Redistribute_Data);
  HostTimer::tic (item_Redistribute_Transfer);
  sender.Isend (ydest1, 3);
  recver.Irecv (ysrc1,  3);
  sender.wait();
  recver.wait();
  ysend1.clearData();
  yrecv1h.add (yrecv1, mask);
  HostTimer::toc (item_Redistribute_Transfer);

  Parallel::Interface::barrier();
  
  HostTimer::tic (item_Redistribute_Data0);
  sender.clearRegistered();
  recver.clearRegistered();
  sender.registerBuff (zsendNum0);
  recver.registerBuff (zrecvNum0);
  sender.build ();
  recver.build ();
  HostTimer::toc (item_Redistribute_Data0);
  HostTimer::tic (item_Redistribute_Transfer);
  sender.Isend (zdest0, 0);
  recver.Irecv (zsrc0,  0);
  sender.wait();
  recver.wait();  
  HostTimer::toc (item_Redistribute_Transfer);

  zsendData0.build ();
  zrecvData0.build ();
  HostTimer::tic (item_Redistribute_Data);
  sender.clearRegistered();
  recver.clearRegistered();
  sender.registerBuff (zsendData0);
  recver.registerBuff (zrecvData0);
  sender.build ();
  recver.build ();
  HostTimer::toc (item_Redistribute_Data);
  HostTimer::tic (item_Redistribute_Transfer);
  sender.Isend (zdest0, 1);
  recver.Irecv (zsrc0,  1);
  sender.wait();
  recver.wait();
  zsend0.clearData ();
  zrecv0h.add (zrecv0, mask);
  HostTimer::toc (item_Redistribute_Transfer);

  HostTimer::tic (item_Redistribute_Data0);
  sender.clearRegistered();
  recver.clearRegistered();
  sender.registerBuff (zsendNum1);
  recver.registerBuff (zrecvNum1);
  sender.build ();
  recver.build ();
  HostTimer::toc (item_Redistribute_Data0);
  HostTimer::tic (item_Redistribute_Transfer);
  sender.Isend (zdest1, 2);
  recver.Irecv (zsrc1,  2);
  sender.wait();
  recver.wait();  
  HostTimer::toc (item_Redistribute_Transfer);

  zsendData1.build ();
  zrecvData1.build ();
  HostTimer::tic (item_Redistribute_Data);
  sender.clearRegistered();
  recver.clearRegistered();
  sender.registerBuff (zsendData1);
  recver.registerBuff (zrecvData1);
  sender.build ();
  recver.build ();
  HostTimer::toc (item_Redistribute_Data);
  HostTimer::tic (item_Redistribute_Transfer);
  sender.Isend (zdest1, 3);
  recver.Irecv (zsrc1,  3);
  sender.wait();
  recver.wait();
  zsend1.clearData();
  zrecv1h.add (zrecv1, mask);
  HostTimer::toc (item_Redistribute_Transfer);
}



Parallel::SystemCollectDataUtil::
SystemCollectDataUtil()
    : ptr_hdata(NULL), ptr_buff(NULL),
      ptr_gdata(NULL),
      mask (MDDataItemMask_Coordinate |
	    MDDataItemMask_CoordinateNoi |
	    MDDataItemMask_Velocity |
	    MDDataItemMask_Force |
	    MDDataItemMask_GlobalIndex)
{
}

void Parallel::SystemCollectDataUtil::
addBuffToGlobal ()
{
  IndexType totalNumCell = ptr_buff->getNumCell().x *
      ptr_buff->getNumCell().y * ptr_buff->getNumCell().z;
  IndexType numThreadsInCell = Parallel::Interface::numThreadsInCell();
  
  for (IndexType i = 0; i < totalNumCell; ++i){
    IndexType cellShift = i * numThreadsInCell;
    if (mask & MDDataItemMask_Coordinate){
      for (IndexType j = cellShift;
	   j < cellShift + ptr_buff->cptr_numAtomInCell()[i];
	   ++j){
	IndexType gIndex = ptr_buff->cptr_globalIndex()[j];
	ptr_gdata->cptr_coordinate()[gIndex] = ptr_buff->cptr_coordinate()[j];
      }
    }

    if (mask & MDDataItemMask_CoordinateNoi){
      for (IndexType j = cellShift;
	   j < cellShift + ptr_buff->cptr_numAtomInCell()[i];
	   ++j){
	IndexType gIndex = ptr_buff->cptr_globalIndex()[j];
	ptr_gdata->cptr_coordinateNoi()[gIndex].x = ptr_buff->cptr_coordinateNoi()[j].x;
	ptr_gdata->cptr_coordinateNoi()[gIndex].y = ptr_buff->cptr_coordinateNoi()[j].y;
	ptr_gdata->cptr_coordinateNoi()[gIndex].z = ptr_buff->cptr_coordinateNoi()[j].z;
      }
    }

    if (mask & MDDataItemMask_Velocity){
      for (IndexType j = cellShift;
	   j < cellShift + ptr_buff->cptr_numAtomInCell()[i];
	   ++j){
	IndexType gIndex = ptr_buff->cptr_globalIndex()[j];
	ptr_gdata->cptr_velocityX()[gIndex] = ptr_buff->cptr_velocityX()[j];
	ptr_gdata->cptr_velocityY()[gIndex] = ptr_buff->cptr_velocityY()[j];
	ptr_gdata->cptr_velocityZ()[gIndex] = ptr_buff->cptr_velocityZ()[j];
      }
    }

    if (mask & MDDataItemMask_Force){
      for (IndexType j = cellShift;
	   j < cellShift + ptr_buff->cptr_numAtomInCell()[i];
	   ++j){
	IndexType gIndex = ptr_buff->cptr_globalIndex()[j];
	ptr_gdata->cptr_forceX()[gIndex] = ptr_buff->cptr_forceX()[j];
	ptr_gdata->cptr_forceY()[gIndex] = ptr_buff->cptr_forceY()[j];
	ptr_gdata->cptr_forceZ()[gIndex] = ptr_buff->cptr_forceZ()[j];
      }
    }
  }
}

      
      

void Parallel::SystemCollectDataUtil::
setHostData (HostCellListedMDData & hdata,
	     HostCellListedMDData & buffdata,
	     HostMDData & globalData)
{
  ptr_hdata = &hdata;
  ptr_buff = &buffdata;
  ptr_gdata = &globalData;

  int myRank = Parallel::Interface::myRank();
  
  if (myRank == 0){
    ptr_buff->buildSubListRealCell (recvlist);
    recvNum.reinit (recvlist);
    recvData.reinit (recvlist, mask);
  }
  else {
    ptr_hdata->buildSubListRealCell (sendlist);
    sendNum.reinit (sendlist);
    sendData.reinit (sendlist, mask);
  }
}


void Parallel::SystemCollectDataUtil::
collect ()
{
  int myRank = Parallel::Interface::myRank();
  if (myRank == 0){
    ptr_buff->copy (*ptr_hdata, mask);
    addBuffToGlobal ();
  }
  Parallel::Interface::barrier();
  
  int Nx, Ny, Nz;
  Parallel::Interface::numProcDim (Nx, Ny, Nz);

  TransferEngine sender;
  TransferEngine recver;
  
  for (int ix = 0; ix < Nx; ++ix){
    for (int iy = 0; iy < Ny; ++iy){
      for (int iz = 0; iz < Nz; ++iz){
	int workingRank;
	Parallel::Interface::cartCoordToRank (ix, iy, iz, workingRank);
	if (workingRank == 0) continue;

	if (myRank == workingRank){
	  sender.clearRegistered();
	  sender.registerBuff (sendNum);
	  sender.build ();
	  sender.Isend (0, workingRank * 10 + 1);
	  sender.wait ();
	  sendData.build ();
	  sender.clearRegistered ();
	  sender.registerBuff (sendData);
	  sender.build ();
	  sender.Isend (0, workingRank * 10 + 2);
	  sender.wait ();
	}
	else if (myRank == 0){
	  recver.clearRegistered();
	  recver.registerBuff (recvNum);
	  recver.build ();
	  recver.Irecv (workingRank, workingRank*10 + 1);
	  recver.wait ();
	  recvData.build ();
	  recver.clearRegistered ();
	  recver.registerBuff (recvData);
	  recver.build ();
	  recver.Irecv (workingRank, workingRank*10 + 2);
	  recver.wait();
	  addBuffToGlobal ();
	}
	Parallel::Interface::barrier();
      }
    }
  }
}



Parallel::SystemTransCoordsTransferUtil::
SystemTransCoordsTransferUtil ()
    : ptr_hdata (NULL), ptr_buff(NULL),
      mask (MDDataItemMask_Coordinate |
	    MDDataItemMask_GlobalIndex)
{
  Parallel::Interface::shiftNeighbor (CoordXIndex,  1, xsrc0, xdest0);
  Parallel::Interface::shiftNeighbor (CoordXIndex, -1, xsrc1, xdest1);
  Parallel::Interface::shiftNeighbor (CoordYIndex,  1, ysrc0, ydest0);
  Parallel::Interface::shiftNeighbor (CoordYIndex, -1, ysrc1, ydest1);
  Parallel::Interface::shiftNeighbor (CoordZIndex,  1, zsrc0, zdest0);
  Parallel::Interface::shiftNeighbor (CoordZIndex, -1, zsrc1, zdest1);
}

void Parallel::SystemTransCoordsTransferUtil::
setHostData (HostCellListedMDData & hdata,
	     HostCellListedMDData & buffdata)
{  
  ptr_hdata = &hdata;
  ptr_buff  = &buffdata;
  
  HostIntVectorType numCell = ptr_hdata->getNumCell ();
  IndexType devideLevel = ptr_hdata->getDevideLevel ();

  ptr_hdata->buildSubList (numCell.x - 2 * devideLevel, numCell.x - devideLevel,
			   devideLevel, numCell.y - devideLevel,
			   devideLevel, numCell.z - devideLevel,
			   xsend0);
  ptr_hdata->buildSubList (devideLevel, 2 * devideLevel,
			   devideLevel, numCell.y - devideLevel,
			   devideLevel, numCell.z - devideLevel,
			   xsend1);
  ptr_hdata->buildSubList (0, devideLevel,
			   devideLevel, numCell.y - devideLevel,
			   devideLevel, numCell.z - devideLevel,
			   xrecv0);
  ptr_hdata->buildSubList (numCell.x - devideLevel, numCell.x,
			   devideLevel, numCell.y - devideLevel,
			   devideLevel, numCell.z - devideLevel,
			   xrecv1);
  
  ptr_hdata->buildSubList (0, numCell.x,
			   numCell.y - 2 * devideLevel, numCell.y - devideLevel,
			   devideLevel, numCell.z - devideLevel,
			   ysend0);
  ptr_hdata->buildSubList (0, numCell.x,
			   devideLevel, 2 * devideLevel,
			   devideLevel, numCell.z - devideLevel,
			   ysend1);
  ptr_hdata->buildSubList (0, numCell.x,
			   0, devideLevel,
			   devideLevel, numCell.z - devideLevel,
			   yrecv0);
  ptr_hdata->buildSubList (0, numCell.x,
			   numCell.y - devideLevel, numCell.y,
			   devideLevel, numCell.z - devideLevel,
			   yrecv1);
  
  ptr_hdata->buildSubList (0, numCell.x,
			   0, numCell.y,
			   numCell.z - 2 * devideLevel, numCell.z - devideLevel,
			   zsend0);
  ptr_hdata->buildSubList (0, numCell.x,
			   0, numCell.y,
			   devideLevel, devideLevel * 2,
			   zsend1);
  ptr_hdata->buildSubList (0, numCell.x,
			   0, numCell.y,
			   0, devideLevel,
			   zrecv0);
  ptr_hdata->buildSubList (0, numCell.x,
			   0, numCell.y,
			   numCell.z - devideLevel, numCell.z,
			   zrecv1);

  xsendNum0.reinit (xsend0);
  xrecvNum0.reinit (xrecv0);
  xsendNum1.reinit (xsend1);
  xrecvNum1.reinit (xrecv1);
  ysendNum0.reinit (ysend0);
  yrecvNum0.reinit (yrecv0);
  ysendNum1.reinit (ysend1);
  yrecvNum1.reinit (yrecv1);
  zsendNum0.reinit (zsend0);
  zrecvNum0.reinit (zrecv0);
  zsendNum1.reinit (zsend1);
  zrecvNum1.reinit (zrecv1);

  xsendData0.reinit (xsend0, mask);
  xrecvData0.reinit (xrecv0, mask);
  xsendData1.reinit (xsend1, mask);
  xrecvData1.reinit (xrecv1, mask);
  ysendData0.reinit (ysend0, mask);
  yrecvData0.reinit (yrecv0, mask);
  ysendData1.reinit (ysend1, mask);
  yrecvData1.reinit (yrecv1, mask);
  zsendData0.reinit (zsend0, mask);
  zrecvData0.reinit (zrecv0, mask);
  zsendData1.reinit (zsend1, mask);
  zrecvData1.reinit (zrecv1, mask);
  
}


void Parallel::SystemTransCoordsTransferUtil::
transCoords ()
{
  Parallel::TransferEngine sender;
  Parallel::TransferEngine recver;
  
  // ptr_buff->clearData();

  sender.clearRegistered();
  recver.clearRegistered();
  sender.registerBuff (xsendNum0);
  recver.registerBuff (xrecvNum0);
  sender.build ();
  recver.build ();
  sender.Isend (xdest0, 0);
  recver.Irecv (xsrc0,  0);
  sender.wait();
  recver.wait();

  xsendData0.build ();
  xrecvData0.build ();
  sender.clearRegistered();
  recver.clearRegistered();
  sender.registerBuff (xsendData0);
  recver.registerBuff (xrecvData0);
  sender.build ();
  recver.build ();
  sender.Isend (xdest0, 1);
  recver.Irecv (xsrc0,  1);
  sender.wait();
  recver.wait();

  sender.clearRegistered();
  recver.clearRegistered();
  sender.registerBuff (xsendNum1);
  recver.registerBuff (xrecvNum1);
  sender.build ();
  recver.build ();
  sender.Isend (xdest1, 2);
  recver.Irecv (xsrc1,  2);
  sender.wait();
  recver.wait();  

  xsendData1.build ();
  xrecvData1.build ();
  sender.clearRegistered();
  recver.clearRegistered();
  sender.registerBuff (xsendData1);
  recver.registerBuff (xrecvData1);
  sender.build ();
  recver.build ();
  sender.Isend (xdest1, 3);
  recver.Irecv (xsrc1,  3);
  sender.wait();
  recver.wait();

  Parallel::Interface::barrier();

  sender.clearRegistered();
  recver.clearRegistered();
  sender.registerBuff (ysendNum0);
  recver.registerBuff (yrecvNum0);
  sender.build ();
  recver.build ();
  sender.Isend (ydest0, 0);
  recver.Irecv (ysrc0,  0);
  sender.wait();
  recver.wait();  

  ysendData0.build ();
  yrecvData0.build ();
  sender.clearRegistered();
  recver.clearRegistered();
  sender.registerBuff (ysendData0);
  recver.registerBuff (yrecvData0);
  sender.build ();
  recver.build ();
  sender.Isend (ydest0, 1);
  recver.Irecv (ysrc0,  1);
  sender.wait();
  recver.wait();

  sender.clearRegistered();
  recver.clearRegistered();
  sender.registerBuff (ysendNum1);
  recver.registerBuff (yrecvNum1);
  sender.build ();
  recver.build ();
  sender.Isend (ydest1, 2);
  recver.Irecv (ysrc1,  2);
  sender.wait();
  recver.wait();  

  ysendData1.build ();
  yrecvData1.build ();
  sender.clearRegistered();
  recver.clearRegistered();
  sender.registerBuff (ysendData1);
  recver.registerBuff (yrecvData1);
  sender.build ();
  recver.build ();
  sender.Isend (ydest1, 3);
  recver.Irecv (ysrc1,  3);
  sender.wait();
  recver.wait();

  Parallel::Interface::barrier();
  
  sender.clearRegistered();
  recver.clearRegistered();
  sender.registerBuff (zsendNum0);
  recver.registerBuff (zrecvNum0);
  sender.build ();
  recver.build ();
  sender.Isend (zdest0, 0);
  recver.Irecv (zsrc0,  0);
  sender.wait();
  recver.wait();  

  zsendData0.build ();
  zrecvData0.build ();
  sender.clearRegistered();
  recver.clearRegistered();
  sender.registerBuff (zsendData0);
  recver.registerBuff (zrecvData0);
  sender.build ();
  recver.build ();
  sender.Isend (zdest0, 1);
  recver.Irecv (zsrc0,  1);
  sender.wait();
  recver.wait();

  sender.clearRegistered();
  recver.clearRegistered();
  sender.registerBuff (zsendNum1);
  recver.registerBuff (zrecvNum1);
  sender.build ();
  recver.build ();
  sender.Isend (zdest1, 2);
  recver.Irecv (zsrc1,  2);
  sender.wait();
  recver.wait();  

  zsendData1.build ();
  zrecvData1.build ();
  sender.clearRegistered();
  recver.clearRegistered();
  sender.registerBuff (zsendData1);
  recver.registerBuff (zrecvData1);
  sender.build ();
  recver.build ();
  sender.Isend (zdest1, 3);
  recver.Irecv (zsrc1,  3);
  sender.wait();
  recver.wait();

  Parallel::Interface::barrier();
}
