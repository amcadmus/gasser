#define MPI_CODE

#include "Parallel_MDSystem.h"
#include "Parallel_TransferEngine.h"
#include "Parallel_CellList.h"
#include "Parallel_Interface.h"
#include "Parallel_Timer.h"

#include "compile_error_mixcode.h"

namespace SystemRedistributeTransferUtil_globalEngine{
  Parallel::TransferEngine engine_xsend0;
  Parallel::TransferEngine engine_xrecv0;
  Parallel::TransferEngine engine_xsend1;
  Parallel::TransferEngine engine_xrecv1;
  Parallel::TransferEngine engine_ysend0;
  Parallel::TransferEngine engine_yrecv0;
  Parallel::TransferEngine engine_ysend1;
  Parallel::TransferEngine engine_yrecv1;
  Parallel::TransferEngine engine_zsend0;
  Parallel::TransferEngine engine_zrecv0;
  Parallel::TransferEngine engine_zsend1;
  Parallel::TransferEngine engine_zrecv1;

  Parallel::TransferEngine thisNum_engine_xsend0;
  Parallel::TransferEngine thisNum_engine_xrecv0;
  Parallel::TransferEngine thisNum_engine_xsend1;
  Parallel::TransferEngine thisNum_engine_xrecv1;
  Parallel::TransferEngine thisNum_engine_ysend0;
  Parallel::TransferEngine thisNum_engine_yrecv0;
  Parallel::TransferEngine thisNum_engine_ysend1;
  Parallel::TransferEngine thisNum_engine_yrecv1;
  Parallel::TransferEngine thisNum_engine_zsend0;
  Parallel::TransferEngine thisNum_engine_zrecv0;
  Parallel::TransferEngine thisNum_engine_zsend1;
  Parallel::TransferEngine thisNum_engine_zrecv1;
}

namespace SystemTransCoordsTransferUtil_globalEngine{
  Parallel::TransferEngine engine_xsend0;
  Parallel::TransferEngine engine_xrecv0;
  Parallel::TransferEngine engine_xsend1;
  Parallel::TransferEngine engine_xrecv1;
  Parallel::TransferEngine engine_ysend0;
  Parallel::TransferEngine engine_yrecv0;
  Parallel::TransferEngine engine_ysend1;
  Parallel::TransferEngine engine_yrecv1;
  Parallel::TransferEngine engine_zsend0;
  Parallel::TransferEngine engine_zrecv0;
  Parallel::TransferEngine engine_zsend1;
  Parallel::TransferEngine engine_zrecv1;

  Parallel::TransferEngine thisNum_engine_xsend0;
  Parallel::TransferEngine thisNum_engine_xrecv0;
  Parallel::TransferEngine thisNum_engine_xsend1;
  Parallel::TransferEngine thisNum_engine_xrecv1;
  Parallel::TransferEngine thisNum_engine_ysend0;
  Parallel::TransferEngine thisNum_engine_yrecv0;
  Parallel::TransferEngine thisNum_engine_ysend1;
  Parallel::TransferEngine thisNum_engine_yrecv1;
  Parallel::TransferEngine thisNum_engine_zsend0;
  Parallel::TransferEngine thisNum_engine_zrecv0;
  Parallel::TransferEngine thisNum_engine_zsend1;
  Parallel::TransferEngine thisNum_engine_zrecv1;
}

using namespace Parallel::Timer;

static void
rebuildEngine (IndexType & maxSend,
	       Parallel::TransNumAtomInSubList & sendNum,
	       Parallel::TransSubListData & sendData,
	       Parallel::TransferEngine & engine_send);
static void
rebuildEnginePersistentSend (IndexType & maxSend,
			     Parallel::TransNumAtomInSubList & sendNum,
			     Parallel::TransSubListData & sendData,
			     const int dest,
			     const int tag,
			     Parallel::TransferEngine & engine_send);
static void
rebuildEnginePersistentRecv (IndexType & maxRecv,
			     Parallel::TransNumAtomInSubList & recvNum,
			     Parallel::TransSubListData & recvData,
			     const int dest,
			     const int tag,
			     Parallel::TransferEngine & engine_recv);


Parallel::SystemRedistributeTransferUtil::
SystemRedistributeTransferUtil ()
    : ptr_hdata (NULL), ptr_buff(NULL), mask (MDDataItemMask_All),
      maxNum_xsend0(0), maxNum_xrecv0(0),
      maxNum_xsend1(0), maxNum_xrecv1(0),
      maxNum_ysend0(0), maxNum_yrecv0(0),
      maxNum_ysend1(0), maxNum_yrecv1(0),
      maxNum_zsend0(0), maxNum_zrecv0(0),
      maxNum_zsend1(0), maxNum_zrecv1(0)
{
  Parallel::Interface::shiftNeighbor (CoordXIndex,  1, xsrc0, xdest0);
  Parallel::Interface::shiftNeighbor (CoordXIndex, -1, xsrc1, xdest1);
  Parallel::Interface::shiftNeighbor (CoordYIndex,  1, ysrc0, ydest0);
  Parallel::Interface::shiftNeighbor (CoordYIndex, -1, ysrc1, ydest1);
  Parallel::Interface::shiftNeighbor (CoordZIndex,  1, zsrc0, zdest0);
  Parallel::Interface::shiftNeighbor (CoordZIndex, -1, zsrc1, zdest1);

  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_xsend0.registerBuff ((void*)&thisNum_xsend0, sizeof(IndexType));
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_xrecv0.registerBuff ((void*)&thisNum_xrecv0, sizeof(IndexType));
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_xsend1.registerBuff ((void*)&thisNum_xsend1, sizeof(IndexType));
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_xrecv1.registerBuff ((void*)&thisNum_xrecv1, sizeof(IndexType));
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_ysend0.registerBuff ((void*)&thisNum_ysend0, sizeof(IndexType));
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_yrecv0.registerBuff ((void*)&thisNum_yrecv0, sizeof(IndexType));
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_ysend1.registerBuff ((void*)&thisNum_ysend1, sizeof(IndexType));
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_yrecv1.registerBuff ((void*)&thisNum_yrecv1, sizeof(IndexType));
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_zsend0.registerBuff ((void*)&thisNum_zsend0, sizeof(IndexType));
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_zrecv0.registerBuff ((void*)&thisNum_zrecv0, sizeof(IndexType));
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_zsend1.registerBuff ((void*)&thisNum_zsend1, sizeof(IndexType));
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_zrecv1.registerBuff ((void*)&thisNum_zrecv1, sizeof(IndexType));
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_xsend0.build();
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_xrecv0.build();
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_xsend1.build();
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_xrecv1.build();
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_ysend0.build();
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_yrecv0.build();
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_ysend1.build();
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_yrecv1.build();
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_zsend0.build();
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_zrecv0.build();
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_zsend1.build();
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_zrecv1.build();
}

void Parallel::SystemRedistributeTransferUtil::
clear ()
{
  SystemRedistributeTransferUtil_globalEngine::engine_xsend0.clearRegistered();
  SystemRedistributeTransferUtil_globalEngine::engine_xrecv0.clearRegistered();
  SystemRedistributeTransferUtil_globalEngine::engine_xsend1.clearRegistered();
  SystemRedistributeTransferUtil_globalEngine::engine_xrecv1.clearRegistered();
  SystemRedistributeTransferUtil_globalEngine::engine_ysend0.clearRegistered();
  SystemRedistributeTransferUtil_globalEngine::engine_yrecv0.clearRegistered();
  SystemRedistributeTransferUtil_globalEngine::engine_ysend1.clearRegistered();
  SystemRedistributeTransferUtil_globalEngine::engine_yrecv1.clearRegistered();
  SystemRedistributeTransferUtil_globalEngine::engine_zsend0.clearRegistered();
  SystemRedistributeTransferUtil_globalEngine::engine_zrecv0.clearRegistered();
  SystemRedistributeTransferUtil_globalEngine::engine_zsend1.clearRegistered();
  SystemRedistributeTransferUtil_globalEngine::engine_zrecv1.clearRegistered();

  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_xsend0.clearRegistered();
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_xrecv0.clearRegistered();
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_xsend1.clearRegistered();
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_xrecv1.clearRegistered();
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_ysend0.clearRegistered();
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_yrecv0.clearRegistered();
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_ysend1.clearRegistered();
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_yrecv1.clearRegistered();
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_zsend0.clearRegistered();
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_zrecv0.clearRegistered();
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_zsend1.clearRegistered();
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_zrecv1.clearRegistered();
}

Parallel::SystemRedistributeTransferUtil::
~SystemRedistributeTransferUtil ()
{
  clear ();
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

  rebuildEnginePersistentSend (maxNum_xsend0, xsendNum0, xsendData0,
			       xdest0, 0, SystemRedistributeTransferUtil_globalEngine::engine_xsend0);
  rebuildEnginePersistentRecv (maxNum_xrecv0, xrecvNum0, xrecvData0,
			       xsrc0,  0, SystemRedistributeTransferUtil_globalEngine::engine_xrecv0);
  rebuildEnginePersistentSend (maxNum_xsend1, xsendNum1, xsendData1,
			       xdest1, 1, SystemRedistributeTransferUtil_globalEngine::engine_xsend1);
  rebuildEnginePersistentRecv (maxNum_xrecv1, xrecvNum1, xrecvData1,
			       xsrc1,  1, SystemRedistributeTransferUtil_globalEngine::engine_xrecv1);
  rebuildEnginePersistentSend (maxNum_ysend0, ysendNum0, ysendData0,
			       ydest0, 0, SystemRedistributeTransferUtil_globalEngine::engine_ysend0);
  rebuildEnginePersistentRecv (maxNum_yrecv0, yrecvNum0, yrecvData0,
			       ysrc0,  0, SystemRedistributeTransferUtil_globalEngine::engine_yrecv0);
  rebuildEnginePersistentSend (maxNum_ysend1, ysendNum1, ysendData1,
			       ydest1, 1, SystemRedistributeTransferUtil_globalEngine::engine_ysend1);
  rebuildEnginePersistentRecv (maxNum_yrecv1, yrecvNum1, yrecvData1,
			       ysrc1,  1, SystemRedistributeTransferUtil_globalEngine::engine_yrecv1);
  rebuildEnginePersistentSend (maxNum_zsend0, zsendNum0, zsendData0,
			       zdest0, 0, SystemRedistributeTransferUtil_globalEngine::engine_zsend0);
  rebuildEnginePersistentRecv (maxNum_zrecv0, zrecvNum0, zrecvData0,
			       zsrc0,  0, SystemRedistributeTransferUtil_globalEngine::engine_zrecv0);
  rebuildEnginePersistentSend (maxNum_zsend1, zsendNum1, zsendData1,
			       zdest1, 1, SystemRedistributeTransferUtil_globalEngine::engine_zsend1);
  rebuildEnginePersistentRecv (maxNum_zrecv1, zrecvNum1, zrecvData1,
			       zsrc1,  1, SystemRedistributeTransferUtil_globalEngine::engine_zrecv1);  
}

void Parallel::SystemRedistributeTransferUtil::
calTransNumX ()
{
  thisNum_xsend0 = xsendNum0.getMaxNum();
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_xsend0.Isend (xdest0, 0);
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_xrecv0.Irecv (xsrc0,  0);
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_xsend0.wait();
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_xrecv0.wait();
  thisNum_xsend1 = xsendNum1.getMaxNum();
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_xsend1.Isend (xdest1, 0);
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_xrecv1.Irecv (xsrc1,  0);
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_xsend1.wait();
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_xrecv1.wait();
}

void Parallel::SystemRedistributeTransferUtil::
calTransNumY ()
{
  thisNum_ysend0 = ysendNum0.getMaxNum();
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_ysend0.Isend (ydest0, 0);
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_yrecv0.Irecv (ysrc0,  0);
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_ysend0.wait();
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_yrecv0.wait();
  thisNum_ysend1 = ysendNum1.getMaxNum();
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_ysend1.Isend (ydest1, 0);
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_yrecv1.Irecv (ysrc1,  0);
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_ysend1.wait();
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_yrecv1.wait();
}

void Parallel::SystemRedistributeTransferUtil::
calTransNumZ ()
{
  thisNum_zsend0 = zsendNum0.getMaxNum();
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_zsend0.Isend (zdest0, 0);
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_zrecv0.Irecv (zsrc0,  0);
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_zsend0.wait();
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_zrecv0.wait();
  thisNum_zsend1 = zsendNum1.getMaxNum();
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_zsend1.Isend (zdest1, 0);
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_zrecv1.Irecv (zsrc1,  0);
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_zsend1.wait();
  SystemRedistributeTransferUtil_globalEngine::thisNum_engine_zrecv1.wait();
}

static void
rebuildEngine (IndexType & maxSend,
	       Parallel::TransNumAtomInSubList & sendNum,
	       Parallel::TransSubListData & sendData,
	       Parallel::TransferEngine & engine_send)
{
  engine_send.clearRegistered();
  engine_send.registerBuff (sendNum);
  sendData.build (maxSend);
  engine_send.registerBuff (sendData);
  engine_send.build ();
}

static void
rebuildEnginePersistentSend (IndexType & maxSend,
			     Parallel::TransNumAtomInSubList & sendNum,
			     Parallel::TransSubListData & sendData,
			     const int dest,
			     const int tag,
			     Parallel::TransferEngine & engine_send)
{
  engine_send.clearRegistered();
  engine_send.registerBuff (sendNum);
  sendData.build (maxSend);
  engine_send.registerBuff (sendData);
  engine_send.buildPersistentSend (dest, tag);
}

static void
rebuildEnginePersistentRecv (IndexType & maxRecv,
			     Parallel::TransNumAtomInSubList & recvNum,
			     Parallel::TransSubListData & recvData,
			     const int src,
			     const int tag,
			     Parallel::TransferEngine & engine_recv)
{
  engine_recv.clearRegistered();
  engine_recv.registerBuff (recvNum);
  recvData.build (maxRecv);
  engine_recv.registerBuff (recvData);
  engine_recv.buildPersistentRecv (src, tag);
}

void Parallel::SystemRedistributeTransferUtil::
redistributeHost ()
{
  HostTimer::tic (item_Redistribute_SyncNum);
  calTransNumX ();
  HostTimer::toc (item_Redistribute_SyncNum);
  HostTimer::tic (item_Redistribute_BuildEngine);
  if (thisNum_xsend0 > maxNum_xsend0){
    maxNum_xsend0 = thisNum_xsend0;
    rebuildEnginePersistentSend (maxNum_xsend0, xsendNum0, xsendData0,
				 xdest0, 0, SystemRedistributeTransferUtil_globalEngine::engine_xsend0);
  }
  if (thisNum_xrecv0 > maxNum_xrecv0){
    maxNum_xrecv0 = thisNum_xrecv0;
    rebuildEnginePersistentRecv (maxNum_xrecv0, xrecvNum0, xrecvData0,
				 xsrc0,  0, SystemRedistributeTransferUtil_globalEngine::engine_xrecv0);
  }
  if (thisNum_xsend1 > maxNum_xsend1){
    maxNum_xsend1 = thisNum_xsend1;
    rebuildEnginePersistentSend (maxNum_xsend1, xsendNum1, xsendData1,
				 xdest1, 1, SystemRedistributeTransferUtil_globalEngine::engine_xsend1);
  }
  if (thisNum_xrecv1 > maxNum_xrecv1){
    maxNum_xrecv1 = thisNum_xrecv1;
    rebuildEnginePersistentRecv (maxNum_xrecv1, xrecvNum1, xrecvData1,
				 xsrc1,  1, SystemRedistributeTransferUtil_globalEngine::engine_xrecv1);
  }
  HostTimer::toc (item_Redistribute_BuildEngine);
  HostTimer::tic (item_Redistribute_Transfer);
  if (thisNum_xsend0 != 0) SystemRedistributeTransferUtil_globalEngine::engine_xsend0.startPersistentRequest();
  if (thisNum_xrecv0 != 0) SystemRedistributeTransferUtil_globalEngine::engine_xrecv0.startPersistentRequest();
  if (thisNum_xsend0 != 0) SystemRedistributeTransferUtil_globalEngine::engine_xsend0.wait();
  if (thisNum_xrecv0 != 0) SystemRedistributeTransferUtil_globalEngine::engine_xrecv0.wait();
  if (thisNum_xsend0 != 0) xsend0.clearData();
  if (thisNum_xrecv0 != 0) xrecv0h.add (xrecv0, mask);
  if (thisNum_xsend1 != 0) SystemRedistributeTransferUtil_globalEngine::engine_xsend1.startPersistentRequest();
  if (thisNum_xrecv1 != 0) SystemRedistributeTransferUtil_globalEngine::engine_xrecv1.startPersistentRequest();
  if (thisNum_xsend1 != 0) SystemRedistributeTransferUtil_globalEngine::engine_xsend1.wait();
  if (thisNum_xrecv1 != 0) SystemRedistributeTransferUtil_globalEngine::engine_xrecv1.wait();
  if (thisNum_xsend1 != 0) xsend1.clearData();
  if (thisNum_xrecv1 != 0) xrecv1h.add (xrecv1, mask);
  HostTimer::toc (item_Redistribute_Transfer);

  Parallel::Interface::barrier();

  HostTimer::tic (item_Redistribute_SyncNum);
  calTransNumY ();
  HostTimer::toc (item_Redistribute_SyncNum);
  HostTimer::tic (item_Redistribute_BuildEngine);
  if (thisNum_ysend0 > maxNum_ysend0){
    maxNum_ysend0 = thisNum_ysend0;
    rebuildEnginePersistentSend (maxNum_ysend0, ysendNum0, ysendData0,
				 ydest0, 0, SystemRedistributeTransferUtil_globalEngine::engine_ysend0);
  }
  if (thisNum_yrecv0 > maxNum_yrecv0){
    maxNum_yrecv0 = thisNum_yrecv0;
    rebuildEnginePersistentRecv (maxNum_yrecv0, yrecvNum0, yrecvData0,
				 ysrc0,  0, SystemRedistributeTransferUtil_globalEngine::engine_yrecv0);
  }
  if (thisNum_ysend1 > maxNum_ysend1){
    maxNum_ysend1 = thisNum_ysend1;
    rebuildEnginePersistentSend (maxNum_ysend1, ysendNum1, ysendData1,
				 ydest1, 1, SystemRedistributeTransferUtil_globalEngine::engine_ysend1);
  }
  if (thisNum_yrecv1 > maxNum_yrecv1){
    maxNum_yrecv1 = thisNum_yrecv1;
    rebuildEnginePersistentRecv (maxNum_yrecv1, yrecvNum1, yrecvData1,
				 ysrc1,  1, SystemRedistributeTransferUtil_globalEngine::engine_yrecv1);
  }
  HostTimer::toc (item_Redistribute_BuildEngine);
  HostTimer::tic (item_Redistribute_Transfer);
  if (thisNum_ysend0 != 0) SystemRedistributeTransferUtil_globalEngine::engine_ysend0.startPersistentRequest();
  if (thisNum_yrecv0 != 0) SystemRedistributeTransferUtil_globalEngine::engine_yrecv0.startPersistentRequest();
  if (thisNum_ysend0 != 0) SystemRedistributeTransferUtil_globalEngine::engine_ysend0.wait();
  if (thisNum_yrecv0 != 0) SystemRedistributeTransferUtil_globalEngine::engine_yrecv0.wait();
  if (thisNum_ysend0 != 0) ysend0.clearData();
  if (thisNum_yrecv0 != 0) yrecv0h.add (yrecv0, mask);
  if (thisNum_ysend1 != 0) SystemRedistributeTransferUtil_globalEngine::engine_ysend1.startPersistentRequest();
  if (thisNum_yrecv1 != 0) SystemRedistributeTransferUtil_globalEngine::engine_yrecv1.startPersistentRequest();
  if (thisNum_ysend1 != 0) SystemRedistributeTransferUtil_globalEngine::engine_ysend1.wait();
  if (thisNum_yrecv1 != 0) SystemRedistributeTransferUtil_globalEngine::engine_yrecv1.wait();
  if (thisNum_ysend1 != 0) ysend1.clearData();
  if (thisNum_yrecv1 != 0) yrecv1h.add (yrecv1, mask);
  HostTimer::toc (item_Redistribute_Transfer);

  Parallel::Interface::barrier();

  HostTimer::tic (item_Redistribute_SyncNum);
  calTransNumZ ();
  HostTimer::toc (item_Redistribute_SyncNum);
  HostTimer::tic (item_Redistribute_BuildEngine);  
  if (thisNum_zsend0 > maxNum_zsend0){
    maxNum_zsend0 = thisNum_zsend0;
    rebuildEnginePersistentSend (maxNum_zsend0, zsendNum0, zsendData0,
				 zdest0, 0, SystemRedistributeTransferUtil_globalEngine::engine_zsend0);
  }
  if (thisNum_zrecv0 > maxNum_zrecv0){
    maxNum_zrecv0 = thisNum_zrecv0;
    rebuildEnginePersistentRecv (maxNum_zrecv0, zrecvNum0, zrecvData0,
				 zsrc0,  0, SystemRedistributeTransferUtil_globalEngine::engine_zrecv0);
  }
  if (thisNum_zsend1 > maxNum_zsend1){
    maxNum_zsend1 = thisNum_zsend1;
    rebuildEnginePersistentSend (maxNum_zsend1, zsendNum1, zsendData1,
				 zdest1, 1, SystemRedistributeTransferUtil_globalEngine::engine_zsend1);
  }
  if (thisNum_zrecv1 > maxNum_zrecv1){
    maxNum_zrecv1 = thisNum_zrecv1;
    rebuildEnginePersistentRecv (maxNum_zrecv1, zrecvNum1, zrecvData1,
				 zsrc1,  1, SystemRedistributeTransferUtil_globalEngine::engine_zrecv1);
  }
  HostTimer::toc (item_Redistribute_BuildEngine);
  HostTimer::tic (item_Redistribute_Transfer);
  if (thisNum_zsend0 != 0) SystemRedistributeTransferUtil_globalEngine::engine_zsend0.startPersistentRequest();
  if (thisNum_zrecv0 != 0) SystemRedistributeTransferUtil_globalEngine::engine_zrecv0.startPersistentRequest();
  if (thisNum_zsend0 != 0) SystemRedistributeTransferUtil_globalEngine::engine_zsend0.wait();
  if (thisNum_zrecv0 != 0) SystemRedistributeTransferUtil_globalEngine::engine_zrecv0.wait();
  if (thisNum_zsend0 != 0) zsend0.clearData();
  if (thisNum_zrecv0 != 0) zrecv0h.add (zrecv0, mask);
  if (thisNum_zsend1 != 0) SystemRedistributeTransferUtil_globalEngine::engine_zsend1.startPersistentRequest();
  if (thisNum_zrecv1 != 0) SystemRedistributeTransferUtil_globalEngine::engine_zrecv1.startPersistentRequest();
  if (thisNum_zsend1 != 0) SystemRedistributeTransferUtil_globalEngine::engine_zsend1.wait();
  if (thisNum_zrecv1 != 0) SystemRedistributeTransferUtil_globalEngine::engine_zrecv1.wait();
  if (thisNum_zsend1 != 0) zsend1.clearData();
  if (thisNum_zrecv1 != 0) zrecv1h.add (zrecv1, mask);
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
	    MDDataItemMask_GlobalIndex),
      maxNum_xsend0(0), maxNum_xrecv0(0),
      maxNum_xsend1(0), maxNum_xrecv1(0),
      maxNum_ysend0(0), maxNum_yrecv0(0),
      maxNum_ysend1(0), maxNum_yrecv1(0),
      maxNum_zsend0(0), maxNum_zrecv0(0),
      maxNum_zsend1(0), maxNum_zrecv1(0)
{
  Parallel::Interface::shiftNeighbor (CoordXIndex,  1, xsrc0, xdest0);
  Parallel::Interface::shiftNeighbor (CoordXIndex, -1, xsrc1, xdest1);
  Parallel::Interface::shiftNeighbor (CoordYIndex,  1, ysrc0, ydest0);
  Parallel::Interface::shiftNeighbor (CoordYIndex, -1, ysrc1, ydest1);
  Parallel::Interface::shiftNeighbor (CoordZIndex,  1, zsrc0, zdest0);
  Parallel::Interface::shiftNeighbor (CoordZIndex, -1, zsrc1, zdest1);

  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_xsend0.registerBuff ((void*)&thisNum_xsend0, sizeof(IndexType));
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_xrecv0.registerBuff ((void*)&thisNum_xrecv0, sizeof(IndexType));
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_xsend1.registerBuff ((void*)&thisNum_xsend1, sizeof(IndexType));
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_xrecv1.registerBuff ((void*)&thisNum_xrecv1, sizeof(IndexType));
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_ysend0.registerBuff ((void*)&thisNum_ysend0, sizeof(IndexType));
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_yrecv0.registerBuff ((void*)&thisNum_yrecv0, sizeof(IndexType));
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_ysend1.registerBuff ((void*)&thisNum_ysend1, sizeof(IndexType));
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_yrecv1.registerBuff ((void*)&thisNum_yrecv1, sizeof(IndexType));
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_zsend0.registerBuff ((void*)&thisNum_zsend0, sizeof(IndexType));
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_zrecv0.registerBuff ((void*)&thisNum_zrecv0, sizeof(IndexType));
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_zsend1.registerBuff ((void*)&thisNum_zsend1, sizeof(IndexType));
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_zrecv1.registerBuff ((void*)&thisNum_zrecv1, sizeof(IndexType));
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_xsend0.build();
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_xrecv0.build();
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_xsend1.build();
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_xrecv1.build();
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_ysend0.build();
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_yrecv0.build();
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_ysend1.build();
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_yrecv1.build();
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_zsend0.build();
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_zrecv0.build();
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_zsend1.build();
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_zrecv1.build();

}

Parallel::SystemTransCoordsTransferUtil::
~SystemTransCoordsTransferUtil()
{
  clear();
}

void Parallel::SystemTransCoordsTransferUtil::
clear ()
{
  SystemTransCoordsTransferUtil_globalEngine::engine_xsend0.clearRegistered();
  SystemTransCoordsTransferUtil_globalEngine::engine_xrecv0.clearRegistered();
  SystemTransCoordsTransferUtil_globalEngine::engine_xsend1.clearRegistered();
  SystemTransCoordsTransferUtil_globalEngine::engine_xrecv1.clearRegistered();
  SystemTransCoordsTransferUtil_globalEngine::engine_ysend0.clearRegistered();
  SystemTransCoordsTransferUtil_globalEngine::engine_yrecv0.clearRegistered();
  SystemTransCoordsTransferUtil_globalEngine::engine_ysend1.clearRegistered();
  SystemTransCoordsTransferUtil_globalEngine::engine_yrecv1.clearRegistered();
  SystemTransCoordsTransferUtil_globalEngine::engine_zsend0.clearRegistered();
  SystemTransCoordsTransferUtil_globalEngine::engine_zrecv0.clearRegistered();
  SystemTransCoordsTransferUtil_globalEngine::engine_zsend1.clearRegistered();
  SystemTransCoordsTransferUtil_globalEngine::engine_zrecv1.clearRegistered();

  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_xsend0.clearRegistered();
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_xrecv0.clearRegistered();
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_xsend1.clearRegistered();
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_xrecv1.clearRegistered();
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_ysend0.clearRegistered();
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_yrecv0.clearRegistered();
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_ysend1.clearRegistered();
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_yrecv1.clearRegistered();
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_zsend0.clearRegistered();
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_zrecv0.clearRegistered();
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_zsend1.clearRegistered();
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_zrecv1.clearRegistered();
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

  rebuildEnginePersistentSend (maxNum_xsend0, xsendNum0, xsendData0,
			       xdest0, 0, SystemTransCoordsTransferUtil_globalEngine::engine_xsend0);
  rebuildEnginePersistentRecv (maxNum_xrecv0, xrecvNum0, xrecvData0,
			       xsrc0,  0, SystemTransCoordsTransferUtil_globalEngine::engine_xrecv0);
  rebuildEnginePersistentSend (maxNum_xsend1, xsendNum1, xsendData1,
			       xdest1, 1, SystemTransCoordsTransferUtil_globalEngine::engine_xsend1);
  rebuildEnginePersistentRecv (maxNum_xrecv1, xrecvNum1, xrecvData1,
			       xsrc1,  1, SystemTransCoordsTransferUtil_globalEngine::engine_xrecv1);
  rebuildEnginePersistentSend (maxNum_ysend0, ysendNum0, ysendData0,
			       ydest0, 0, SystemTransCoordsTransferUtil_globalEngine::engine_ysend0);
  rebuildEnginePersistentRecv (maxNum_yrecv0, yrecvNum0, yrecvData0,
			       ysrc0,  0, SystemTransCoordsTransferUtil_globalEngine::engine_yrecv0);
  rebuildEnginePersistentSend (maxNum_ysend1, ysendNum1, ysendData1,
			       ydest1, 1, SystemTransCoordsTransferUtil_globalEngine::engine_ysend1);
  rebuildEnginePersistentRecv (maxNum_yrecv1, yrecvNum1, yrecvData1,
			       ysrc1,  1, SystemTransCoordsTransferUtil_globalEngine::engine_yrecv1);
  rebuildEnginePersistentSend (maxNum_zsend0, zsendNum0, zsendData0,
			       zdest0, 0, SystemTransCoordsTransferUtil_globalEngine::engine_zsend0);
  rebuildEnginePersistentRecv (maxNum_zrecv0, zrecvNum0, zrecvData0,
			       zsrc0,  0, SystemTransCoordsTransferUtil_globalEngine::engine_zrecv0);
  rebuildEnginePersistentSend (maxNum_zsend1, zsendNum1, zsendData1,
			       zdest1, 1, SystemTransCoordsTransferUtil_globalEngine::engine_zsend1);
  rebuildEnginePersistentRecv (maxNum_zrecv1, zrecvNum1, zrecvData1,
			       zsrc1,  1, SystemTransCoordsTransferUtil_globalEngine::engine_zrecv1);  
  
}

void Parallel::SystemTransCoordsTransferUtil::
calTransNumX ()
{
  thisNum_xsend0 = xsendNum0.getMaxNum();
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_xsend0.Isend (xdest0, 0);
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_xrecv0.Irecv (xsrc0,  0);
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_xsend0.wait();
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_xrecv0.wait();
  thisNum_xsend1 = xsendNum1.getMaxNum();
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_xsend1.Isend (xdest1, 0);
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_xrecv1.Irecv (xsrc1,  0);
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_xsend1.wait();
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_xrecv1.wait();
}

void Parallel::SystemTransCoordsTransferUtil::
calTransNumY ()
{
  thisNum_ysend0 = ysendNum0.getMaxNum();
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_ysend0.Isend (ydest0, 0);
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_yrecv0.Irecv (ysrc0,  0);
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_ysend0.wait();
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_yrecv0.wait();
  thisNum_ysend1 = ysendNum1.getMaxNum();
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_ysend1.Isend (ydest1, 0);
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_yrecv1.Irecv (ysrc1,  0);
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_ysend1.wait();
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_yrecv1.wait();
}

void Parallel::SystemTransCoordsTransferUtil::
calTransNumZ ()
{
  thisNum_zsend0 = zsendNum0.getMaxNum();
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_zsend0.Isend (zdest0, 0);
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_zrecv0.Irecv (zsrc0,  0);
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_zsend0.wait();
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_zrecv0.wait();
  thisNum_zsend1 = zsendNum1.getMaxNum();
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_zsend1.Isend (zdest1, 0);
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_zrecv1.Irecv (zsrc1,  0);
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_zsend1.wait();
  SystemTransCoordsTransferUtil_globalEngine::thisNum_engine_zrecv1.wait();
}


void Parallel::SystemTransCoordsTransferUtil::
transCoords ()
{
  calTransNumX ();
  if (thisNum_xsend0 > maxNum_xsend0){
    maxNum_xsend0 = thisNum_xsend0;
    rebuildEnginePersistentSend (maxNum_xsend0, xsendNum0, xsendData0,
				 xdest0, 0, SystemTransCoordsTransferUtil_globalEngine::engine_xsend0);
  }
  if (thisNum_xrecv0 > maxNum_xrecv0){
    maxNum_xrecv0 = thisNum_xrecv0;
    rebuildEnginePersistentRecv (maxNum_xrecv0, xrecvNum0, xrecvData0,
				 xsrc0,  0, SystemTransCoordsTransferUtil_globalEngine::engine_xrecv0);
  }
  if (thisNum_xsend1 > maxNum_xsend1){
    maxNum_xsend1 = thisNum_xsend1;
    rebuildEnginePersistentSend (maxNum_xsend1, xsendNum1, xsendData1,
				 xdest1, 1, SystemTransCoordsTransferUtil_globalEngine::engine_xsend1);
  }
  if (thisNum_xrecv1 > maxNum_xrecv1){
    maxNum_xrecv1 = thisNum_xrecv1;
    rebuildEnginePersistentRecv (maxNum_xrecv1, xrecvNum1, xrecvData1,
				 xsrc1,  1, SystemTransCoordsTransferUtil_globalEngine::engine_xrecv1);
  }
  SystemTransCoordsTransferUtil_globalEngine::engine_xsend0.startPersistentRequest();
  SystemTransCoordsTransferUtil_globalEngine::engine_xrecv0.startPersistentRequest();
  SystemTransCoordsTransferUtil_globalEngine::engine_xsend0.wait();
  SystemTransCoordsTransferUtil_globalEngine::engine_xrecv0.wait();
  SystemTransCoordsTransferUtil_globalEngine::engine_xsend1.startPersistentRequest();
  SystemTransCoordsTransferUtil_globalEngine::engine_xrecv1.startPersistentRequest();
  SystemTransCoordsTransferUtil_globalEngine::engine_xsend1.wait();
  SystemTransCoordsTransferUtil_globalEngine::engine_xrecv1.wait();

  Parallel::Interface::barrier();
  
  calTransNumY ();
  if (thisNum_ysend0 > maxNum_ysend0){
    maxNum_ysend0 = thisNum_ysend0;
    rebuildEnginePersistentSend (maxNum_ysend0, ysendNum0, ysendData0,
				 ydest0, 0, SystemTransCoordsTransferUtil_globalEngine::engine_ysend0);
  }
  if (thisNum_yrecv0 > maxNum_yrecv0){
    maxNum_yrecv0 = thisNum_yrecv0;
    rebuildEnginePersistentRecv (maxNum_yrecv0, yrecvNum0, yrecvData0,
				 ysrc0,  0, SystemTransCoordsTransferUtil_globalEngine::engine_yrecv0);
  }
  if (thisNum_ysend1 > maxNum_ysend1){
    maxNum_ysend1 = thisNum_ysend1;
    rebuildEnginePersistentSend (maxNum_ysend1, ysendNum1, ysendData1,
				 ydest1, 1, SystemTransCoordsTransferUtil_globalEngine::engine_ysend1);
  }
  if (thisNum_yrecv1 > maxNum_yrecv1){
    maxNum_yrecv1 = thisNum_yrecv1;
    rebuildEnginePersistentRecv (maxNum_yrecv1, yrecvNum1, yrecvData1,
				 ysrc1,  1, SystemTransCoordsTransferUtil_globalEngine::engine_yrecv1);
  }
  SystemTransCoordsTransferUtil_globalEngine::engine_ysend0.startPersistentRequest();
  SystemTransCoordsTransferUtil_globalEngine::engine_yrecv0.startPersistentRequest();
  SystemTransCoordsTransferUtil_globalEngine::engine_ysend0.wait();
  SystemTransCoordsTransferUtil_globalEngine::engine_yrecv0.wait();
  SystemTransCoordsTransferUtil_globalEngine::engine_ysend1.startPersistentRequest();
  SystemTransCoordsTransferUtil_globalEngine::engine_yrecv1.startPersistentRequest();
  SystemTransCoordsTransferUtil_globalEngine::engine_ysend1.wait();
  SystemTransCoordsTransferUtil_globalEngine::engine_yrecv1.wait();

  Parallel::Interface::barrier();

  calTransNumZ ();
  if (thisNum_zsend0 > maxNum_zsend0){
    maxNum_zsend0 = thisNum_zsend0;
    rebuildEnginePersistentSend (maxNum_zsend0, zsendNum0, zsendData0,
				 zdest0, 0, SystemTransCoordsTransferUtil_globalEngine::engine_zsend0);
  }
  if (thisNum_zrecv0 > maxNum_zrecv0){
    maxNum_zrecv0 = thisNum_zrecv0;
    rebuildEnginePersistentRecv (maxNum_zrecv0, zrecvNum0, zrecvData0,
				 zsrc0,  0, SystemTransCoordsTransferUtil_globalEngine::engine_zrecv0);
  }
  if (thisNum_zsend1 > maxNum_zsend1){
    maxNum_zsend1 = thisNum_zsend1;
    rebuildEnginePersistentSend (maxNum_zsend1, zsendNum1, zsendData1,
				 zdest1, 1, SystemTransCoordsTransferUtil_globalEngine::engine_zsend1);
  }
  if (thisNum_zrecv1 > maxNum_zrecv1){
    maxNum_zrecv1 = thisNum_zrecv1;
    rebuildEnginePersistentRecv (maxNum_zrecv1, zrecvNum1, zrecvData1,
				 zsrc1,  1, SystemTransCoordsTransferUtil_globalEngine::engine_zrecv1);
  }
  SystemTransCoordsTransferUtil_globalEngine::engine_zsend0.startPersistentRequest();
  SystemTransCoordsTransferUtil_globalEngine::engine_zrecv0.startPersistentRequest();
  SystemTransCoordsTransferUtil_globalEngine::engine_zsend0.wait();
  SystemTransCoordsTransferUtil_globalEngine::engine_zrecv0.wait();
  SystemTransCoordsTransferUtil_globalEngine::engine_zsend1.startPersistentRequest();
  SystemTransCoordsTransferUtil_globalEngine::engine_zrecv1.startPersistentRequest();
  SystemTransCoordsTransferUtil_globalEngine::engine_zsend1.wait();
  SystemTransCoordsTransferUtil_globalEngine::engine_zrecv1.wait();


  // sender.clearRegistered();
  // recver.clearRegistered();
  // sender.registerBuff (xsendNum0);
  // recver.registerBuff (xrecvNum0);
  // sender.build ();
  // recver.build ();
  // sender.Isend (xdest0, 0);
  // recver.Irecv (xsrc0,  0);
  // sender.wait();
  // recver.wait();

  // xsendData0.build ();
  // xrecvData0.build ();
  // sender.clearRegistered();
  // recver.clearRegistered();
  // sender.registerBuff (xsendData0);
  // recver.registerBuff (xrecvData0);
  // sender.build ();
  // recver.build ();
  // sender.Isend (xdest0, 1);
  // recver.Irecv (xsrc0,  1);
  // sender.wait();
  // recver.wait();

  // sender.clearRegistered();
  // recver.clearRegistered();
  // sender.registerBuff (xsendNum1);
  // recver.registerBuff (xrecvNum1);
  // sender.build ();
  // recver.build ();
  // sender.Isend (xdest1, 2);
  // recver.Irecv (xsrc1,  2);
  // sender.wait();
  // recver.wait();  

  // xsendData1.build ();
  // xrecvData1.build ();
  // sender.clearRegistered();
  // recver.clearRegistered();
  // sender.registerBuff (xsendData1);
  // recver.registerBuff (xrecvData1);
  // sender.build ();
  // recver.build ();
  // sender.Isend (xdest1, 3);
  // recver.Irecv (xsrc1,  3);
  // sender.wait();
  // recver.wait();

  // Parallel::Interface::barrier();

  // sender.clearRegistered();
  // recver.clearRegistered();
  // sender.registerBuff (ysendNum0);
  // recver.registerBuff (yrecvNum0);
  // sender.build ();
  // recver.build ();
  // sender.Isend (ydest0, 0);
  // recver.Irecv (ysrc0,  0);
  // sender.wait();
  // recver.wait();  

  // ysendData0.build ();
  // yrecvData0.build ();
  // sender.clearRegistered();
  // recver.clearRegistered();
  // sender.registerBuff (ysendData0);
  // recver.registerBuff (yrecvData0);
  // sender.build ();
  // recver.build ();
  // sender.Isend (ydest0, 1);
  // recver.Irecv (ysrc0,  1);
  // sender.wait();
  // recver.wait();

  // sender.clearRegistered();
  // recver.clearRegistered();
  // sender.registerBuff (ysendNum1);
  // recver.registerBuff (yrecvNum1);
  // sender.build ();
  // recver.build ();
  // sender.Isend (ydest1, 2);
  // recver.Irecv (ysrc1,  2);
  // sender.wait();
  // recver.wait();  

  // ysendData1.build ();
  // yrecvData1.build ();
  // sender.clearRegistered();
  // recver.clearRegistered();
  // sender.registerBuff (ysendData1);
  // recver.registerBuff (yrecvData1);
  // sender.build ();
  // recver.build ();
  // sender.Isend (ydest1, 3);
  // recver.Irecv (ysrc1,  3);
  // sender.wait();
  // recver.wait();

  // Parallel::Interface::barrier();
  
  // sender.clearRegistered();
  // recver.clearRegistered();
  // sender.registerBuff (zsendNum0);
  // recver.registerBuff (zrecvNum0);
  // sender.build ();
  // recver.build ();
  // sender.Isend (zdest0, 0);
  // recver.Irecv (zsrc0,  0);
  // sender.wait();
  // recver.wait();  

  // zsendData0.build ();
  // zrecvData0.build ();
  // sender.clearRegistered();
  // recver.clearRegistered();
  // sender.registerBuff (zsendData0);
  // recver.registerBuff (zrecvData0);
  // sender.build ();
  // recver.build ();
  // sender.Isend (zdest0, 1);
  // recver.Irecv (zsrc0,  1);
  // sender.wait();
  // recver.wait();

  // sender.clearRegistered();
  // recver.clearRegistered();
  // sender.registerBuff (zsendNum1);
  // recver.registerBuff (zrecvNum1);
  // sender.build ();
  // recver.build ();
  // sender.Isend (zdest1, 2);
  // recver.Irecv (zsrc1,  2);
  // sender.wait();
  // recver.wait();  

  // zsendData1.build ();
  // zrecvData1.build ();
  // sender.clearRegistered();
  // recver.clearRegistered();
  // sender.registerBuff (zsendData1);
  // recver.registerBuff (zrecvData1);
  // sender.build ();
  // recver.build ();
  // sender.Isend (zdest1, 3);
  // recver.Irecv (zsrc1,  3);
  // sender.wait();
  // recver.wait();

  // Parallel::Interface::barrier();
}
