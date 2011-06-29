#ifndef __CellList_interface_h_wanghan__
#define __CellList_interface_h_wanghan__

#include "common.h"
#include "BoxGeometry.h"

#include "MDSystem.h"
#include "MDSystem_interface.h"

#include "common.h"
#include "BoxGeometry.h"
#include "Reshufflable.h"
#include "MDError_interface.h"
#include "MDTimer_interface.h"

using namespace RectangularBoxGeometry;

class CellList : public Reshufflable
{
  dim3			cellGridDim;
  dim3			atomGridDim;
  dim3			cellBlockDim;
  dim3			atomBlockDim;
  IndexType		bitDeepth;
  IndexType *		mySendBuff;
  IndexType *		myTargetBuff;
  bool			mallocedDeviceCellList;
  ScalorType		mycellSize;
  BoxDirection_t	mybdir;
  IndexType		mydivide;
  IndexType		buildDeviceCellList_step1_sbuffSize;
  MDError		err;
public:
  DeviceCellList	dclist;
private:
  void clearDeviceCellList();
  void mallocDeviceCellList (const IntVectorType & NCell,
			     const HostVectorType & boxSize);
  bool DecideNumCell (const MDSystem &		sys,
		      const ScalorType &	cellSize,
		      const BoxDirection_t &	bdir,
		      const IndexType &		divide,
		      IntVectorType &		NCell);
  void naivelyBuildDeviceCellList (const MDSystem & sys);
  void buildDeviceCellList (const MDSystem & sys);
public:
  CellList (const MDSystem &		sys,
	    const ScalorType &		cellSize,
	    const IndexType &		NTreadCell,
	    const IndexType &		NTreadAtom,
	    const BoxDirection_t &	bdir = 7,
	    const IndexType &		divide = 1);
  ~CellList();
  void reinit (const MDSystem &		sys,
	       const ScalorType &	cellSize,
	       const IndexType &	NTreadCell,
	       const IndexType &	NTreadAtom,
	       const BoxDirection_t &	bdir = 7,
	       const IndexType &	divide = 1);
  
  void rebuild (const MDSystem & sys,
		MDTimer * timer = NULL);
  bool isempty () const {return !mallocedDeviceCellList;}
  virtual void reshuffle (const IndexType * indexTable,
			  const IndexType & numAtom,
			  MDTimer *timer = NULL);
public:
  const dim3 & getCellGrimDim () const {return cellGridDim;}
  const dim3 & getAtomGrimDim () const {return atomGridDim;}
  const dim3 & getCellBlockDim () const {return cellBlockDim;}
  const dim3 & getAtomBlockDim () const {return atomBlockDim;}
}
    ;

#endif
