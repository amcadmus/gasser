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

/// The cell list for the short-range non-bonded interaction.

/**
 * The cell list can be built on specified directions. The sub-cell is
 * also supported.
 * 
 */


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
  /** 
   * Constructor
   * 
   * @param sys MDSystem for which the cell list is built.
   * @param cellSize Minimum Cell size. The actural size can be larger
   * than this one. If the 4*cellSize is larger than box size, then
   * no cell list is built.
   * @param NTreadCell Number of threads in a block for cell operations,
   * which should be larger than the maximum possible number of atom
   * in a cell.
   * @param NTreadAtom Number of threads in a block for atom operations
   * @param bdir On which direction the cell list is build. the default
   * is mdRectBoxDirectionX | mdRectBoxDirectionY | mdRectBoxDirectionZ.
   * @param divide Number of subcell on each direction.
   * 
   * @return 
   */
  CellList (const MDSystem &		sys,
	    const ScalorType &		cellSize,
	    const IndexType &		NTreadCell,
	    const IndexType &		NTreadAtom,
	    const BoxDirection_t &	bdir = 7,
	    const IndexType &		divide = 1);
  /** 
   * Destructor.
   */
  ~CellList();
  /** 
   * Reinitializer.
   * 
   * @param sys MDSystem for which the cell list is built.
   * @param cellSize Minimum Cell size. The actural size can be larger
   * than this one. If the 4*cellSize is larger than box size, then
   * no cell list is built.
   * @param NTreadCell Number of threads in a block for cell operations,
   * which should be larger than the maximum possible number of atom
   * in a cell.
   * @param NTreadAtom Number of threads in a block for atom operations
   * @param bdir On which direction the cell list is build. the default
   * is mdRectBoxDirectionX | mdRectBoxDirectionY | mdRectBoxDirectionZ.
   * @param divide Number of sub-cell on each direction.
   * 
   */
  void reinit (const MDSystem &		sys,
	       const ScalorType &	cellSize,
	       const IndexType &	NTreadCell,
	       const IndexType &	NTreadAtom,
	       const BoxDirection_t &	bdir = 7,
	       const IndexType &	divide = 1);
  /** 
   * Rebuild cell list
   * 
   * @param sys MDSystem for which the cell list is built.
   * @param timer Timer measuring the performance.
   */
  void rebuild (const MDSystem & sys,
		MDTimer * timer = NULL);
  /** 
   * Check if the cell list is empty.
   * 
   * 
   * @return true if the cell list is empty.
   */
  bool isempty () const {return !mallocedDeviceCellList;}
  /** 
   * Reshuffle the cell list.
   * 
   * @param indexTable Index table to reshuffle.
   * @param numAtom Number of atom in system, which is the same as the
   * length of the index table
   * @param timer Timer measuring the performance.
   */
  virtual void reshuffle (const IndexType * indexTable,
			  const IndexType & numAtom,
			  MDTimer *timer = NULL);
public:
  /** 
   * Get gridDim of cell operations.
   * 
   * 
   * @return gridDim of cell operations.
   */
  const dim3 & getCellGrimDim () const {return cellGridDim;}
  /** 
   * Get gridDim of atom operations.
   * 
   * 
   * @return gridDim of atom operations.
   */
  const dim3 & getAtomGrimDim () const {return atomGridDim;}
  /** 
   * Get blockDim of cell operations.
   * 
   * 
   * @return blockDim of cell operations.
   */
  const dim3 & getCellBlockDim () const {return cellBlockDim;}
  /** 
   * Get blockDim of atom operations.
   * 
   * 
   * @return blockDim of atom operations.
   */
  const dim3 & getAtomBlockDim () const {return atomBlockDim;}
}
    ;

#endif
