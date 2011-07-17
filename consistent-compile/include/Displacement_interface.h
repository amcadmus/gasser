#ifndef __Displacement_interface_h_wanghan__
#define __Displacement_interface_h_wanghan__

#include "common.h"
#include "BoxGeometry.h"
#include "Reshufflable.h"
#include "Reshuffle_interface.h"
#include "MDSystem.h"
#include "MDSystem_interface.h"
#include "MDError_interface.h"
#include "MDTimer_interface.h"
#include "SumVector.h"
#include "MaxVector.h"

/// Calculate the maximum displacement of atom in the system.

class Displacement_max : public Reshufflable
{
  dim3			atomGridDim;
  dim3			myBlockDim;
  IndexType		displacement_sbuffSize;
  CoordType *		backupCoord;
  CoordType *		bkbackupCoord;
  ScalorType *		dresult;
  ScalorType		hresult;
  MaxVector<ScalorType>	max;
  bool			malloced;
  void mallocDisplacemant (const MDSystem & sys);
  void clearDisplacement ();
public:
  /** 
   * Constructor.
   * 
   * @param sys The MDSystem, which should be the same as those used
   * in the following calls.   
   * @param NThread Number of threads in a block.
   */
  Displacement_max (const MDSystem & sys,
		    const IndexType & NThread);
  ~Displacement_max ();
  /** 
   * Rinitializer.
   * 
   * @param sys The MDSystem, which should be the same as those used
   * in the following calls.   
   * @param NThread Number of threads in a block.
   */
  void reinit (const MDSystem & sys,
	       const IndexType & NThread);
  /** 
   * Record the coordinates of atoms in the system.
   * 
   * @param sys The MDSystem.
   * @param timer Timer measuring the performance.
   */
  void recordCoord (const MDSystem & sys,
		    MDTimer * timer = NULL);
  /** 
   * Calculate the maximum displacement according to the recorded
   * coordinates.
   * 
   * @param sys The MDSystem.
   * @param timer Timer measuring the performance.
   */
  ScalorType calMaxDisplacemant (const MDSystem & sys,
				 MDTimer * timer = NULL);
  /** 
   * Get the maximum displacement.
   * 
   * @return The maximum displacement.
   */
  ScalorType getResult () const {return hresult;}
  virtual void reshuffle (const IndexType * indexTable,
			  const IndexType & numAtom,
			  MDTimer * timer = NULL);
}
    ;

/// Calculate the mean displacement of atom in the system.

class Displacement_mean : public Reshufflable
{
  dim3			atomGridDim;
  dim3			myBlockDim;
  IndexType		displacement_sbuffSize;
  CoordType *		backupCoord;
  CoordType *		bkbackupCoord;
  ScalorType *		dresult;
  ScalorType		hresult;
  SumVector<ScalorType>	sum;
  bool			malloced;
  void mallocDisplacemant (const MDSystem & sys);
  void clearDisplacement ();
public:
  /** 
   * Constructor.
   * 
   * @param sys The MDSystem, which should be the same as those used
   * in the following calls.   
   * @param NThread Number of threads in a block.
   */
  Displacement_mean (const MDSystem & sys,
		    const IndexType & NThread);
  ~Displacement_mean ();
  /** 
   * Reinitializer
   * 
   * @param sys The MDSystem, which should be the same as those used
   * in the following calls.   
   * @param NThread Number of threads in a block.
   */
  void reinit (const MDSystem & sys,
	       const IndexType & NThread);
  /** 
   * Record the coordinates of atoms in the system.
   * 
   * @param sys The MDSystem.
   * @param timer Timer measuring the performance.
   */
  void recordCoord (const MDSystem & sys,
		    MDTimer * timer = NULL);
  /** 
   * Calculate the mean displacement according to the recorded
   * coordinates.
   * 
   * @param sys The MDSystem.
   * @param timer Timer measuring the performance.
   */
  ScalorType calMeanDisplacemant (const MDSystem & sys,
				 MDTimer * timer = NULL);
  /** 
   * Get the mean displacement.
   * 
   * @return The mean displacement.
   */
  ScalorType getResult () const {return hresult;}
  virtual void reshuffle (const IndexType * indexTable,
			  const IndexType & numAtom,
			  MDTimer * timer = NULL);
}
    ;



#endif

