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
  Displacement_max (const MDSystem & sys,
		    const IndexType & NThread);
  ~Displacement_max ();
  void reinit (const MDSystem & sys,
	       const IndexType & NThread);
  void recordCoord (const MDSystem & sys,
		    MDTimer * timer = NULL);
  ScalorType calMaxDisplacemant (const MDSystem & sys,
				 MDTimer * timer = NULL);
  ScalorType getResult () const {return hresult;}
  virtual void reshuffle (const IndexType * indexTable,
			  const IndexType & numAtom,
			  MDTimer * timer = NULL);
}
    ;


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
  Displacement_mean (const MDSystem & sys,
		    const IndexType & NThread);
  ~Displacement_mean ();
  void reinit (const MDSystem & sys,
	       const IndexType & NThread);
  void recordCoord (const MDSystem & sys,
		    MDTimer * timer = NULL);
  ScalorType calMeanDisplacemant (const MDSystem & sys,
				 MDTimer * timer = NULL);
  ScalorType getResult () const {return hresult;}
  virtual void reshuffle (const IndexType * indexTable,
			  const IndexType & numAtom,
			  MDTimer * timer = NULL);
}
    ;



#endif

