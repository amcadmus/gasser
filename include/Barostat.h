#ifndef __Barostat_h_wanghan__
#define __Barostat_h_wanghan__

#include"common.h"
#include"BoxGeometry.h"

using namespace RectangularBoxGeometry;

class Barostat
{
};

class Barostat_XRescale : public Barostat
{
public:
  virtual void calScale (const ScalorType * nowP,
			 ScalorType * scale) const = 0;
};

class Barostat_Berendsen : public Barostat_XRescale
{
  ScalorType refP[3];
  ScalorType dt;
  ScalorType tau;
  ScalorType beta[3];

  ScalorType scalor;

  IndexType numPcoupleGroup;
  BoxDirection_t PcoupleDirections[3];
  IndexType numDirInGroup [3];
  IndexType dirIndex[3][3];
public:
  void reinit (const ScalorType & dt,
	       const ScalorType & tau_P);
  void assignGroup (const BoxDirection_t & dirs,
		    const ScalorType & refP,
		    const ScalorType & beta);
  void clearGroups ();
  void calScale (const ScalorType * nowP,
		 ScalorType * scale) const;
};

#endif
