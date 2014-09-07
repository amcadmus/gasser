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

class Barostat_VCouple : public Barostat
{
public:
  // virtual void calCouple (const ScalorType * nowP,
  // 			  ScalorType * lambda,
  // 			  RectangularBox & box) const = 0;
  virtual void getCouple (ScalorType * lambda,
			  RectangularBox & box) const = 0;
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

class Barostat_ParrinelloRahman : public Barostat_VCouple
{
  mutable ScalorType b[3];
  mutable ScalorType c[3];
  
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
  Barostat_ParrinelloRahman ();
  void reinit (const ScalorType & dt,
	       const ScalorType & tau_P,
	       const RectangularBox & box);
  void assignGroup (const BoxDirection_t & dirs,
		    const ScalorType & refP,
		    const ScalorType & beta);
  void clearGroups ();
  
  // void calCouple (const ScalorType * nowP,
  // 		  ScalorType * lambda,
  // 		  RectangularBox & box) const;
  void integrate_LeapFrog (const ScalorType * nowP);
  void calCouple (ScalorType * lambda,
		  RectangularBox & box) const = 0;

};


#endif

