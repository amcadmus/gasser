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
  virtual void calCouple (const ScalorType * nowP,
			  ScalorType * lambda,
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
  void calCouple (const ScalorType * nowP,
		  ScalorType * lambda,
		  RectangularBox & box) const;
};


#ifdef DEVICE_CODE

#include "SumVector.h"
#include "MDSystem_interface.h"
#include "Statistic_interface.h"
#include "MDTimer_interface.h"

class NoseHoover_Chains2_Isobaric 
{
  dim3 atomGridDim;
  dim3 myBlockDim;
  IndexType sharedBuffSize;
// private:
public:
  ScalorType ep;
  ScalorType vep;
  // ScalorType NN;
  ScalorType Nf;
  ScalorType refP;
  ScalorType WW;

  // ScalorType xi1, xi2;
  // ScalorType vxi1, vxi2;
  ScalorType xi1;
  ScalorType vxi1;
  ScalorType Q1, Q2;
  ScalorType refT;
  // ScalorType LL;

  MDStatistic tmp_st;
  mdStatisticItem_t virial_array[3];
  mdStatisticItem_t kinetic_array[3];
  SumVector<ScalorType> sum_kxx;
  SumVector<ScalorType> sum_kyy;
  SumVector<ScalorType> sum_kzz;
public:
  void reinit (const MDSystem &sys,
	       const IndexType & NThread,
	       const ScalorType & ref_T,
	       const ScalorType & tau_T,
	       const ScalorType & ref_P,
	       const ScalorType & tau_P,
	       const ScalorType & pressureCorr);  
public:
  void operator_L_box (const ScalorType & dt,
  		       RectangularBox & box,
		       MDTimer * timer = NULL);
  void operator_L_ep  (const ScalorType & dt);
  void operator_L_r   (const ScalorType & dt,
		       MDSystem & sys,
		       MDTimer * timer = NULL);
  // output kinetic energy
  void operator_L_v (const ScalorType & dt,
		     MDSystem & sys,
		     MDStatistic & output_st,
		     MDTimer * timer = NULL);
  void operator_L_v (const ScalorType & dt,
		     MDSystem & sys,
		     MDTimer * timer = NULL);
  void operator_L_xi  (const ScalorType & dt);
  // output kinetic energy
  void operator_L_Cv (const ScalorType & dt,
		      MDSystem & sys,
		      MDStatistic & output_st);
  // input viral and kinetic energy
  void operator_L_Gep (const ScalorType & dt,
		       const MDSystem & sys,
		       const MDStatistic & input_st);
  void operator_L_vep (const ScalorType & dt);
  // input kinetic energy
  void operator_L_G1 (const ScalorType & dt,
		      const MDSystem & sys,
		      const MDStatistic & input_st);
  void operator_L_vxi1 (const ScalorType & dt);
  void operator_L_G2 (const ScalorType & dt);

  // input kinetic energy and viral
  // output kinetic energy
  void operator_L_CP (const ScalorType & dt,
		      MDSystem & sys,
		      const MDStatistic & input_st,
		      MDStatistic & output_st,
		      MDTimer * timer = NULL);
  // input kinetic energy and viral
  void operator_L_CP (const ScalorType & dt,
		      MDSystem & sys,
		      const MDStatistic & input_st,
		      MDTimer * timer = NULL);

  ScalorType HamiltonianContribution (const RectangularBox & box) const;
}
    ;

#endif


		       

#endif

