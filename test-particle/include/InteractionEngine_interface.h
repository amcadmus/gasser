#ifndef __InteractionEngine_interface_h_wanghan__
#define __InteractionEngine_interface_h_wanghan__

#include "InteractionEngine.h"
#include "NeighborList_interface.h"
#include "Statistic_interface.h"
#include "MDSystem_interface.h"
#include "SystemNonBondedInteraction.h"
#include "SystemBondedInteraction.h"
#include "BondedInteractionList.h"

using namespace RectangularBoxGeometry;

class InteractionEngine_interface
{
  dim3 atomGridDim;
  dim3 myBlockDim;
  IndexType calBondInteraction_sbuffSize;
  IndexType calAngleInteraction_sbuffSize;
  bool hasBond;
  bool hasAngle;
  IndexType applyNonBondedInteraction_CellList_sbuffSize;
  // ScalorType * statistic_nb_buff0;
  // ScalorType * statistic_nb_buff1;
  // ScalorType * statistic_nb_buff2;
  // ScalorType * statistic_nb_buff3;  
  // ScalorType * statistic_b_buff0;
  // ScalorType * statistic_b_buff1;
  // ScalorType * statistic_b_buff2;
  // ScalorType * statistic_b_buff3;  
  SumVector<ScalorType> sum_nb_p;
  SumVector<ScalorType> sum_nb_vxx;
  SumVector<ScalorType> sum_nb_vyy;
  SumVector<ScalorType> sum_nb_vzz;
  SumVector<ScalorType> sum_b_p;
  SumVector<ScalorType> sum_b_vxx;
  SumVector<ScalorType> sum_b_vyy;
  SumVector<ScalorType> sum_b_vzz;
  SumVector<ScalorType> sum_angle_p;
  cudaStream_t sum_stream[8];
  MDError err;
private:
  void initNonBondedInteraction (const MDSystem & sys);
public:
  InteractionEngine_interface(const MDSystem  & sys, 
			      const IndexType & NThread)
      {init (sys, NThread);}
  ~InteractionEngine_interface();
  void init (const MDSystem  & sys, 
	     const IndexType & NThread);
  void registNonBondedInteraction (const SystemNonBondedInteraction & sysNbInter);
  void registBondedInteraction    (const SystemBondedInteraction    & sysBdInter);
public:
  void clearInteraction (MDSystem & sys);
  void applyNonBondedInteraction (MDSystem & sys,
				  const NeighborList & nlist,
				  MDTimer *timer = NULL);
  void applyNonBondedInteraction (MDSystem & sys,
				  const NeighborList & nlist,
				  MDStatistic & st,
				  MDTimer *timer = NULL);
  void applyNonBondedInteractionCell  (MDSystem & sys,
				       const NeighborList & nlist,
				       MDTimer *timer );
  void applyNonBondedInteractionCell  (MDSystem & sys,
				       const NeighborList & nlist,
				       MDStatistic & st,
				       MDTimer *timer );
  void applyBondedInteraction (MDSystem & sys,
			       const BondedInteractionList & bdlist,
			       MDTimer *timer = NULL);
  void applyBondedInteraction (MDSystem & sys,
			       const BondedInteractionList & bdlist,
			       MDStatistic & st,
			       MDTimer *timer = NULL);
}
    ;



#endif