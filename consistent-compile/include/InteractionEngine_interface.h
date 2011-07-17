/**
 * @file   InteractionEngine_interface.h
 * @author Han Wang <han_wang@math.pku.edu.cn>
 * @date   Thu Jun 30 01:24:49 2011
 * 
 * @brief  Description of the class InteractionEngine_interface
 * 
 * 
 */

#ifndef __InteractionEngine_interface_h_wanghan__
#define __InteractionEngine_interface_h_wanghan__

#include "InteractionEngine.h"
#include "NeighborList_interface.h"
#include "Statistic_interface.h"
#include "MDSystem_interface.h"
#include "SystemNonBondedInteraction.h"
#include "SystemBondedInteraction.h"
#include "BondedInteractionList.h"
#include "WidomTestParticleInsertion.h"
#include "TwinRangeCorrectionRecorder.h"
#include "ExclusionList.h"

using namespace RectangularBoxGeometry;

/// InteractionEngine help applying both bonded and non-bonded interactions to atoms in the MDSystem.

/**
 * The bonded interaction is applied
 * (InteractionEngine::applyBondedInteraction) with the help of
 * BondedInteractionList.  Three methods are implemented to apply
 * non-bonded interaction
 * (InteractionEngine::applyNonBondedInteraction): all-pair comparing,
 * cell list and neighbor list. The cell list method makes use of
 * CellList, while the neighbor list method makes use of
 * NeighborList. All interactions are \b added to the force of
 * atoms. Therefore, at every MD time step, one should first run
 * InteractionEngine::clearInteraction to clear the interactions of
 * the last time step and then run
 * InteractionEngine::applyNonBondedInteraction and
 * InteractionEngine::applyBondedInteraction subsequently to apply new
 * interactions if necessary.
 * 
 */

class InteractionEngine
{
  dim3 atomGridDim;
  dim3 myBlockDim;
  IndexType calBondInteraction_sbuffSize;
  IndexType calAngleInteraction_sbuffSize;
  bool hasBond;
  bool hasAngle;
  bool hasExclusion;
  // IndexType applyNonBondedInteraction_CellList_sbuffSize;
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
  SumVector<ScalorType> sum_angle_vxx;
  SumVector<ScalorType> sum_angle_vyy;
  SumVector<ScalorType> sum_angle_vzz;
  cudaStream_t sum_stream[8];
  MDError err;
  ScalorType energyCorr;
  ScalorType pressureCorr;
private:
  IndexType maxNumExclusion;
  bool sharedExclusionList;
  size_t exclusion_sbuffSize;
private:
  void initNonBondedInteraction (const MDSystem & sys);
public:
  /** 
   * Constructor
   * 
   * @param sys Initialize a object of InteractionEngine for sys.
   * @param NThread Number of threads in a block.
   */
  InteractionEngine(const MDSystem  & sys, 
		    const IndexType & NThread)
      {init (sys, NThread);}
  /** 
   * Destructor
   * 
   */
  ~InteractionEngine();
  /** 
   * Initialize the class
   * 
   * @param sys Initialize a object of InteractionEngine for sys.
   * @param NThread Number of threads in a block.
   */
  void init (const MDSystem  & sys, 
	     const IndexType & NThread);
  /** 
   * Register non-bonded interaction.
   * 
   * @param sysNbInter Non-bonded interaction of the system.
   */
  void registNonBondedInteraction (const SystemNonBondedInteraction & sysNbInter);
  /** 
   * Register bonded interaction.
   * 
   * @param sysNbInter Bonded interaction of the system.
   */
  void registBondedInteraction    (const SystemBondedInteraction    & sysBdInter);
public:
  /** 
   * Set all interactions to zero.
   * 
   * @param sys 
   */
  void clearInteraction (MDSystem & sys);
  // apply nonbonded interaction and make nlist
  // void applyNonBondedInteraction (MDSystem & sys,
  // 				  const CellList & clist,
  // 				  const ScalorType & rcut,
  // 				  NeighborList & nlist,
  // 				  MDTimer *timer = NULL);
  // void applyNonBondedInteraction (MDSystem & sys,
  // 				  const CellList & clist,
  // 				  const ScalorType & rcut,
  // 				  NeighborList & nlist,
  // 				  MDStatistic & st,
  // 				  MDTimer *timer = NULL);
  /** 
   * Apply non-bonded interaction. All-parir compare method.
   * 
   * @param sys System to apply, which should be the same as the one
   * passed to InteractionEngine::init.
   * @param rcut Cut-off raidus for non-bonded interaction.
   * @param timer Timer measuring the performance.
   */
  void applyNonBondedInteraction (MDSystem & sys,
				  const ScalorType & rcut,
				  // const ExclusionList * excllist = NULL,
				  MDTimer *timer = NULL);
  /** 
   * Apply non-bonded interaction with statistic. All-parir compare method.
   * 
   * @param sys System to apply, which should be the same as the one
   * passed to InteractionEngine::init.
   * @param rcut Cut-off raidus for non-bonded interaction.
   * @param st Statistic.
   * @param timer Timer measuring the performance.
   */
  void applyNonBondedInteraction (MDSystem & sys,
				  const ScalorType & rcut,
				  MDStatistic & st,
				  // const ExclusionList * excllist = NULL,
				  MDTimer *timer = NULL);
  /** 
   * Apply non-bonded interaction. Cell list method.
   * 
   * @param sys System to apply, which should be the same as the one
   * passed to InteractionEngine::init.
   * @param clist Cell list
   * @param rcut Cut-off raidus
   * @param timer Timer measuring the performance.
   */
  void applyNonBondedInteraction (MDSystem & sys,
				  const CellList & clist,
				  const ScalorType & rcut,
				  // const ExclusionList * excllist = NULL,
				  MDTimer *timer = NULL);
  /** 
   * Apply non-bonded interaction with statistic. Cell list method.
   * 
   * @param sys System to apply, which should be the same as the one
   * passed to InteractionEngine::init.
   * @param clist Cell list
   * @param rcut Cut-off raidus
   * @param st Statistic.
   * @param timer Timer measuring the performance.
   */
  void applyNonBondedInteraction (MDSystem & sys,
				  const CellList & clist,
				  const ScalorType & rcut,
				  MDStatistic & st,
				  // const ExclusionList * excllist = NULL,
				  MDTimer *timer = NULL);
  /** 
   * Apply non-bonded interaction. Neighbor list method.
   * 
   * @param sys System to apply, which should be the same as the one
   * passed to InteractionEngine::init.
   * @param nlist Neighbor list
   * @param excllist Exclusion list
   * @param timer Timer measuring the performance.
   */
  void applyNonBondedInteraction (MDSystem & sys,
				  const NeighborList & nlist,
				  const ExclusionList * excllist = NULL,
				  MDTimer *timer = NULL);  
  /** 
   * Apply non-bonded interaction with statistic. Neighbor list method.
   * 
   * @param sys System to apply, which should be the same as the one
   * passed to InteractionEngine::init.
   * @param nlist Neighbor list
   * @param st Statistic.
   * @param excllist Exclusion list
   * @param timer Timer measuring the performance.
   */
  void applyNonBondedInteraction (MDSystem & sys,
				  const NeighborList & nlist,
				  MDStatistic & st,
				  const ExclusionList * excllist = NULL,
				  MDTimer *timer = NULL);
  
  void calTwinRangeCorrection (const MDSystem &			sys,
			       const CellList &			clist,
			       const ScalorType &		rcut1,
			       const ScalorType &		rcut2,
			       TwinRangeCorrectionRecorder &	twrec,
			       // const ExclusionList *		excllist = NULL,
			       MDTimer *			timer = NULL);
  void buildNeighborListCalTwinRangeCorrection (const MDSystem &	sys,
						const CellList &	clist,
						const ScalorType &	rcut1,
						const ScalorType &	rcut2,
						NeighborList &		nlist,
						TwinRangeCorrectionRecorder & twrec,
						// const ExclusionList *	excllist = NULL,
						MDTimer *		timer = NULL);
			       
  void calculateWidomDeltaEnergy (const MDSystem & sys,
				  const NeighborList & nlist,
				  WidomTestParticleInsertion_NVT & wtest,
				  MDTimer * timer = NULL);
  void calculateWidomDeltaEnergy (const MDSystem & sys,
				  const NeighborList & nlist,
				  WidomTestParticleInsertion_NVT2 & wtest,
				  MDTimer * timer = NULL);
  void calculateWidomDeltaEnergy (const MDSystem & sys,
				  const NeighborList & nlist,
				  WidomTestParticleInsertion_NPT & wtest,
				  MDTimer * timer = NULL);
  // void applyNonBondedInteractionCell  (MDSystem & sys,
  // 				       const NeighborList & nlist,
  // 				       MDTimer *timer );
  // void applyNonBondedInteractionCell  (MDSystem & sys,
  // 				       const NeighborList & nlist,
  // 				       MDStatistic & st,
  // 				       MDTimer *timer );
  /** 
   * Apply bonded interaction.
   * 
   * @param sys System to apply, which should be the same as the one
   * passed to InteractionEngine::init.   
   * @param bdlist Bonded interaction list.
   * @param timer Timer measuring the performance.
   */
  void applyBondedInteraction (MDSystem & sys,
			       const BondedInteractionList & bdlist,
			       MDTimer *timer = NULL);
  /** 
   * Apply bonded interaction with statistic.
   * 
   * @param sys System to apply, which should be the same as the one
   * passed to InteractionEngine::init.   
   * @param bdlist Bonded interaction list.
   * @param st Statistic.
   * @param timer Timer measuring the performance.
   */
  void applyBondedInteraction (MDSystem & sys,
			       const BondedInteractionList & bdlist,
			       MDStatistic & st,
			       MDTimer *timer = NULL);
}
    ;



#endif
