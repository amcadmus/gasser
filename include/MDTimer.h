#ifndef __MDTimer_h_wanghan__
#define __MDTimer_h_wanghan__

#define NumberOfMemberInMdTime 32

typedef double		TimeType;


enum mdTimeItem {
  mdTimeTotal			= 0,
  mdTimeNormalizeSys		= 1,
  mdTimeBuildCellList		= 2,
  mdTimeBuildNeighborList	= 3,
  mdTimeReshuffleSystem		= 4,
  mdTimeIntegrator		= 5,
  mdTimeRemoveTransFreedom	= 6,
  mdTimeNonBondedInteraction	= 7,
  mdTimeNBInterStatistic	= 8,
  mdTimeNBInterTwinRange	= 14,
  mdTimeBondedInteraction	= 9,
  mdTimeBInterStatistic		= 10,
  mdTimeAngleInteraction	= 11,
  mdTimeAngleInterStatistic	= 12,
  mdTimeDataTransfer		= 13,
  mdTimeSPMERecMeshNeighborList	= 20,
  mdTimeSPMECalQFromNList	= 21,
  mdTimeSPMERecCalQ		= 15,
  mdTimeSPMERecFFT		= 16,
  mdTimeSPMERecTimeMatrix	= 17,
  mdTimeSPMERecForce		= 18,
  mdTimeSPMERecEnergy		= 19,
  mdTimeJudgeRebuild		= 29,
  mdTimeDataIO			= 30
};

typedef enum mdTimeItem mdTimeItem_t;

#endif
