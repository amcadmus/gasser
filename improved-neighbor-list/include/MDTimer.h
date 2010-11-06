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
  mdTimeJudgeRebuild		= 19,
  mdTimeDataIO			= 20
};

typedef enum mdTimeItem mdTimeItem_t;

#endif
