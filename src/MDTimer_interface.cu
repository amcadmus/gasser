#include "MDTimer_interface.h"
#include "common.h"

MDTimer::MDTimer()
{
  for (IndexType i = 0; i < NumberOfMemberInMdTime; ++i){
    cudaEventCreate (&start[i]);
    cudaEventCreate (&stop [i]);
    timeRecord[i] = 0.f;
  }
}

MDTimer::~MDTimer()
{
  for (IndexType i = 0; i < NumberOfMemberInMdTime; ++i){
    cudaEventDestroy (start[i]);
    cudaEventDestroy (stop [i]);
  }
}

void MDTimer::printRecord (FILE * fp)
{
  fprintf(fp, "Normalize system:                    %1.3e s   %3.1f %\n",
	  0.001 * timeRecord[mdTimeNormalizeSys],
	  100 * timeRecord[mdTimeNormalizeSys] / timeRecord[mdTimeTotal]);
  fprintf(fp, "Build cell list:                     %1.3e s   %3.1f %\n",
	  0.001 * timeRecord[mdTimeBuildCellList],
	  100 * timeRecord[mdTimeBuildCellList] / timeRecord[mdTimeTotal]);
  fprintf(fp, "Build neighbor list:                 %1.3e s   %3.1f %\n",
	  0.001 * timeRecord[mdTimeBuildNeighborList],
	  100 * timeRecord[mdTimeBuildNeighborList] / timeRecord[mdTimeTotal]);
  fprintf(fp, "Judge rebuild:                       %1.3e s   %3.1f %\n",
	  0.001 * timeRecord[mdTimeJudgeRebuild],
	  100 * timeRecord[mdTimeJudgeRebuild] / timeRecord[mdTimeTotal]);
  fprintf(fp, "Reshuffle system:                    %1.3e s   %3.1f %\n",
	  0.001 * timeRecord[mdTimeReshuffleSystem],
	  100 * timeRecord[mdTimeReshuffleSystem] / timeRecord[mdTimeTotal]);
  fprintf(fp, "Integrate:                           %1.3e s   %3.1f %\n",
	  0.001 * timeRecord[mdTimeIntegrator],
	  100 * timeRecord[mdTimeIntegrator] / timeRecord[mdTimeTotal]);
  fprintf(fp, "Remove translational freedom:        %1.3e s   %3.1f %\n",
	  0.001 * timeRecord[mdTimeRemoveTransFreedom],
	  100 * timeRecord[mdTimeRemoveTransFreedom] / timeRecord[mdTimeTotal]);
  fprintf(fp, "Cal. non-bonded interaction:         %1.3e s   %3.1f %\n",
	  0.001 * timeRecord[mdTimeNonBondedInteraction],
	  100 * timeRecord[mdTimeNonBondedInteraction] / timeRecord[mdTimeTotal]);
  fprintf(fp, "Cal. non-bonded interaction(st):     %1.3e s   %3.1f %\n",
	  0.001 * timeRecord[mdTimeNBInterStatistic],
	  100 * timeRecord[mdTimeNBInterStatistic] / timeRecord[mdTimeTotal]);
  fprintf(fp, "Cal. non-bonded twin range corr:     %1.3e s   %3.1f %\n",
	  0.001 * timeRecord[mdTimeNBInterTwinRange],
	  100 * timeRecord[mdTimeNBInterTwinRange] / timeRecord[mdTimeTotal]);
  fprintf(fp, "Cal. bond interaction:               %1.3e s   %3.1f %\n",
	  0.001 * timeRecord[mdTimeBondedInteraction],
	  100 * timeRecord[mdTimeBondedInteraction] / timeRecord[mdTimeTotal]);
  fprintf(fp, "Cal. bond interaction(st):           %1.3e s   %3.1f %\n",
	  0.001 * timeRecord[mdTimeBInterStatistic],
	  100 * timeRecord[mdTimeBInterStatistic] / timeRecord[mdTimeTotal]);
  fprintf(fp, "Cal. angle interaction:              %1.3e s   %3.1f %\n",
	  0.001 * timeRecord[mdTimeAngleInteraction],
	  100 * timeRecord[mdTimeAngleInteraction] / timeRecord[mdTimeTotal]);
  fprintf(fp, "Cal. angle interaction(st):          %1.3e s   %3.1f %\n",
	  0.001 * timeRecord[mdTimeAngleInterStatistic],
	  100 * timeRecord[mdTimeAngleInterStatistic] / timeRecord[mdTimeTotal]);
  fprintf(fp, "SPME cal Q: build mesh nlist:        %1.3e s   %3.1f %\n",
	  0.001 * timeRecord[mdTimeSPMERecMeshNeighborList],
	  100 * timeRecord[mdTimeSPMERecMeshNeighborList] / timeRecord[mdTimeTotal]);
  fprintf(fp, "SPME cal Q: from mesh nlist:         %1.3e s   %3.1f %\n",
	  0.001 * timeRecord[mdTimeSPMECalQFromNList],
	  100 * timeRecord[mdTimeSPMECalQFromNList] / timeRecord[mdTimeTotal]);
  fprintf(fp, "SPME cal Q:                          %1.3e s   %3.1f %\n",
	  0.001 * timeRecord[mdTimeSPMERecCalQ],
	  100 * timeRecord[mdTimeSPMERecCalQ] / timeRecord[mdTimeTotal]);
  fprintf(fp, "SPME FFT:                            %1.3e s   %3.1f %\n",
	  0.001 * timeRecord[mdTimeSPMERecFFT],
	  100 * timeRecord[mdTimeSPMERecFFT] / timeRecord[mdTimeTotal]);
  fprintf(fp, "SPME time matrix:                    %1.3e s   %3.1f %\n",
	  0.001 * timeRecord[mdTimeSPMERecTimeMatrix],
	  100 * timeRecord[mdTimeSPMERecTimeMatrix] / timeRecord[mdTimeTotal]);
  fprintf(fp, "SPME cal force:                      %1.3e s   %3.1f %\n",
	  0.001 * timeRecord[mdTimeSPMERecForce],
	  100 * timeRecord[mdTimeSPMERecForce] / timeRecord[mdTimeTotal]);
  fprintf(fp, "SPME cal energy:                     %1.3e s   %3.1f %\n",
	  0.001 * timeRecord[mdTimeSPMERecEnergy],
	  100 * timeRecord[mdTimeSPMERecEnergy] / timeRecord[mdTimeTotal]);
  fprintf(fp, "Data transfer:                       %1.3e s   %3.1f %\n",
	  0.001 * timeRecord[mdTimeDataTransfer],
	  100 * timeRecord[mdTimeDataTransfer] / timeRecord[mdTimeTotal]);
  fprintf(fp, "Data IO:                             %1.3e s   %3.1f %\n",
	  0.001 * timeRecord[mdTimeDataIO],
	  100 * timeRecord[mdTimeDataIO] / timeRecord[mdTimeTotal]);
  TimeType otherTime = TimeType(0.);
  for (IndexType i = 1; i < NumberOfMemberInMdTime; ++i){
    otherTime += timeRecord[i];
  }
  otherTime = timeRecord[mdTimeTotal] - otherTime;
  fprintf(fp, "Others uncounted:                    %1.3e s   %3.1f %\n",
	  0.001 * otherTime,
	  100 * otherTime / timeRecord[mdTimeTotal]);
  fprintf(fp, "Total time:                          %1.3e s   %3.1f %\n",
	  0.001 * timeRecord[mdTimeTotal],
	  100 * timeRecord[mdTimeTotal] / timeRecord[mdTimeTotal]);
}

