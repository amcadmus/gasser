#define DEVICE_CODE

#include "systemDefines.h"
#include "InteractionEngine_interface.h"
#include "NonBondedInteraction.h"
#include "BondInteraction.h"
#include "AngleInteraction.h"
// #include "CellList_interface.h"
#include "Auxiliary.h"

texture<CoordType,  1, cudaReadModeElementType> global_texRef_interaction_coord;
texture<TypeType ,  1, cudaReadModeElementType> global_texRef_interaction_type;

__constant__
InteractionType nonBondedInteractionType [MaxNumberNonBondedInteraction];
__constant__
ScalorType nonBondedInteractionParameter [MaxNumberNonBondedInteractionParameter];
__constant__
IndexType nonBondedInteractionParameterPosition [MaxNumberNonBondedInteraction];
__constant__
InteractionType bondedInteractionType [MaxNumberBondedInteraction];
__constant__
IndexType bondedInteractionParameterPosition [MaxNumberBondedInteraction];
__constant__
ScalorType bondedInteractionParameter [MaxNumberBondedInteractionParamemter];
__constant__
IndexType const_nonBondedInteractionTableLength[1];
__constant__
IndexType const_numAtomType[1];
__constant__
IndexType const_nonBondedInteractionTable [MaxLengthNonBondedInteractionTable];


void InteractionEngine::init (const MDSystem  & sys,
			      const IndexType & NTread)
{
  hasBond = false;
  hasAngle = false;
  myBlockDim.y = 1;
  myBlockDim.z = 1;
  myBlockDim.x = NTread;
  IndexType nob;
  if (sys.ddata.numAtom % myBlockDim.x == 0){
    nob = sys.ddata.numAtom / myBlockDim.x;
  } else {
    nob = sys.ddata.numAtom / myBlockDim.x + 1;
  }
  atomGridDim = toGridDim (nob);

  // size_t sizetype = sizeof(TypeType)*sys.ddata.numMem;
  cudaBindTexture(0, global_texRef_interaction_coord, sys.ddata.coord,
		  sizeof(CoordType) * sys.ddata.numMem);
  cudaBindTexture(0, global_texRef_interaction_type, sys.ddata.type,
		  sizeof(TypeType)  * sys.ddata.numMem);
  checkCUDAError ("InteractionEngine::init, bind texture");
  
  // init sum vectors
  sum_nb_p.reinit (sys.ddata.numAtom, NThreadForSum);
  sum_nb_vxx.reinit (sys.ddata.numAtom, NThreadForSum);
  sum_nb_vyy.reinit (sys.ddata.numAtom, NThreadForSum);
  sum_nb_vzz.reinit (sys.ddata.numAtom, NThreadForSum);
  sum_b_p.reinit (nob, NThreadForSum);
  sum_b_vxx.reinit (nob, NThreadForSum);
  sum_b_vyy.reinit (nob, NThreadForSum);
  sum_b_vzz.reinit (nob, NThreadForSum);
  sum_angle_p.reinit (nob, NThreadForSum);
  sum_angle_vxx.reinit (nob, NThreadForSum);
  sum_angle_vyy.reinit (nob, NThreadForSum);
  sum_angle_vzz.reinit (nob, NThreadForSum);
  for (IndexType i = 0; i < 8; ++i){
    cudaStreamCreate(&sum_stream[i]);
  }
  checkCUDAError ("InteractionEngine::init init sum statistic");

  // exclusion list
  maxNumExclusion = 0;
  sharedExclusionList = false;
  exclusion_sbuffSize = size_t(0);
}


static IndexType hroundUp4 (const IndexType x)
{
  if (x & 3 == 0){
    return x;
  }
  else {
    return ((x >> 2) + 1) << 2;
  }
}

void InteractionEngine::
registNonBondedInteraction (const SystemNonBondedInteraction & sysNbInter)
{
  if (! sysNbInter.beBuilt()) {
    throw MDExcptUnbuiltNonBondedInteraction ("InteractionEngine");
  }
  if (sysNbInter.numberOfInteraction() > MaxNumberBondedInteraction ){
    throw MDExcptExceedConstantMemLimit (
	"InteractionEngine::registNonBondedInteraction",
	"nonBonedInteractionType",
	MaxNumberNonBondedInteraction * sizeof(InteractionType));
  }
  if (sysNbInter.numberOfParameter() > MaxNumberNonBondedInteractionParameter ){
    throw MDExcptExceedConstantMemLimit (
	"InteractionEngine::registNonBondedInteraction",
	"nonBondedInteractionParameter",
	MaxNumberNonBondedInteractionParameter * sizeof(ScalorType));
  }

  cudaMemcpyToSymbol (nonBondedInteractionType,
		      sysNbInter.interactionType(), 
  		      sizeof(InteractionType) * sysNbInter.numberOfInteraction());
  cudaMemcpyToSymbol (nonBondedInteractionParameterPosition,
		      sysNbInter.interactionParameterPosition(),
  		      sizeof(ScalorType) * sysNbInter.numberOfInteraction());
  cudaMemcpyToSymbol (nonBondedInteractionParameter,
		      sysNbInter.interactionParameter(),
		      sizeof(IndexType) * sysNbInter.numberOfParameter());
  checkCUDAError ("InteractionEngine::init, init NB force setting");

  IndexType tableSize = sysNbInter.interactionTableSize();
  IndexType tmpNumAtomType = sysNbInter.numberOfAtomTypes();
  if (tableSize > MaxLengthNonBondedInteractionTable){
    throw MDExcptExceedConstantMemLimit(
	"InteractionEngine::registNonBondedInteraction",
	"nonBondedInteractionTable",
	MaxLengthNonBondedInteractionTable * sizeof (ScalorType));
  }
  cudaMemcpyToSymbol (const_nonBondedInteractionTableLength,
  		      &tableSize,
  		      sizeof (IndexType));
  checkCUDAError ("InteractionEngine::init, const_nonBondedInteractionTableLength");
  cudaMemcpyToSymbol (const_numAtomType,
		      &tmpNumAtomType,
		      sizeof (IndexType));
  checkCUDAError ("InteractionEngine::init, const_numAtomType");
  cudaMemcpyToSymbol (const_nonBondedInteractionTable,
  		      sysNbInter.interactionTable(),
  		      sizeof (IndexType) * tableSize);
  checkCUDAError ("InteractionEngine::init, const_nonBondedInteractionTable");

  // applyNonBondedInteraction_CellList_sbuffSize =
  //     sizeof(IndexType) *	hroundUp4(myBlockDim.x) +
  //     sizeof(CoordType) *	hroundUp4(myBlockDim.x) +
  //     sizeof(TypeType)  *	hroundUp4(myBlockDim.x);
  // printf ("total %d\npart1 %d\npart2 %d\npart3 %d\nround %d\n",
  // 	  applyNonBondedInteraction_CellList_sbuffSize,
  // 	  sizeof(IndexType) *	hroundUp4(myBlockDim.x),
  // 	  sizeof(CoordType) *	hroundUp4(myBlockDim.x),
  // 	  sizeof(TypeType)  *	hroundUp4(myBlockDim.x),
  // 	  hroundUp4(myBlockDim.x));
  // checkCUDAError ("InteractionEngine::init, init nonBondedInteractionTable");

  energyCorr = sysNbInter.energyCorrection ();
  pressureCorr = sysNbInter.pressureCorrection ();

  maxNumExclusion = sysNbInter.maxNumberOfExclusion();
  if (maxNumExclusion != 0){
    sharedExclusionList = true;
    exclusion_sbuffSize = myBlockDim.x * maxNumExclusion * sizeof(IndexType);
    if (exclusion_sbuffSize > SystemSharedBuffSize){
      sharedExclusionList = false;
    }
  }
}


void InteractionEngine::
registBondedInteraction (const SystemBondedInteraction & sysBdInter)
{
  if (sysBdInter.hasBond() ){
    hasBond = true;
  }
  if (sysBdInter.hasAngle()){
    hasAngle = true;
  }

  if (sysBdInter.numberOfInteraction() > MaxNumberBondedInteraction ){
    throw MDExcptExceedConstantMemLimit (
	"InteractionEngine::registBondedInteraction",
	"bondedInteractionType",
	MaxNumberBondedInteraction * sizeof(InteractionType));
  }
  if (sysBdInter.numberOfParameter() > MaxNumberBondedInteractionParamemter ){
    throw MDExcptExceedConstantMemLimit (
	"InteractionEngine::registBondedInteraction",
	"bondedInteractionParameter",
	MaxNumberBondedInteractionParamemter * sizeof(ScalorType));
  }

  if (hasBond || hasAngle){
    cudaMemcpyToSymbol (bondedInteractionType,
			sysBdInter.interactionType(),
			sizeof(InteractionType) * sysBdInter.numberOfInteraction());
    cudaMemcpyToSymbol (bondedInteractionParameterPosition,
			sysBdInter.interactionParameterPosition(),
			sizeof(ScalorType) * sysBdInter.numberOfInteraction());
    cudaMemcpyToSymbol (bondedInteractionParameter,
			sysBdInter.interactionParameter(),
			sizeof(IndexType) * sysBdInter.numberOfParameter());
    checkCUDAError ("InteractionEngine::init, init bond force setting");
    // cal shared buff size
    calBondInteraction_sbuffSize  = myBlockDim.x * sizeof(ScalorType);
    calAngleInteraction_sbuffSize = myBlockDim.x * sizeof(ScalorType);
  }
}

InteractionEngine::~InteractionEngine()
{
  cudaUnbindTexture(global_texRef_interaction_coord);
  cudaUnbindTexture(global_texRef_interaction_type);
  for (IndexType i = 0; i < 8; ++i){
    cudaStreamDestroy(sum_stream[i]);
  }
}

void InteractionEngine::clearInteraction (MDSystem & sys)
{
  clearForce
      <<<atomGridDim, myBlockDim>>>(
	  sys.ddata.numAtom,
	  sys.ddata.forcx, sys.ddata.forcy, sys.ddata.forcz);
  checkCUDAError ("InteractionEngine::clearInteraction");
}


// nblock should be 1 and block size should be 1
__global__ void
applyEnergyPressureCorrection (ScalorType * ddata,
			       ScalorType energyCorr,
			       ScalorType pressureCorr)
{
  ddata[mdStatisticEnergyCorrection] = energyCorr;
  ddata[mdStatisticPressureCorrection] = pressureCorr;
}



void InteractionEngine::
applyNonBondedInteraction  (MDSystem & sys,
			    const NeighborList & nlist,
			    const ExclusionList * excllist,
			    MDTimer *timer )
{
  if (timer != NULL) timer->tic(mdTimeNonBondedInteraction);
  if (excllist == NULL){
  calNonBondedInteraction_neighbor
      <<<atomGridDim, myBlockDim>>> (
	  sys.ddata.numAtom,
	  sys.ddata.coord,
	  sys.ddata.forcx,
	  sys.ddata.forcy,
	  sys.ddata.forcz,
	  sys.ddata.type, 
	  sys.box,
	  nlist.dnlist);
  }
  else{
    calNonBondedInteraction_neighbor
	<<<atomGridDim, myBlockDim,
	exclusion_sbuffSize>>> (
	    sys.ddata.numAtom,
	    sys.ddata.coord,
	    sys.ddata.forcx,
	    sys.ddata.forcy,
	    sys.ddata.forcz,
	    sys.ddata.type, 
	    sys.box,
	    nlist.dnlist,
	    excllist->dexcllist,
	    sharedExclusionList
	    );
  }
  checkCUDAError ("InteractionEngine::applyInteraction nb");
  err.check ("interaction engine nb");	
  if (timer != NULL) timer->toc(mdTimeNonBondedInteraction);
}

// void InteractionEngine::
// applyNonBondedInteraction  (MDSystem & sys,
// 			    const CellList & clist,
// 			    const ScalorType & rcut,
// 			    NeighborList & nlist,
// 			    MDTimer *timer )
// {
//   if (timer != NULL) timer->tic(mdTimeBuildNeighborList);
//   size_t applyNonBondedInteraction_CellList_sbuffSize =
//       (sizeof(IndexType) + sizeof(CoordType) + sizeof(TypeType)) *
//       hroundUp4(clist.getBlockDim().x);
//       // sizeof(IndexType) *	hroundUp4(myBlockDim.x) +
//       // sizeof(CoordType) *	hroundUp4(myBlockDim.x) +
//       // sizeof(TypeType)  *	hroundUp4(myBlockDim.x);
//   calNonBondedInteraction
//       <<<clist.getCellGrimDim(), clist.getBlockDim(),
//       applyNonBondedInteraction_CellList_sbuffSize>>> (
// 	  sys.ddata.numAtom,
// 	  sys.ddata.coord,
// 	  sys.ddata.forcx,
// 	  sys.ddata.forcy,
// 	  sys.ddata.forcz,
// 	  sys.ddata.type, 
// 	  sys.box,
// 	  clist.dclist,
// 	  rcut,
// 	  nlist.dnlist,
// 	  err.ptr_de);
//   checkCUDAError ("InteractionEngine::applyInteraction nb");
//   err.check ("interaction engine nb");
//   if (timer != NULL) timer->toc(mdTimeBuildNeighborList);
// }



void InteractionEngine::
applyNonBondedInteraction (MDSystem & sys,
			   const NeighborList & nlist,
			   MDStatistic & st,
			   const ExclusionList * excllist,
			   MDTimer *timer )
{
  if (timer != NULL) timer->tic(mdTimeNBInterStatistic);
  if (excllist == NULL){
    calNonBondedInteraction_neighbor
	<<<atomGridDim, myBlockDim>>> (
	    sys.ddata.numAtom,
	    sys.ddata.coord,
	    sys.ddata.forcx,
	    sys.ddata.forcy,
	    sys.ddata.forcz,
	    sys.ddata.type, 
	    sys.box,
	    nlist.dnlist
	    ,
	    sum_nb_p.buff,
	    sum_nb_vxx.buff,
	    sum_nb_vyy.buff,
	    sum_nb_vzz.buff
	    );
  }
  else {
    calNonBondedInteraction_neighbor
	<<<atomGridDim, myBlockDim,
	exclusion_sbuffSize>>> (
	    sys.ddata.numAtom,
	    sys.ddata.coord,
	    sys.ddata.forcx,
	    sys.ddata.forcy,
	    sys.ddata.forcz,
	    sys.ddata.type, 
	    sys.box,
	    nlist.dnlist
	    ,
	    excllist->dexcllist,
	    sharedExclusionList,
	    sum_nb_p.buff,
	    sum_nb_vxx.buff,
	    sum_nb_vyy.buff,
	    sum_nb_vzz.buff
	    );
  }
  checkCUDAError ("InteractionEngine::applyInteraction nb (with statistic)");
  err.check ("interaction engine nb");	
  cudaThreadSynchronize();
  sum_nb_p.sumBuffAdd(st.ddata, mdStatisticNonBondedPotential, 0);
  sum_nb_vxx.sumBuffAdd(st.ddata, mdStatisticVirialXX);
  sum_nb_vyy.sumBuffAdd(st.ddata, mdStatisticVirialYY);
  sum_nb_vzz.sumBuffAdd(st.ddata, mdStatisticVirialZZ);
  ScalorType volumei = sys.box.size.x * sys.box.size.y * sys.box.size.z;
  volumei = 1.f / volumei;
  // printf ("apply Ec %f, Pc %f\n",
  // 	  energyCorr * volumei,
  // 	  pressureCorr * volumei * volumei);
  applyEnergyPressureCorrection
      <<<1, 1, 0>>> (st.ddata,
			energyCorr * volumei,
			pressureCorr * volumei * volumei);
  cudaThreadSynchronize();
  if (timer != NULL) timer->toc(mdTimeNBInterStatistic);
}


void InteractionEngine::
applyNonBondedInteraction  (MDSystem & sys,
			    const CellList & clist,
			    const ScalorType & rcut,
			    MDTimer *timer )
{
  if (!clist.isempty()){
    if (timer != NULL) timer->tic(mdTimeNonBondedInteraction);
    size_t applyNonBondedInteraction_CellList_sbuffSize =
	(sizeof(IndexType) + sizeof(CoordType) + sizeof(TypeType)) *
	hroundUp4(clist.getCellBlockDim().x);
    calNonBondedInteraction_cell
	<<<clist.getCellGrimDim(), clist.getCellBlockDim(),
	applyNonBondedInteraction_CellList_sbuffSize>>> (
	    sys.ddata.numAtom,
	    sys.ddata.coord,
	    sys.ddata.forcx,
	    sys.ddata.forcy,
	    sys.ddata.forcz,
	    sys.ddata.type, 
	    sys.box,
	    clist.dclist,
	    rcut,
	    err.ptr_de);
    checkCUDAError ("InteractionEngine::applyInteraction nb");
    err.check ("interaction engine nb");    
    if (timer != NULL) timer->toc(mdTimeNonBondedInteraction);
  }
  else {
    applyNonBondedInteraction (sys, rcut, timer);
  }
}


void InteractionEngine::
applyNonBondedInteraction (MDSystem & sys,
			   const CellList & clist,
			   const ScalorType & rcut,
			   MDStatistic & st,
			   MDTimer *timer )
{
  if (!clist.isempty()){
    if (timer != NULL) timer->tic(mdTimeNBInterStatistic);
    size_t applyNonBondedInteraction_CellList_sbuffSize =
	(sizeof(IndexType) + sizeof(CoordType) + sizeof(TypeType)) *
	hroundUp4(clist.getCellBlockDim().x);
    calNonBondedInteraction_cell
	<<<clist.getCellGrimDim(), clist.getCellBlockDim(),
	applyNonBondedInteraction_CellList_sbuffSize>>> (
	    sys.ddata.numAtom,
	    sys.ddata.coord,
	    sys.ddata.forcx,
	    sys.ddata.forcy,
	    sys.ddata.forcz,
	    sys.ddata.type, 
	    sys.box,
	    clist.dclist,
	    rcut,
	    sum_nb_p.buff,
	    sum_nb_vxx.buff,
	    sum_nb_vyy.buff,
	    sum_nb_vzz.buff,
	    err.ptr_de
	    );
    checkCUDAError ("InteractionEngine::applyInteraction nb (with statistic)");
    err.check ("interaction engine nb");	
    cudaThreadSynchronize();
    sum_nb_p.sumBuffAdd(st.ddata, mdStatisticNonBondedPotential, 0);
    sum_nb_vxx.sumBuffAdd(st.ddata, mdStatisticVirialXX);
    sum_nb_vyy.sumBuffAdd(st.ddata, mdStatisticVirialYY);
    sum_nb_vzz.sumBuffAdd(st.ddata, mdStatisticVirialZZ);
    ScalorType volumei = sys.box.size.x * sys.box.size.y * sys.box.size.z;
    volumei = 1.f / volumei;
    // printf ("apply Ec %f, Pc %f\n",
    // 	  energyCorr * volumei,
    // 	  pressureCorr * volumei * volumei);
    applyEnergyPressureCorrection
	<<<1, 1, 0>>> (st.ddata,
			  energyCorr * volumei,
			  pressureCorr * volumei * volumei);
    cudaThreadSynchronize();
    if (timer != NULL) timer->toc(mdTimeNBInterStatistic);
  }
  else {
    applyNonBondedInteraction (sys, rcut, st, timer);
  }
}


void InteractionEngine::
applyNonBondedInteraction  (MDSystem & sys,
			    const ScalorType & rcut,
			    MDTimer *timer )
{
  if (timer != NULL) timer->tic(mdTimeNonBondedInteraction);
  size_t applyNonBondedInteraction_AllPair_sbuffSize =
      (sizeof(CoordType) + sizeof(TypeType)) *
      hroundUp4(myBlockDim.x);
  calNonBondedInteraction_all
      <<<atomGridDim, myBlockDim,
      applyNonBondedInteraction_AllPair_sbuffSize>>> (
	  sys.ddata.numAtom,
	  sys.ddata.coord,
	  sys.ddata.forcx,
	  sys.ddata.forcy,
	  sys.ddata.forcz,
	  sys.ddata.type, 
	  sys.box,
	  rcut,
	  err.ptr_de);
  checkCUDAError ("InteractionEngine::applyInteraction nb");
  err.check ("interaction engine nb");    
  if (timer != NULL) timer->toc(mdTimeNonBondedInteraction);
}


void InteractionEngine::
applyNonBondedInteraction (MDSystem & sys,
			   const ScalorType & rcut,
			   MDStatistic & st,
			   MDTimer *timer )
{
  if (timer != NULL) timer->tic(mdTimeNBInterStatistic);
  size_t applyNonBondedInteraction_AllPair_sbuffSize =
      (sizeof(CoordType) + sizeof(TypeType)) *
      hroundUp4(myBlockDim.x);
  calNonBondedInteraction_all
      <<<atomGridDim, myBlockDim,
      applyNonBondedInteraction_AllPair_sbuffSize>>> (
	  sys.ddata.numAtom,
	  sys.ddata.coord,
	  sys.ddata.forcx,
	  sys.ddata.forcy,
	  sys.ddata.forcz,
	  sys.ddata.type, 
	  sys.box,
	  rcut,
	  sum_nb_p.buff,
	  sum_nb_vxx.buff,
	  sum_nb_vyy.buff,
	  sum_nb_vzz.buff,
	  err.ptr_de);
  checkCUDAError ("InteractionEngine::applyInteraction nb (with statistic)");
  err.check ("interaction engine nb");	
  cudaThreadSynchronize();
  sum_nb_p.sumBuffAdd(st.ddata, mdStatisticNonBondedPotential, 0);
  sum_nb_vxx.sumBuffAdd(st.ddata, mdStatisticVirialXX);
  sum_nb_vyy.sumBuffAdd(st.ddata, mdStatisticVirialYY);
  sum_nb_vzz.sumBuffAdd(st.ddata, mdStatisticVirialZZ);
  ScalorType volumei = sys.box.size.x * sys.box.size.y * sys.box.size.z;
  volumei = 1.f / volumei;
  // printf ("apply Ec %f, Pc %f\n",
  // 	  energyCorr * volumei,
  // 	  pressureCorr * volumei * volumei);
  applyEnergyPressureCorrection
      <<<1, 1, 0>>> (st.ddata,
			energyCorr * volumei,
			pressureCorr * volumei * volumei);
  cudaThreadSynchronize();
  if (timer != NULL) timer->toc(mdTimeNBInterStatistic);
}



// void InteractionEngine::
// applyNonBondedInteraction (MDSystem & sys,
// 			   const CellList & clist,
// 			   const ScalorType & rcut,
// 			   NeighborList & nlist,
// 			   MDStatistic & st,
// 			   MDTimer *timer )
// {
//   if (timer != NULL) timer->tic(mdTimeBuildNeighborList);
//   if (!clist.isempty()){
//     size_t applyNonBondedInteraction_CellList_sbuffSize =
// 	(sizeof(IndexType) + sizeof(CoordType) + sizeof(TypeType)) *
// 	hroundUp4(clist.getBlockDim().x);
//     calNonBondedInteraction
// 	<<<clist.getCellGrimDim(), clist.getBlockDim(),
// 	applyNonBondedInteraction_CellList_sbuffSize>>> (
// 	    sys.ddata.numAtom,
// 	    sys.ddata.coord,
// 	    sys.ddata.forcx,
// 	    sys.ddata.forcy,
// 	    sys.ddata.forcz,
// 	    sys.ddata.type, 
// 	    sys.box,
// 	    clist.dclist,
// 	    rcut,
// 	    nlist.dnlist,
// 	    sum_nb_p.buff,
// 	    sum_nb_vxx.buff,
// 	    sum_nb_vyy.buff,
// 	    sum_nb_vzz.buff,
// 	    err.ptr_de
// 	    );
//   }
//   checkCUDAError ("InteractionEngine::applyInteraction nb (with statistic)");
//   err.check ("interaction engine nb");	
//   cudaThreadSynchronize();
//   sum_nb_p.sumBuffAdd(st.ddata, mdStatisticNonBondedPotential, 0);
//   sum_nb_vxx.sumBuffAdd(st.ddata, mdStatisticVirialXX, 1);
//   sum_nb_vyy.sumBuffAdd(st.ddata, mdStatisticVirialYY, 2);
//   sum_nb_vzz.sumBuffAdd(st.ddata, mdStatisticVirialZZ, 3);
//   ScalorType volumei = sys.box.size.x * sys.box.size.y * sys.box.size.z;
//   volumei = 1.f / volumei;
//   // printf ("apply Ec %f, Pc %f\n",
//   // 	  energyCorr * volumei,
//   // 	  pressureCorr * volumei * volumei);
//   applyEnergyPressureCorrection
//       <<<1, 1, 0, 4>>> (st.ddata,
// 			energyCorr * volumei,
// 			pressureCorr * volumei * volumei);
//   cudaThreadSynchronize();
//   if (timer != NULL) timer->toc(mdTimeBuildNeighborList);
// }



void InteractionEngine::
calTwinRangeCorrection (const MDSystem &		sys,
			const CellList &		clist,
			const ScalorType &		rcut1,
			const ScalorType &		rcut2,
			TwinRangeCorrectionRecorder &	twrec,
			MDTimer *			timer)
{
  if (timer != NULL) timer->tic(mdTimeNBInterTwinRange);
  if (clist.isempty()){
    size_t applyNonBondedInteraction_AllPair_sbuffSize =
	(sizeof(CoordType) + sizeof(TypeType)) *
	hroundUp4(myBlockDim.x);
    calTwinRangeCorrection_all
	<<<atomGridDim, myBlockDim,
	applyNonBondedInteraction_AllPair_sbuffSize>>> (
	    sys.ddata.numAtom,
	    sys.ddata.coord,
	    twrec.forcx,
	    twrec.forcy,
	    twrec.forcz,
	    sys.ddata.type,
	    sys.box,
	    rcut1,
	    rcut2,
	    sum_nb_p.buff,
	    sum_nb_vxx.buff,
	    sum_nb_vyy.buff,
	    sum_nb_vzz.buff,
	    err.ptr_de);
  }
  else {
    size_t applyNonBondedInteraction_CellList_sbuffSize =
	(sizeof(IndexType) + sizeof(CoordType) + sizeof(TypeType)) *
	hroundUp4(clist.getCellBlockDim().x);
    calTwinRangeCorrection_cell
	<<<clist.getCellGrimDim(), clist.getCellBlockDim(),
	applyNonBondedInteraction_CellList_sbuffSize>>> (
	    sys.ddata.numAtom,
	    sys.ddata.coord,
	    twrec.forcx,
	    twrec.forcy,
	    twrec.forcz,
	    sys.ddata.type,
	    sys.box,
	    clist.dclist,
	    rcut1,
	    rcut2,
	    sum_nb_p.buff,
	    sum_nb_vxx.buff,
	    sum_nb_vyy.buff,
	    sum_nb_vzz.buff,
	    err.ptr_de);
  } 
  checkCUDAError ("TwinRangeCorrectionRecorder::calTwinRangeCorrection");
  err.check ("TwinRangeCorrectionRecorder::calTwinRangeCorrection");	
  cudaThreadSynchronize();
  MDStatistic st (sys);
  sum_nb_p.sumBuffAdd(st.ddata, mdStatisticNonBondedPotential, 0);
  sum_nb_vxx.sumBuffAdd(st.ddata, mdStatisticVirialXX);
  sum_nb_vyy.sumBuffAdd(st.ddata, mdStatisticVirialYY);
  sum_nb_vzz.sumBuffAdd(st.ddata, mdStatisticVirialZZ);
  cudaThreadSynchronize();
  st.updateHost ();
  twrec.energyCorrection() = st.nonBondedEnergy();
  twrec.pressureCorrection() = st.pressure(sys.box);
  if (timer != NULL) timer->toc(mdTimeNBInterTwinRange);
}


void InteractionEngine::
buildNeighborListCalTwinRangeCorrection (const MDSystem &		sys,
					 const CellList &		clist,
					 const ScalorType &		rcut1,
					 const ScalorType &		rcut2,
					 NeighborList &			nlist,
					 TwinRangeCorrectionRecorder &	twrec,
					 MDTimer *			timer)
{
  if (timer != NULL) timer->tic(mdTimeBuildNeighborList);
  if (clist.isempty()){
    size_t applyNonBondedInteraction_AllPair_sbuffSize =
	(sizeof(CoordType) + sizeof(TypeType)) *
	hroundUp4(myBlockDim.x);
    buildNeighborListCalTwinRangeCorr_all
	<<<atomGridDim, myBlockDim,
	applyNonBondedInteraction_AllPair_sbuffSize>>> (
	    sys.ddata.numAtom,
	    sys.ddata.coord,
	    twrec.forcx,
	    twrec.forcy,
	    twrec.forcz,
	    sys.ddata.type,
	    sys.box,
	    rcut1,
	    rcut2,
	    nlist.dnlist,
	    sum_nb_p.buff,
	    sum_nb_vxx.buff,
	    sum_nb_vyy.buff,
	    sum_nb_vzz.buff,
	    err.ptr_de);
  }
  else {
    size_t applyNonBondedInteraction_CellList_sbuffSize =
	(sizeof(IndexType) + sizeof(CoordType) + sizeof(TypeType)) *
	hroundUp4(clist.getCellBlockDim().x);
    buildNeighborListCalTwinRangeCorr_cell
	<<<clist.getCellGrimDim(), clist.getCellBlockDim(),
	applyNonBondedInteraction_CellList_sbuffSize>>> (
	    sys.ddata.numAtom,
	    sys.ddata.coord,
	    twrec.forcx,
	    twrec.forcy,
	    twrec.forcz,
	    sys.ddata.type,
	    sys.box,
	    clist.dclist,
	    rcut1,
	    rcut2,
	    nlist.dnlist,
	    sum_nb_p.buff,
	    sum_nb_vxx.buff,
	    sum_nb_vyy.buff,
	    sum_nb_vzz.buff,
	    err.ptr_de);
  } 
  checkCUDAError ("TwinRangeCorrectionRecorder::calTwinRangeCorrection");
  err.check ("TwinRangeCorrectionRecorder::calTwinRangeCorrection");	
  cudaThreadSynchronize();
  MDStatistic st (sys);
  sum_nb_p.sumBuffAdd(st.ddata, mdStatisticNonBondedPotential, 0);
  sum_nb_vxx.sumBuffAdd(st.ddata, mdStatisticVirialXX);
  sum_nb_vyy.sumBuffAdd(st.ddata, mdStatisticVirialYY);
  sum_nb_vzz.sumBuffAdd(st.ddata, mdStatisticVirialZZ);
  cudaThreadSynchronize();
  st.updateHost ();
  twrec.energyCorrection() = st.nonBondedEnergy();
  twrec.pressureCorrection() = st.pressure(sys.box);
  if (timer != NULL) timer->toc(mdTimeBuildNeighborList);
}




void InteractionEngine::
applyBondedInteraction (MDSystem & sys,
			const BondedInteractionList & bdlist,
			MDTimer *timer )
{
  if (hasBond) {
    if (timer != NULL) timer->tic(mdTimeBondedInteraction);
    calBondInteraction
	<<<atomGridDim, myBlockDim>>> (
	    sys.ddata.numAtom,
	    sys.ddata.coord,
	    sys.ddata.forcx,  sys.ddata.forcy,  sys.ddata.forcz,
	    sys.box,
	    bdlist.dbondlist);
    checkCUDAError ("InteractionEngine::applyInteraction bonded");
    err.check ("interaction engine b");	
    if (timer != NULL) timer->toc(mdTimeBondedInteraction);
  }
  if (hasAngle){
    if (timer != NULL) timer->tic(mdTimeAngleInteraction);
    calAngleInteraction
	<<<atomGridDim, myBlockDim>>> (
	    sys.ddata.numAtom,
	    sys.ddata.coord,
	    sys.ddata.forcx,  sys.ddata.forcy,  sys.ddata.forcz,
	    sys.box,
	    bdlist.danglelist);
    checkCUDAError ("InteractionEngine::applyInteraction angle");
    err.check ("interaction engine angle");	
    if (timer != NULL) timer->toc(mdTimeAngleInteraction);
  }
}

void InteractionEngine::
applyBondedInteraction (MDSystem & sys,
			const BondedInteractionList & bdlist,
			MDStatistic & st,
			MDTimer *timer)
{
  if (hasBond) {
    if (timer != NULL) timer->tic(mdTimeBInterStatistic);
    calBondInteraction
	<<<atomGridDim, myBlockDim,
	calBondInteraction_sbuffSize>>> (
	    sys.ddata.numAtom,
	    sys.ddata.coord,
	    sys.ddata.forcx,  sys.ddata.forcy,  sys.ddata.forcz,
	    sys.box,
	    bdlist.dbondlist
	    ,
	    sum_b_p.buff,
	    sum_b_vxx.buff,
	    sum_b_vyy.buff,
	    sum_b_vzz.buff,
	    err.ptr_de
	    );
    checkCUDAError ("InteractionEngine::applyInteraction bonded (with statistic)");
    err.check ("interaction engine");	
    if (timer != NULL) timer->toc(mdTimeBInterStatistic);
  }
  if (hasBond) {
    if (timer != NULL) timer->tic(mdTimeBInterStatistic);
    cudaThreadSynchronize();
    sum_b_p.sumBuffAdd(st.ddata, mdStatisticBondedPotential);
    sum_b_vxx.sumBuffAdd(st.ddata, mdStatisticVirialXX);
    sum_b_vyy.sumBuffAdd(st.ddata, mdStatisticVirialYY);
    sum_b_vzz.sumBuffAdd(st.ddata, mdStatisticVirialZZ);
    cudaThreadSynchronize();
    if (timer != NULL) timer->toc(mdTimeBInterStatistic);
    checkCUDAError ("InteractionEngine::applyInteraction sum bond statistic (with statistic)");
  }
  if (hasAngle){
    if (timer != NULL) timer->tic(mdTimeAngleInterStatistic);
    calAngleInteraction
	<<<atomGridDim, myBlockDim,
	calAngleInteraction_sbuffSize>>> (
	    sys.ddata.numAtom,
	    sys.ddata.coord,
	    sys.ddata.forcx,  sys.ddata.forcy,  sys.ddata.forcz,
	    sys.box,
	    bdlist.danglelist,
	    sum_angle_p.buff,
	    sum_angle_vxx.buff,
	    sum_angle_vyy.buff,
	    sum_angle_vzz.buff,
	    err.ptr_de);
    checkCUDAError ("InteractionEngine::applyInteraction angle");
    err.check ("interaction engine angle");	
    if (timer != NULL) timer->toc(mdTimeAngleInterStatistic);
  }
  if (hasAngle){
    if (timer != NULL) timer->tic(mdTimeAngleInterStatistic);
    sum_angle_p.sumBuffAdd(st.ddata, mdStatisticBondedPotential);
    sum_angle_vxx.sumBuffAdd(st.ddata, mdStatisticVirialXX);
    sum_angle_vyy.sumBuffAdd(st.ddata, mdStatisticVirialYY);
    sum_angle_vzz.sumBuffAdd(st.ddata, mdStatisticVirialZZ);
    cudaThreadSynchronize();
    if (timer != NULL) timer->toc(mdTimeAngleInterStatistic);
    checkCUDAError ("InteractionEngine::applyInteraction sum angle statistic (with statistic)");
  }
}


// void InteractionEngine::
// calculateWidomDeltaEnergy (const MDSystem & sys,
// 			   const NeighborList & nlist,
// 			   WidomTestParticleInsertion_NVT & wtest,
// 			   MDTimer * timer )
// {
//   if (timer != NULL) timer->tic(mdTimeNBInterStatistic);
//   // printf ("### %d\n", nlist.mode);
//   if (nlist.mode == CellListBuilt){
//     // printf ("### here %f\n", wtest.energyCorrection());
//     widomDeltaPoten_NVT
// 	<<<toGridDim(wtest.numTestParticle()),
// 	nlist.myBlockDim.x,
// 	nlist.myBlockDim.x * sizeof(ScalorType)>>> (
// 	    wtest.numTestParticle(),
// 	    wtest.coordTestParticle,
// 	    wtest.typeTestParticle,
// 	    sys.ddata.numAtom,
// 	    sys.ddata.coord,
// 	    sys.ddata.type,
// 	    sys.box,
// 	    nlist.dclist,
// 	    wtest.sumExpDeltaU.buff,
// 	    err.ptr_de);
//   }
//   else if (nlist.mode == AllPairBuilt){
//     // printf ("### here %f\n", wtest.energyCorrection());
//     widomDeltaPoten_allPair_NVT
// 	<<<toGridDim(wtest.numTestParticle()),
// 	DefaultNThreadPerBlock,
// 	DefaultNThreadPerBlock * sizeof(ScalorType)>>> (
// 	    wtest.numTestParticle(),
// 	    wtest.coordTestParticle,
// 	    wtest.typeTestParticle,
// 	    sys.ddata.numAtom,
// 	    sys.ddata.coord,
// 	    sys.ddata.type,
// 	    sys.box,
// 	    nlist.myrlist,
// 	    wtest.sumExpDeltaU.buff,
// 	    err.ptr_de);
//   }
//   if (timer != NULL) timer->toc(mdTimeNBInterStatistic);
// }

// void InteractionEngine::
// calculateWidomDeltaEnergy (const MDSystem & sys,
// 			   const NeighborList & nlist,
// 			   WidomTestParticleInsertion_NVT2 & wtest,
// 			   MDTimer * timer )
// {
//   if (timer != NULL) timer->tic(mdTimeNBInterStatistic);
//   // printf ("### %d\n", nlist.mode);
//   if (nlist.mode == CellListBuilt){
//     // printf ("### here %f\n", wtest.energyCorrection());
//     widomDeltaPoten_NVT
// 	<<<toGridDim(wtest.numTestParticle()),
// 	nlist.myBlockDim.x,
// 	nlist.myBlockDim.x * sizeof(ScalorType)>>> (
// 	    wtest.numTestParticle(),
// 	    wtest.coordTestParticle,
// 	    wtest.typeTestParticle,
// 	    sys.ddata.numAtom,
// 	    sys.ddata.coord,
// 	    sys.ddata.type,
// 	    sys.box,
// 	    nlist.dclist,
// 	    wtest.sumExpDeltaU.buff,
// 	    err.ptr_de);
//   }
//   if (timer != NULL) timer->toc(mdTimeNBInterStatistic);
// }

// void InteractionEngine::
// calculateWidomDeltaEnergy (const MDSystem & sys,
// 			   const NeighborList & nlist,
// 			   WidomTestParticleInsertion_NPT & wtest,
// 			   MDTimer * timer )
// {
//   if (timer != NULL) timer->tic(mdTimeNBInterStatistic);
//   // printf ("### %d\n", nlist.mode);
//   if (nlist.mode == CellListBuilt){
//     // printf ("### here %f, n: %d\n", wtest.energyCorrection(), wtest.numTestParticle());
//     widomDeltaPoten_NVT
// 	<<<toGridDim(wtest.numTestParticle()),
// 	nlist.myBlockDim.x,
// 	nlist.myBlockDim.x * sizeof(ScalorType)>>> (
// 	    wtest.numTestParticle(),
// 	    wtest.coordTestParticle,
// 	    wtest.typeTestParticle,
// 	    sys.ddata.numAtom,
// 	    sys.ddata.coord,
// 	    sys.ddata.type,
// 	    sys.box,
// 	    nlist.dclist,
// 	    wtest.sumExpDeltaU.buff,
// 	    err.ptr_de);
//   }
//   else if (nlist.mode == AllPairBuilt){
//     // printf ("### here %f\n", wtest.energyCorrection());
//     widomDeltaPoten_allPair_NVT
// 	<<<toGridDim(wtest.numTestParticle()),
// 	DefaultNThreadPerBlock,
// 	DefaultNThreadPerBlock * sizeof(ScalorType)>>> (
// 	    wtest.numTestParticle(),
// 	    wtest.coordTestParticle,
// 	    wtest.typeTestParticle,
// 	    sys.ddata.numAtom,
// 	    sys.ddata.coord,
// 	    sys.ddata.type,
// 	    sys.box,
// 	    nlist.myrlist,
// 	    wtest.sumExpDeltaU.buff,
// 	    err.ptr_de);
//   }
//   // for (unsigned i = 0; i < wtest.numTestParticle(); ++i){
//   //   printf ("%d %f  (%f %f %f)\n", i,
//   // 	    wtest.sumExpDeltaU.buff[i],
//   // 	    wtest.coordTestParticle[i].x,
//   // 	    wtest.coordTestParticle[i].y,
//   // 	    wtest.coordTestParticle[i].z
//   // 	);
//   // }
  
//   if (timer != NULL) timer->toc(mdTimeNBInterStatistic);
// }


__global__ void clearForce (const IndexType numAtom,
			    ScalorType * forcx,
			    ScalorType * forcy, 
			    ScalorType * forcz)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  if (ii < numAtom) {
    forcx[ii] = 0.0f;
    forcy[ii] = 0.0f;
    forcz[ii] = 0.0f;
  }
}


// __global__ void
// calNonBondedInteraction (const CoordType * coord,
// 			 const TypeType * type,
// 			 DeviceCellListData clist,
// 			 DeviceCellListProperty clistPro,
// 			 ScalorType * forcx,
// 			 ScalorType * forcy,
// 			 ScalorType * forcz,
// 			 bool sharednbForceTable)
// {
//   IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
//   IndexType tid = threadIdx.x;
//   ScalorType fsumx(0.f), fsumy(0.f), fsumz(0.f);
  
//   extern __shared__ volatile char pub_sbuff[];
//   volatile IndexType * targetIndex = (volatile IndexType *) pub_sbuff;
//   CoordType * targetCoord = (CoordType *) &targetIndex[roundUp4(blockDim.x)];
//   volatile TypeType * targetType = (volatile TypeType *) &targetCoord[roundUp4(blockDim.x)];
//   __syncthreads();

//   IndexType ii = get (clist, bid, tid);
//   CoordType ref;
//   TypeType refType;
//   if (ii != MaxIndexValue){
//     ref = tex1Dfetch (global_texRef_interaction_coord, ii);
//     refType = tex1Dfetch(global_texRef_interaction_type, ii);
//   }
  
//   for (unsigned i = 0; i < numNeighborCell(clistPro, bid); ++i){
//     __syncthreads();
//     IndexType targetCellIndex = getTargetCellIndex (clistPro, bid, i);
//     CoordType shift = getShiftValue (clistPro, bid, i);
//     IndexType targetIndex[tid] = get (clist, targetCellIndex, tid);
//     if (targetIndex[tid] != MaxIndexValue){
//       targetCoord[tid] = tex1Dfetch (global_texRef_interaction_coord, targetIndexes[tid]);
//       targetType[tid] = tex1Dfetch (global_texRef_interaction_type, targetIndexes[tid]);
//     }
//     __syncthreads ();
//     if (ii != MaxIndexValue){
//       for (IndexType jj = 0; jj < blockDim.x; ++jj){
// 	if (targetIndex[jj] == MaxIndexValue) continue;
// 	ScalorType diffx = targetCoord[jj].x + shift.x - ref.x;
// 	ScalorType diffy = targetCoord[jj].y + shift.y - ref.y;
// 	ScalorType diffz = targetCoord[jj].z + shift.z - ref.z;
// 	if ((diffx*diffx+diffy*diffy+diffz*diffz) < rlist2 &&
// 	    targetIndex[jj] != ii){
// 	  ForceIndexType fidx;
// 	  if (sharednbForceTable){
// 	    fidx = nonBondedInteractionTableItem (
// 		nonBondedInteractionTable, const_numAtomType, refType, targetType[jj]);
// 	  }
// 	  else {
// 	    fidx = nonBondedInteractionTableItem (
// 		nonBondedInteractionTable, const_numAtomType, refType, targetType[jj]);
// 	  }
// 	  ScalorType fx, fy, fz;
// 	  nbforce (nonBondedInteractionType[fidx],
// 		   &nonBondedInteractionParameter
// 		   [nonBondedInteractionParameterPosition[fidx]],
// 		   diffx, diffy, diffz,
// 		   &fx, &fy, &fz);
// 	  fsumx += fx;
// 	  fsumy += fy;
// 	  fsumz += fz;
// 	}
//       }
//     }
//   }

//   if (ii != MaxIndexValue){
//     forcx[ii] += fsumx;
//     forcy[ii] += fsumy;
//     forcz[ii] += fsumz;
//   }
// }


__global__ void
calNonBondedInteraction_neighbor (const IndexType		numAtom,
				  const CoordType *		coord,
				  ScalorType *			forcx,
				  ScalorType *			forcy, 
				  ScalorType *			forcz,
				  const TypeType *		type,
				  const RectangularBox		box,
				  const DeviceNeighborList	nlist)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  ScalorType fsumx = 0.0f;
  ScalorType fsumy = 0.0f;
  ScalorType fsumz = 0.0f;
  IndexType ii = tid + bid * blockDim.x;

  if (ii < numAtom) {
    CoordType ref (tex1Dfetch(global_texRef_interaction_coord, ii));
    ScalorType fx(0.f), fy(0.f), fz(0.f);
    for (IndexType jj = 0, nlistPosi = ii;
    	 jj < nlist.Nneighbor[ii];
    	 ++jj, nlistPosi += nlist.stride){
      IndexType targetIdx ( nlist.data [nlistPosi] );
      IndexType nbForceIndex ( nlist.forceIndex [nlistPosi] );
      CoordType target ( tex1Dfetch(global_texRef_interaction_coord, targetIdx) );
      ScalorType diffx ( target.x - ref.x );
      ScalorType diffy ( target.y - ref.y );
      ScalorType diffz ( target.z - ref.z );
      shortestImage (box, &diffx, &diffy, &diffz);
      nbForce (nonBondedInteractionType[nbForceIndex],
	       &nonBondedInteractionParameter
	       [nonBondedInteractionParameterPosition[nbForceIndex]],
      	       diffx, diffy, diffz, 
      	       &fx, &fy, &fz);
      fsumx += fx;
      fsumy += fy;
      fsumz += fz;
    }
    forcx[ii] += fsumx;
    forcy[ii] += fsumy;
    forcz[ii] += fsumz;
  }
}

__global__ void
calNonBondedInteraction_neighbor (const IndexType		numAtom,
				  const CoordType *		coord,
				  ScalorType *			forcx,
				  ScalorType *			forcy, 
				  ScalorType *			forcz,
				  const TypeType *		type,
				  const RectangularBox		box,
				  const DeviceNeighborList	nlist,
				  const DeviceExclusionList	dexcllist,
				  const bool			sharedExclusionList)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  ScalorType fsumx = 0.0f;
  ScalorType fsumy = 0.0f;
  ScalorType fsumz = 0.0f;
  IndexType ii = tid + bid * blockDim.x;

  IndexType * ptr_excllist;
  IndexType myNumExclusion (0);
  extern __shared__ char excl_sbuff[];
  
  if (dexcllist.maxNumExclusion != 0 && ii < numAtom){
    myNumExclusion  = dexcllist.numExclusion[ii];
    if (sharedExclusionList){
      ptr_excllist = (IndexType *) excl_sbuff;
      for (IndexType jj = 0; jj < myNumExclusion; ++jj){
	ptr_excllist[jj*blockDim.x+tid] =
	    dexcllist.exclusionNeighborIndex[jj*dexcllist.stride+ii];
      }
    }
  }

  if (ii < numAtom) {
    CoordType ref = tex1Dfetch(global_texRef_interaction_coord, ii);
    ScalorType fx(0.f), fy(0.f), fz(0.f);
    for (IndexType jj = 0, nlistPosi = ii;
    	 jj < nlist.Nneighbor[ii];
    	 ++jj, nlistPosi += nlist.stride){
      IndexType targetIdx ( nlist.data [nlistPosi] );
      IndexType nbForceIndex;
      CoordType target;
      ScalorType diffx, diffy, diffz;
      if (sharedExclusionList){
	for (IndexType kk = 0; kk < myNumExclusion; ++kk){
	  if (ptr_excllist[kk*blockDim.x+tid] == targetIdx) {
	    goto skipInter;
	  }
	}
      }
      else {
	for (IndexType kk = 0; kk < myNumExclusion; ++kk){
	  if (dexcllist.exclusionNeighborIndex[kk*dexcllist.stride+ii] == targetIdx) {
	    goto skipInter;
	  }	  
	}
      }
      nbForceIndex = ( nlist.forceIndex [nlistPosi] );
      target = ( tex1Dfetch(global_texRef_interaction_coord, targetIdx) );
      diffx = ( target.x - ref.x );
      diffy = ( target.y - ref.y );
      diffz = ( target.z - ref.z );
      shortestImage (box, &diffx, &diffy, &diffz);
      nbForce (nonBondedInteractionType[nbForceIndex],
	       &nonBondedInteractionParameter
	       [nonBondedInteractionParameterPosition[nbForceIndex]],
      	       diffx, diffy, diffz, 
      	       &fx, &fy, &fz);
      fsumx += fx;
      fsumy += fy;
      fsumz += fz;
      skipInter:
      {
      }
    }
    forcx[ii] += fsumx;
    forcy[ii] += fsumy;
    forcz[ii] += fsumz;
  }
}

__global__ void
calNonBondedInteraction_neighbor (const IndexType		numAtom,
				  const CoordType *		coord,
				  ScalorType *			forcx,
				  ScalorType *			forcy, 
				  ScalorType *			forcz,
				  const TypeType *		type,
				  const RectangularBox		box,
				  const DeviceNeighborList	nlist,
				  ScalorType *			statistic_nb_buff0,
				  ScalorType *			statistic_nb_buff1,
				  ScalorType *			statistic_nb_buff2,
				  ScalorType *			statistic_nb_buff3)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;

  ScalorType fsumx = 0.0f;
  ScalorType fsumy = 0.0f;
  ScalorType fsumz = 0.0f;
  IndexType ii = tid + bid * blockDim.x;
  ScalorType myPoten = 0.0f, myVxx = 0.0f, myVyy = 0.0f, myVzz = 0.0f;

  if (ii < numAtom) {
    CoordType ref;
    ref = tex1Dfetch(global_texRef_interaction_coord, ii);
    ScalorType fx(0.f), fy(0.f), fz(0.f);
    ScalorType dp;
    for (IndexType jj = 0, nlistPosi = ii;
	 jj < nlist.Nneighbor[ii];
	 ++jj, nlistPosi += nlist.stride){
      IndexType targetIdx ( nlist.data[nlistPosi] );
      IndexType nbForceIndex ( nlist.forceIndex [nlistPosi] );
      CoordType target ( tex1Dfetch(global_texRef_interaction_coord, targetIdx) );
      ScalorType diffx ( target.x - ref.x );
      ScalorType diffy ( target.y - ref.y );
      ScalorType diffz ( target.z - ref.z );
      shortestImage (box, &diffx, &diffy, &diffz);
      nbForcePoten (nonBondedInteractionType[nbForceIndex],
		    &nonBondedInteractionParameter
		    [nonBondedInteractionParameterPosition[nbForceIndex]],
      		    diffx, diffy, diffz, 
      		    &fx, &fy, &fz, &dp);
      // printf ("## %d\t%d\t%f\t%f\t%f\n",
      // 	      ii, targetIdx,
      // 	      ref.z, target.z, fz);
      // printf ("%f, %f %f %f,  %f %f %f,  %f %f %f, %f\n",
      // 	      sqrtf(diffx*diffx+diffy*diffy+diffz*diffz),
      // 	      ref.x, ref.y, ref.z,
      // 	      target.x, target.y, target.z,
      // 	      diffx, diffy, diffz,
      // 	      dp
      // 	  );
      myPoten += dp;
      myVxx += fx * diffx;
      myVyy += fy * diffy;
      myVzz += fz * diffz;
      fsumx += fx;
      fsumy += fy;
      fsumz += fz;
    }
    forcx[ii] += fsumx;
    forcy[ii] += fsumy;
    forcz[ii] += fsumz;
  }
  
  if (ii < numAtom){
    statistic_nb_buff0[ii] = myPoten * 0.5f;
    statistic_nb_buff1[ii] = myVxx * 0.5f;
    statistic_nb_buff2[ii] = myVyy * 0.5f;
    statistic_nb_buff3[ii] = myVzz * 0.5f;
  }  
}


__global__ void
calNonBondedInteraction_neighbor (const IndexType		numAtom,
				  const CoordType *		coord,
				  ScalorType *			forcx,
				  ScalorType *			forcy, 
				  ScalorType *			forcz,
				  const TypeType *		type,
				  const RectangularBox		box,
				  const DeviceNeighborList	nlist,
				  const DeviceExclusionList	dexcllist,
				  const bool			sharedExclusionList,
				  ScalorType *			statistic_nb_buff0,
				  ScalorType *			statistic_nb_buff1,
				  ScalorType *			statistic_nb_buff2,
				  ScalorType *			statistic_nb_buff3)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;

  ScalorType fsumx = 0.0f;
  ScalorType fsumy = 0.0f;
  ScalorType fsumz = 0.0f;
  IndexType ii = tid + bid * blockDim.x;
  ScalorType myPoten = 0.0f, myVxx = 0.0f, myVyy = 0.0f, myVzz = 0.0f;

  IndexType * ptr_excllist;
  IndexType myNumExclusion (0);
  extern __shared__ char excl_sbuff[];
  
  if (dexcllist.maxNumExclusion != 0 && ii < numAtom){
    myNumExclusion  = dexcllist.numExclusion[ii];
    if (sharedExclusionList){
      ptr_excllist = (IndexType *) excl_sbuff;
      for (IndexType jj = 0; jj < myNumExclusion; ++jj){
	ptr_excllist[jj*blockDim.x+tid] =
	    dexcllist.exclusionNeighborIndex[jj*dexcllist.stride+ii];
      }
    }
  }

  if (ii < numAtom) {
    CoordType ref;
    ref = tex1Dfetch(global_texRef_interaction_coord, ii);
    ScalorType fx(0.f), fy(0.f), fz(0.f);
    ScalorType dp;
    for (IndexType jj = 0, nlistPosi = ii;
	 jj < nlist.Nneighbor[ii];
	 ++jj, nlistPosi += nlist.stride){
      IndexType targetIdx ( nlist.data[nlistPosi] );
      IndexType nbForceIndex;
      CoordType target;
      ScalorType diffx, diffy, diffz;
      if (sharedExclusionList){
	for (IndexType kk = 0; kk < myNumExclusion; ++kk){
	  if (ptr_excllist[kk*blockDim.x+tid] == targetIdx) {
	    goto skipInter;
	  }
	}
      }
      else {
	for (IndexType kk = 0; kk < myNumExclusion; ++kk){
	  if (dexcllist.exclusionNeighborIndex[kk*dexcllist.stride+ii] == targetIdx) {
	    goto skipInter;
	  }	  
	}
      }
      nbForceIndex = ( nlist.forceIndex [nlistPosi] );
      target = ( tex1Dfetch(global_texRef_interaction_coord, targetIdx) );
      diffx = ( target.x - ref.x );
      diffy = ( target.y - ref.y );
      diffz = ( target.z - ref.z );
      shortestImage (box, &diffx, &diffy, &diffz);
      nbForcePoten (nonBondedInteractionType[nbForceIndex],
		    &nonBondedInteractionParameter
		    [nonBondedInteractionParameterPosition[nbForceIndex]],
      		    diffx, diffy, diffz, 
      		    &fx, &fy, &fz, &dp);
      // printf ("## %d\t%d\t%f\t%f\t%f\n",
      // 	      ii, targetIdx,
      // 	      ref.z, target.z, fz);
      // printf ("%f, %f %f %f,  %f %f %f,  %f %f %f, %f\n",
      // 	      sqrtf(diffx*diffx+diffy*diffy+diffz*diffz),
      // 	      ref.x, ref.y, ref.z,
      // 	      target.x, target.y, target.z,
      // 	      diffx, diffy, diffz,
      // 	      dp
      // 	  );
      myPoten += dp;
      myVxx += fx * diffx;
      myVyy += fy * diffy;
      myVzz += fz * diffz;
      fsumx += fx;
      fsumy += fy;
      fsumz += fz;
      skipInter:
      {
      }
    }
    forcx[ii] += fsumx;
    forcy[ii] += fsumy;
    forcz[ii] += fsumz;
  }
  
  if (ii < numAtom){
    statistic_nb_buff0[ii] = myPoten * 0.5f;
    statistic_nb_buff1[ii] = myVxx * 0.5f;
    statistic_nb_buff2[ii] = myVyy * 0.5f;
    statistic_nb_buff3[ii] = myVzz * 0.5f;
  }  
}





__global__ void calBondInteraction (const IndexType numAtom,
				    const CoordType * coord,
				    ScalorType * forcx,
				    ScalorType * forcy, 
				    ScalorType * forcz,
				    const RectangularBox box,
				    const DeviceBondList bdlist)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  ScalorType fsumx = 0.0f;
  ScalorType fsumy = 0.0f;
  ScalorType fsumz = 0.0f;
  IndexType ii = tid + bid * blockDim.x;
  
  if (ii >= numAtom) return;
  CoordType ref;
#ifdef COMPILE_NO_TEX
  ref = coord[ii];
#else
  ref = tex1Dfetch(global_texRef_interaction_coord, ii);
#endif  
      
  IndexType myNumBond = bdlist.numBond[ii];
  
  for (IndexType jj = 0; jj < bdlist.maxNumBond; ++jj){
    if (jj == myNumBond) break;
    IndexType targetIdx = bdlist.bondNeighborIndex[jj * bdlist.stride + ii];
    CoordType target;
#ifdef COMPILE_NO_TEX
    target = coord[targetIdx];
#else
    target = tex1Dfetch(global_texRef_interaction_coord, targetIdx);
#endif 
    ScalorType diffx, diffy, diffz;
    diffx = target.x - ref.x;
    diffy = target.y - ref.y;
    diffz = target.z - ref.z;
    shortestImage (box, &diffx, &diffy, &diffz);
    ScalorType fx, fy, fz;
    IndexType bondFindex = bdlist.bondIndex[jj * bdlist.stride + ii];
    bondForce (bondedInteractionType[bondFindex],
	       &bondedInteractionParameter
	       [bondedInteractionParameterPosition[bondFindex]],
	       diffx, diffy, diffz, &fx, &fy, &fz);
    fsumx += fx;
    fsumy += fy;
    fsumz += fz;
  }
  forcx[ii] += fsumx;
  forcy[ii] += fsumy;
  forcz[ii] += fsumz;
}


__global__ void calBondInteraction (const IndexType numAtom,
				    const CoordType * coord,
				    ScalorType * forcx,
				    ScalorType * forcy, 
				    ScalorType * forcz,
				    const RectangularBox box,
				    const DeviceBondList bdlist,
				    ScalorType * statistic_b_buff0,
				    ScalorType * statistic_b_buff1,
				    ScalorType * statistic_b_buff2,
				    ScalorType * statistic_b_buff3,
				    mdError_t * ptr_de)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;

  extern __shared__ volatile ScalorType buff[];
  buff[tid] = 0.f;
  __syncthreads();
  
  ScalorType fsumx = 0.0f;
  ScalorType fsumy = 0.0f;
  ScalorType fsumz = 0.0f;
  IndexType ii = tid + bid * blockDim.x;
  ScalorType myPoten = 0.0f, myVxx = 0.0f, myVyy = 0.0f, myVzz = 0.0f;
  if (ii < numAtom) {
    CoordType ref;
#ifdef COMPILE_NO_TEX
    ref = coord[ii];
#else 
    ref = tex1Dfetch(global_texRef_interaction_coord, ii);
#endif
    IndexType myNumBond = bdlist.numBond[ii];
    for (IndexType jj = 0; jj < bdlist.maxNumBond; ++jj){
      if (jj == myNumBond) break;
      IndexType targetIdx = bdlist.bondNeighborIndex[jj * bdlist.stride + ii];
      CoordType target;
#ifdef COMPILE_NO_TEX
      target = coord[targetIdx];
#else
      target = tex1Dfetch(global_texRef_interaction_coord, targetIdx);
#endif
      ScalorType diffx, diffy, diffz;
      diffx = target.x - ref.x;
      diffy = target.y - ref.y;
      diffz = target.z - ref.z;
      shortestImage (box, &diffx, &diffy, &diffz);
      ScalorType fx, fy, fz;
      IndexType bondFindex = bdlist.bondIndex[jj * bdlist.stride + ii];
      ScalorType dp;
      bondForcePoten (bondedInteractionType[bondFindex],
		      &bondedInteractionParameter
		      [bondedInteractionParameterPosition[bondFindex]],
		      diffx, diffy, diffz, &fx, &fy, &fz, &dp);
      myPoten += dp;
      myVxx += fx * diffx;
      myVyy += fy * diffy;
      myVzz += fz * diffz;
      fsumx += fx;
      fsumy += fy;
      fsumz += fz;
    }
    forcx[ii] += fsumx;
    forcy[ii] += fsumy;
    forcz[ii] += fsumz;
  }

  buff[tid] = myPoten * 0.5f;
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) statistic_b_buff0[bid] = buff[0];
  __syncthreads();
  buff[tid] = myVxx * 0.5f;
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) statistic_b_buff1[bid] = buff[0];
  __syncthreads();
  buff[tid] = myVyy * 0.5f;
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) statistic_b_buff2[bid] = buff[0];
  __syncthreads();
  buff[tid] = myVzz * 0.5f;
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) statistic_b_buff3[bid] = buff[0];
  __syncthreads();
}



__global__ void calAngleInteraction (const IndexType numAtom,
				     const CoordType * coord,
				     ScalorType * forcx,
				     ScalorType * forcy, 
				     ScalorType * forcz,
				     const RectangularBox box,
				     const DeviceAngleList anglelist)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  ScalorType fsumx = 0.0f;
  ScalorType fsumy = 0.0f;
  ScalorType fsumz = 0.0f;
  IndexType ii = tid + bid * blockDim.x;
  IndexType myNumAngle;
  
  if (ii < numAtom){
    myNumAngle = anglelist.numAngle[ii];  
  }
  else {
    myNumAngle = 0;
    return ;
  }
  // if (__all(myNumAngle == 0)) return ;
  
  if (ii < numAtom) {
    CoordType ref;
#ifdef COMPILE_NO_TEX
    ref = coord[ii];
#else
    ref = tex1Dfetch(global_texRef_interaction_coord, ii);
#endif
    for (IndexType jj = 0; jj < myNumAngle; ++jj){
      IndexType targetIdx0 =
	  anglelist.angleNeighborIndex[((jj<<1)  ) * anglelist.stride + ii];
      IndexType targetIdx1 =
	  anglelist.angleNeighborIndex[((jj<<1)+1) * anglelist.stride + ii];
      IndexType myPosi     =
	  anglelist.anglePosi[jj * anglelist.stride + ii];
      CoordType target0, target1;
#ifdef COMPILE_NO_TEX
      target0 = coord[targetIdx0];
      target1 = coord[targetIdx1];
#else
      target0 = tex1Dfetch(global_texRef_interaction_coord, targetIdx0);
      target1 = tex1Dfetch(global_texRef_interaction_coord, targetIdx1);
#endif 
      ScalorType diff0x, diff0y, diff0z;
      ScalorType diff1x, diff1y, diff1z;
      bool center (myPosi == 1);
      if (center){
	diff0x = ref.x - target0.x;
	diff0y = ref.y - target0.y;
	diff0z = ref.z - target0.z;
	diff1x = target1.x -  ref.x;
	diff1y = target1.y -  ref.y;
	diff1z = target1.z -  ref.z;
      } else {
	diff0x = target0.x - ref.x;
	diff0y = target0.y - ref.y;
	diff0z = target0.z - ref.z;
	diff1x = target1.x - target0.x;
	diff1y = target1.y - target0.y;
	diff1z = target1.z - target0.z;
      }      
      shortestImage (box, &diff0x, &diff0y, &diff0z);
      shortestImage (box, &diff1x, &diff1y, &diff1z);
      ScalorType f0x, f0y, f0z;
      ScalorType f1x, f1y, f1z;
      IndexType angleFindex = anglelist.angleIndex[jj * anglelist.stride + ii];
      angleForce (center,
		  bondedInteractionType[angleFindex],
		  &bondedInteractionParameter
		  [bondedInteractionParameterPosition[angleFindex]],
		  diff0x, diff0y, diff0z,
		  diff1x, diff1y, diff1z,
		  &f0x, &f0y, &f0z,
		  &f1x, &f1y, &f1z);
      if (center){
	fsumx += f0x + f1x;
	fsumy += f0y + f1y;
	fsumz += f0z + f1z;
      }
      else {
	fsumx -= f0x;
	fsumy -= f0y;
	fsumz -= f0z;
      }
    }
    forcx[ii] += fsumx;
    forcy[ii] += fsumy;
    forcz[ii] += fsumz;
  }
}


__global__ void calAngleInteraction (const IndexType numAtom,
				     const CoordType * coord,
				     ScalorType * forcx,
				     ScalorType * forcy, 
				     ScalorType * forcz,
				     const RectangularBox box,
				     const DeviceAngleList anglelist,
				     ScalorType * statistic_b_buff0,
				     ScalorType * statistic_b_buff1,
				     ScalorType * statistic_b_buff2,
				     ScalorType * statistic_b_buff3,
				     mdError_t * ptr_de)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  ScalorType fsumx = 0.0f;
  ScalorType fsumy = 0.0f;
  ScalorType fsumz = 0.0f;
  ScalorType myPoten = 0.0f, myVxx = 0.0f, myVyy = 0.0f, myVzz = 0.0f;
  IndexType ii = tid + bid * blockDim.x;
  IndexType myNumAngle;
  extern __shared__ volatile ScalorType buff[];
  buff[tid] = 0.f;
  __syncthreads();
  
  if (ii < numAtom) {
    CoordType ref;
#ifdef COMPILE_NO_TEX
    ref = coord[ii];
#else
    ref = tex1Dfetch(global_texRef_interaction_coord, ii);
#endif
    myNumAngle = anglelist.numAngle[ii];  
    for (IndexType jj = 0; jj < myNumAngle; ++jj){
      IndexType targetIdx0 =
	  anglelist.angleNeighborIndex[((jj<<1)  ) * anglelist.stride + ii];
      IndexType targetIdx1 =
	  anglelist.angleNeighborIndex[((jj<<1)+1) * anglelist.stride + ii];
      IndexType myPosi     =
	  anglelist.anglePosi[jj * anglelist.stride + ii];
      CoordType target0, target1;
#ifdef COMPILE_NO_TEX
      target0 = coord[targetIdx0];
      target1 = coord[targetIdx1];
#else
      target0 = tex1Dfetch(global_texRef_interaction_coord, targetIdx0);
      target1 = tex1Dfetch(global_texRef_interaction_coord, targetIdx1);
#endif 
      ScalorType diff0x, diff0y, diff0z;
      ScalorType diff1x, diff1y, diff1z;
      bool center = (myPosi == 1);
      if (center){
	diff0x = ref.x - target0.x;
	diff0y = ref.y - target0.y;
	diff0z = ref.z - target0.z;
	diff1x = target1.x -  ref.x;
	diff1y = target1.y -  ref.y;
	diff1z = target1.z -  ref.z;
      } else {
	diff0x = target0.x - ref.x;
	diff0y = target0.y - ref.y;
	diff0z = target0.z - ref.z;
	diff1x = target1.x - target0.x;
	diff1y = target1.y - target0.y;
	diff1z = target1.z - target0.z;
      }      
      shortestImage (box, &diff0x, &diff0y, &diff0z);
      shortestImage (box, &diff1x, &diff1y, &diff1z);
      ScalorType f0x, f0y, f0z;
      ScalorType f1x, f1y, f1z;
      IndexType angleFindex = anglelist.angleIndex[jj * anglelist.stride + ii];
      ScalorType dp;
      angleForcePoten (center,
		       bondedInteractionType[angleFindex],
		       &bondedInteractionParameter
		       [bondedInteractionParameterPosition[angleFindex]],
		       diff0x, diff0y, diff0z,
		       diff1x, diff1y, diff1z,
		       &f0x, &f0y, &f0z,
		       &f1x, &f1y, &f1z,
		       &dp);
      myPoten += dp;
      if (center){
	fsumx += f0x + f1x;
	fsumy += f0y + f1y;
	fsumz += f0z + f1z;
	myVxx -= f0x * diff0x - f1x * diff1x;
	myVyy -= f0y * diff0y - f1y * diff1y;
	myVzz -= f0z * diff0z - f1z * diff1z;
      }
      else {
	fsumx -= f0x;
	fsumy -= f0y;
	fsumz -= f0z;
      }
    }
    forcx[ii] += fsumx;
    forcy[ii] += fsumy;
    forcz[ii] += fsumz;
  }

  buff[tid] = myPoten * 0.33333333333333333f;
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) statistic_b_buff0[bid] = buff[0];
  __syncthreads();
  buff[tid] = myVxx;
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) statistic_b_buff1[bid] = buff[0];
  __syncthreads();
  buff[tid] = myVyy;
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) statistic_b_buff2[bid] = buff[0];
  __syncthreads();
  buff[tid] = myVzz;
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) statistic_b_buff3[bid] = buff[0];
  __syncthreads();
}



// static __device__ IndexType shiftedD3toD1 (
//     DeviceCellList clist,
//     RectangularBoxGeometry::RectangularBox box,
//     int ix, int iy, int iz,
//     ScalorType * shiftx , ScalorType * shifty, ScalorType * shiftz)
// {
//   int tmp;
//   ix += (tmp = -int(floorf(ix * clist.NCelli.x))) * clist.NCell.x;
//   *shiftx = tmp * box.size.x;
//   iy += (tmp = -int(floorf(iy * clist.NCelli.y))) * clist.NCell.y;
//   *shifty = tmp * box.size.y;
//   iz += (tmp = -int(floorf(iz * clist.NCelli.z))) * clist.NCell.z;
//   *shiftz = tmp * box.size.z;
//   return D3toD1 (clist.NCell, ix, iy, iz);
// }



// __global__ void calNonBondedInteraction (
//     const IndexType numAtom,
//     const CoordType * coord,
//     ScalorType * forcx,
//     ScalorType * forcy, 
//     ScalorType * forcz,
//     const TypeType * type,
//     const RectangularBox box,
//     DeviceCellList clist,
//     mdError_t * ptr_de)
// {
//   // RectangularBoxGeometry::normalizeSystem (box, &ddata);
//   IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
//   IndexType tid = threadIdx.x;
//   IndexType bidx, bidy, bidz;
//   D1toD3 (clist.NCell, bid, bidx, bidy, bidz);
  
//   // load index
//   IndexType ii = getDeviceCellListData (clist, bid, tid);
//   // load iith coordinate // use texturefetch instead
//   CoordType ref;
//   TypeType reftype;
//   ScalorType fsumx (0.f), fsumy(0.f), fsumz(0.f);
//   if (ii != MaxIndexValue){
// #ifdef COMPILE_NO_TEX
//     ref = coord[ii];
//     reftype = type[ii];
// #else
//     ref = tex1Dfetch (global_texRef_interaction_coord, ii);
//     // reftype = tex1Dfetch(global_texRef_interaction_type, ii);
// #endif
//   }
//   ScalorType rlist = clist.rlist;

//   // the target index and coordinates are shared

//   extern __shared__ volatile char pub_sbuff[];
  
//   volatile IndexType * targetIndexes =
//       (volatile IndexType *) pub_sbuff;
//   CoordType * target =
//       (CoordType *) &targetIndexes[roundUp4(blockDim.x)];
//   volatile TypeType * targettype =
//       (volatile TypeType *) &target[roundUp4(blockDim.x)];
//   __syncthreads();

//   // bool oneCellX(false), oneCellY(false), oneCellZ(false);
//   // if (clist.NCell.x == 1) oneCellX = true;
//   // if (clist.NCell.y == 1) oneCellY = true;
//   // if (clist.NCell.z == 1) oneCellZ = true;
//   // int upperx(1), lowerx(-1);
//   // int uppery(1), lowery(-1);
//   // int upperz(1), lowerz(-1);
//   // if (oneCellX) {lowerx =  0; upperx = 0;}
//   // if (oneCellY) {lowery =  0; uppery = 0;}
//   // if (oneCellZ) {lowerz =  0; upperz = 0;}
//   ScalorType rlist2 = rlist * rlist;
  
//   // loop over 27 neighbor cells
// #pragma unroll 3
//   // for (int nci = bidx + lowerx; nci <= bidx + upperx; ++nci){
//   //   for (int ncj = bidy + lowery; ncj <= bidy + uppery; ++ncj){
//   //     for (int nck = bidz + lowerz; nck <= bidz + upperz; ++nck){
//   for (int nci = int(bidx) - 1; nci <= int(bidx) + 1; ++nci){
//     for (int ncj = int(bidy) - 1; ncj <= int(bidy) + 1; ++ncj){
//       for (int nck = int(bidz) - 1; nck <= int(bidz) + 1; ++nck){
//   // for (int di = lowerx; di <= upperx; ++di){
//   //   for (int dj = lowery; dj <= uppery; ++dj){
//   //     for (int dk = lowerz; dk <= upperz; ++dk){
// 	__syncthreads();
// 	// the shift value of a cell is pre-computed
// 	ScalorType xshift, yshift, zshift;
// 	// int nci = di + bidx;
// 	// int ncj = dj + bidy;
// 	// int nck = dk + bidz;
// 	IndexType targetCellIdx = shiftedD3toD1 (clist, box, 
// 						 nci, ncj, nck, 
// 						 &xshift, &yshift, &zshift);
// 	// load target index and coordinates form global memary
// 	// IndexType tmp = (targetIndexes[tid] = 
// 	// 		 getDeviceCellListData(clist, targetCellIdx, tid));
// 	targetIndexes[tid] = getDeviceCellListData(clist, targetCellIdx, tid);
// 	if (targetIndexes[tid] != MaxIndexValue){
// // #ifdef COMPILE_NO_TEX
// // 	  target[tid] = coord[tmp];
// // 	  // targettype[tid] = type[tmp];
// // #else
// 	  target[tid] = tex1Dfetch(global_texRef_interaction_coord, targetIndexes[tid]);
// 	  // targettype[tid] = tex1Dfetch(global_texRef_interaction_type, tmp);
// // #endif
// 	}
// 	__syncthreads();
// 	// find neighbor
// 	if (ii != MaxIndexValue){
// 	  for (IndexType jj = 0; jj < clist.numbers[targetCellIdx]; ++jj){
// 	    // if (targetIndexes[jj] == MaxIndexValue) break;
// 	    ScalorType diffx = target[jj].x - xshift - ref.x;
// 	    ScalorType diffy = target[jj].y - yshift - ref.y;
// 	    ScalorType diffz = target[jj].z - zshift - ref.z;
// 	    // if (oneCellX) shortestImage (box.size.x, box.sizei.x, &diffx);
// 	    // if (oneCellY) shortestImage (box.size.y, box.sizei.y, &diffy);
// 	    // if (oneCellZ) shortestImage (box.size.z, box.sizei.z, &diffz);
// 	    //printf ("%d\t%d\t%f\t%f\n", ii,
// 	    ScalorType dr2;
// 	    if ((dr2 = (diffx*diffx+diffy*diffy+diffz*diffz)) < rlist2 &&
// 		targetIndexes[jj] != ii){
// 	      IndexType fidx(0);
// 	      // fidx = AtomNBForceTable::calForceIndex (
// 	      // 	  nonBondedInteractionTable,
// 	      // 	  const_numAtomType[0],
// 	      // 	  reftype,
// 	      // 	  targettype[jj]);
// 	      // if (fidx != mdForceNULL) {
// 	      ScalorType fx, fy, fz;
// 	      nbForce (nonBondedInteractionType[fidx],
// 		       &nonBondedInteractionParameter
// 		       [nonBondedInteractionParameterPosition[fidx]],
// 		       diffx, diffy, diffz,
// 		       dr2,
// 		       &fx, &fy, &fz);
// 	      fsumx += fx;
// 	      fsumy += fy;
// 	      fsumz += fz;
// 	      // }
// 	    }
// 	  }
// 	}
//       }
//     }
//   }
//   if (ii != MaxIndexValue){
//     forcx[ii] += fsumx;
//     forcy[ii] += fsumy;
//     forcz[ii] += fsumz;
//   }
// }







// __global__ void calNonBondedInteraction (
//     const IndexType numAtom,
//     const CoordType * coord,
//     ScalorType * forcx,
//     ScalorType * forcy, 
//     ScalorType * forcz,
//     const TypeType * type,
//     const RectangularBox box,
//     DeviceCellList clist,
//     ScalorType * statistic_nb_buff0,
//     ScalorType * statistic_nb_buff1,
//     ScalorType * statistic_nb_buff2,
//     ScalorType * statistic_nb_buff3,
//     mdError_t * ptr_de)
// {
//   // RectangularBoxGeometry::normalizeSystem (box, &ddata);
//   IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
//   IndexType tid = threadIdx.x;
//   IndexType bidx, bidy, bidz;
//   D1toD3 (clist.NCell, bid, bidx, bidy, bidz);
  
//   // load index
//   IndexType ii = getDeviceCellListData (clist, bid, tid);
//   // load iith coordinate // use texturefetch instead
//   CoordType ref;
//   TypeType reftype;
//   ScalorType fsumx (0.f), fsumy(0.f), fsumz(0.f);
//   ScalorType myPoten (0.0f), myVxx (0.0f), myVyy (0.0f), myVzz (0.0f);
//   if (ii != MaxIndexValue){
// #ifdef COMPILE_NO_TEX
//     ref = coord[ii];
//     reftype = type[ii];
// #else
//     ref = tex1Dfetch (global_texRef_interaction_coord, ii);
//     reftype = tex1Dfetch(global_texRef_interaction_type, ii);
// #endif
//   }
//   ScalorType rlist = clist.rlist;

//   // the target index and coordinates are shared

//   extern __shared__ volatile char pub_sbuff[];
  
//   volatile IndexType * targetIndexes =
//       (volatile IndexType *) pub_sbuff;
//   CoordType * target =
//       (CoordType *) &targetIndexes[roundUp4(blockDim.x)];
//   volatile TypeType * targettype =
//       (volatile TypeType *) &target[roundUp4(blockDim.x)];
//   __syncthreads();

//   // bool oneCellX(false), oneCellY(false), oneCellZ(false);
//   // if (clist.NCell.x == 1) oneCellX = true;
//   // if (clist.NCell.y == 1) oneCellY = true;
//   // if (clist.NCell.z == 1) oneCellZ = true;
//   // int upperx(1), lowerx(-1);
//   // int uppery(1), lowery(-1);
//   // int upperz(1), lowerz(-1);
//   // if (oneCellX) {lowerx =  0; upperx = 0;}
//   // if (oneCellY) {lowery =  0; uppery = 0;}
//   // if (oneCellZ) {lowerz =  0; upperz = 0;}
//   ScalorType rlist2 = rlist * rlist;
  
//   // loop over 27 neighbor cells
// #pragma unroll 3
//   for (int nci = int(bidx) - 1; nci <= int(bidx) + 1; ++nci){
//     for (int ncj = int(bidy) - 1; ncj <= int(bidy) + 1; ++ncj){
//       for (int nck = int(bidz) - 1; nck <= int(bidz) + 1; ++nck){
//   // for (int di = lowerx; di <= upperx; ++di){
//   //   for (int dj = lowery; dj <= uppery; ++dj){
//   //     for (int dk = lowerz; dk <= upperz; ++dk){
// 	__syncthreads();
// 	// the shift value of a cell is pre-computed
// 	ScalorType xshift, yshift, zshift;
// 	// int nci = di + bidx;
// 	// int ncj = dj + bidy;
// 	// int nck = dk + bidz;
// 	IndexType targetCellIdx = shiftedD3toD1 (clist, box, 
// 						 nci, ncj, nck, 
// 						 &xshift, &yshift, &zshift);
// 	// load target index and coordinates form global memary
// 	IndexType tmp = (targetIndexes[tid] = 
// 			 getDeviceCellListData(clist, targetCellIdx, tid));
// 	if (tmp != MaxIndexValue){
// #ifdef COMPILE_NO_TEX
// 	  target[tid] = coord[tmp];
// 	  targettype[tid] = type[tmp];
// #else
// 	  target[tid] = tex1Dfetch(global_texRef_interaction_coord, tmp);
// 	  targettype[tid] = tex1Dfetch(global_texRef_interaction_type, tmp);
// #endif
// 	}
// 	__syncthreads();
// 	// find neighbor
// 	if (ii != MaxIndexValue){
// 	  for (IndexType jj = 0; jj < clist.numbers[targetCellIdx]; ++jj){
// 	    ScalorType diffx = target[jj].x - xshift - ref.x;
// 	    ScalorType diffy = target[jj].y - yshift - ref.y;
// 	    ScalorType diffz = target[jj].z - zshift - ref.z;
// 	    // if (oneCellX) shortestImage (box.size.x, box.sizei.x, &diffx);
// 	    // if (oneCellY) shortestImage (box.size.y, box.sizei.y, &diffy);
// 	    // if (oneCellZ) shortestImage (box.size.z, box.sizei.z, &diffz);
// 	    //printf ("%d\t%d\t%f\t%f\n", ii, 
// 	    if ((diffx*diffx+diffy*diffy+diffz*diffz) < rlist2 &&
// 		targetIndexes[jj] != ii){
// 	      IndexType fidx(0);
// 	      // fidx = AtomNBForceTable::calForceIndex (
// 	      // 	  nonBondedInteractionTable,
// 	      // 	  const_numAtomType[0],
// 	      // 	  reftype,
// 	      // 	  targettype[jj]);
// 	      // if (fidx != mdForceNULL) {
// 	      ScalorType fx, fy, fz, dp;
// 	      nbForcePoten (nonBondedInteractionType[fidx],
// 			    &nonBondedInteractionParameter
// 			    [nonBondedInteractionParameterPosition[fidx]],
// 			    diffx, diffy, diffz, 
// 			    &fx, &fy, &fz, &dp);
// 	      myPoten += dp;
// 	      myVxx += fx * diffx;
// 	      myVyy += fy * diffy;
// 	      myVzz += fz * diffz;
// 	      fsumx += fx;
// 	      fsumy += fy;
// 	      fsumz += fz;
// 	      // }
// 	    }
// 	  }
// 	}
//       }
//     }
//   }
//   if (ii != MaxIndexValue){
//     forcx[ii] += fsumx;
//     forcy[ii] += fsumy;
//     forcz[ii] += fsumz;
//     statistic_nb_buff0[ii] = myPoten * 0.5f;
//     statistic_nb_buff1[ii] = myVxx * 0.5f;
//     statistic_nb_buff2[ii] = myVyy * 0.5f;
//     statistic_nb_buff3[ii] = myVzz * 0.5f;
//   }  
// }



__global__ void
calNonBondedInteraction_cell (const IndexType		numAtom,
			      const CoordType *		coord,
			      ScalorType *		forcx,
			      ScalorType *		forcy, 
			      ScalorType *		forcz,
			      const TypeType *		type,
			      const RectangularBox	box,
			      const DeviceCellList	clist,
			      const ScalorType		rcut,
			      mdError_t *		ptr_de)
{
  // RectangularBoxGeometry::normalizeSystem (box, &ddata);
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType bidx, bidy, bidz;
  D1toD3 (clist.NCell, bid, bidx, bidy, bidz);
  
  // load index
  IndexType ii = getDeviceCellListData (clist, bid, tid);
  // load iith coordinate // use texturefetch instead
  CoordType ref;
  TypeType reftype;
  ScalorType fsumx (0.f), fsumy(0.f), fsumz(0.f);
  if (ii != MaxIndexValue){
#ifdef COMPILE_NO_TEX
    ref = coord[ii];
    reftype = type[ii];
#else
    ref = tex1Dfetch (global_texRef_interaction_coord, ii);
    reftype = tex1Dfetch(global_texRef_interaction_type, ii);
#endif
  }

  // the target index and coordinates are shared

  extern __shared__ volatile char pub_sbuff[];
  
  volatile IndexType * targetIndexes =
      (volatile IndexType *) pub_sbuff;
  CoordType * target =
      (CoordType *) &targetIndexes[roundUp4(blockDim.x)];
  volatile TypeType * targettype =
      (volatile TypeType *) &target[roundUp4(blockDim.x)];
  
  ScalorType rcut2 = rcut * rcut;
  for (IndexType i = 0; i < clist.numNeighborCell[bid]; ++i){
    __syncthreads();
    IndexType targetCellIdx = getNeighborCellIndex (clist, bid, i);
    CoordNoiType shiftNoi   = getNeighborCellShiftNoi (clist, bid, i);
    CoordType shift;
    shift.x = shiftNoi.x * box.size.x;
    shift.y = shiftNoi.y * box.size.y;
    shift.z = shiftNoi.z * box.size.z;    
    targetIndexes[tid] = getDeviceCellListData(clist, targetCellIdx, tid);  
    if (targetIndexes[tid] != MaxIndexValue){
      target[tid] = tex1Dfetch(global_texRef_interaction_coord, targetIndexes[tid]);
      targettype[tid] = tex1Dfetch(global_texRef_interaction_type, targetIndexes[tid]);
    }
    bool oneCellX(false), oneCellY(false), oneCellZ(false);
    if (clist.NCell.x == 1) oneCellX = true;
    if (clist.NCell.y == 1) oneCellY = true;
    if (clist.NCell.z == 1) oneCellZ = true;

    __syncthreads();
    // find neighbor
    if (ii != MaxIndexValue){
      for (IndexType jj = 0; jj < clist.numbers[targetCellIdx]; ++jj){
	// if (targetIndexes[jj] == MaxIndexValue) break;
	ScalorType diffx = target[jj].x - shift.x - ref.x;
	ScalorType diffy = target[jj].y - shift.y - ref.y;
	ScalorType diffz = target[jj].z - shift.z - ref.z;
	if (oneCellX) shortestImage (box.size.x, box.sizei.x, &diffx);
	if (oneCellY) shortestImage (box.size.y, box.sizei.y, &diffy);
	if (oneCellZ) shortestImage (box.size.z, box.sizei.z, &diffz);
	//printf ("%d\t%d\t%f\t%f\n", ii,
	// ScalorType dr2;
	if (((diffx*diffx+diffy*diffy+diffz*diffz)) < rcut2 &&
	    targetIndexes[jj] != ii){
	  IndexType fidx(0);
	  fidx = AtomNBForceTable::calForceIndex (
	      const_nonBondedInteractionTable,
	      const_numAtomType[0],
	      reftype,
	      targettype[jj]);
	  // if (fidx != mdForceNULL) {
	  ScalorType fx, fy, fz;
	  nbForce (nonBondedInteractionType[fidx],
		   &nonBondedInteractionParameter
		   [nonBondedInteractionParameterPosition[fidx]],
		   diffx, diffy, diffz,
		   // dr2,
		   &fx, &fy, &fz);
	  fsumx += fx;
	  fsumy += fy;
	  fsumz += fz;
	  // }
	}
      }
    }
  }
  if (ii != MaxIndexValue){
    forcx[ii] += fsumx;
    forcy[ii] += fsumy;
    forcz[ii] += fsumz;
  }
}


__global__ void
calNonBondedInteraction_cell (const IndexType		numAtom,
			      const CoordType *		coord,
			      ScalorType *		forcx,
			      ScalorType *		forcy, 
			      ScalorType *		forcz,
			      const TypeType *		type,
			      const RectangularBox	box,
			      const DeviceCellList	clist,
			      const ScalorType		rcut,
			      ScalorType *		statistic_nb_buff0,
			      ScalorType *		statistic_nb_buff1,
			      ScalorType *		statistic_nb_buff2,
			      ScalorType *		statistic_nb_buff3,
			      mdError_t *		ptr_de)
{
  // RectangularBoxGeometry::normalizeSystem (box, &ddata);
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType bidx, bidy, bidz;
  D1toD3 (clist.NCell, bid, bidx, bidy, bidz);
  
  // load index
  IndexType ii = getDeviceCellListData (clist, bid, tid);
  // load iith coordinate // use texturefetch instead
  CoordType ref;
  TypeType reftype;
  ScalorType fsumx (0.f), fsumy(0.f), fsumz(0.f);
  ScalorType myPoten (0.0f), myVxx (0.0f), myVyy (0.0f), myVzz (0.0f);
  if (ii != MaxIndexValue){
#ifdef COMPILE_NO_TEX
    ref = coord[ii];
    reftype = type[ii];
#else
    ref = tex1Dfetch (global_texRef_interaction_coord, ii);
    reftype = tex1Dfetch(global_texRef_interaction_type, ii);
#endif
  }
  
  // the target index and coordinates are shared

  extern __shared__ volatile char pub_sbuff[];
  
  volatile IndexType * targetIndexes =
      (volatile IndexType *) pub_sbuff;
  CoordType * target =
      (CoordType *) &targetIndexes[roundUp4(blockDim.x)];
  volatile TypeType * targettype =
      (volatile TypeType *) &target[roundUp4(blockDim.x)];

  ScalorType rcut2 = rcut * rcut;
  for (IndexType i = 0; i < clist.numNeighborCell[bid]; ++i){
    __syncthreads();
    IndexType targetCellIdx = getNeighborCellIndex (clist, bid, i);
    CoordNoiType shiftNoi   = getNeighborCellShiftNoi (clist, bid, i);
    CoordType shift;
    shift.x = shiftNoi.x * box.size.x;
    shift.y = shiftNoi.y * box.size.y;
    shift.z = shiftNoi.z * box.size.z;
    targetIndexes[tid] = getDeviceCellListData(clist, targetCellIdx, tid);  
    if (targetIndexes[tid] != MaxIndexValue){
      target[tid] = tex1Dfetch(global_texRef_interaction_coord, targetIndexes[tid]);
      targettype[tid] = tex1Dfetch(global_texRef_interaction_type, targetIndexes[tid]);
    }
    bool oneCellX(false), oneCellY(false), oneCellZ(false);
    if (clist.NCell.x == 1) oneCellX = true;
    if (clist.NCell.y == 1) oneCellY = true;
    if (clist.NCell.z == 1) oneCellZ = true;

    __syncthreads();
    // find neighbor
    if (ii != MaxIndexValue){
      for (IndexType jj = 0; jj < clist.numbers[targetCellIdx]; ++jj){
	// if (targetIndexes[jj] == MaxIndexValue) break;
	ScalorType diffx = target[jj].x - shift.x - ref.x;
	ScalorType diffy = target[jj].y - shift.y - ref.y;
	ScalorType diffz = target[jj].z - shift.z - ref.z;
	if (oneCellX) shortestImage (box.size.x, box.sizei.x, &diffx);
	if (oneCellY) shortestImage (box.size.y, box.sizei.y, &diffy);
	if (oneCellZ) shortestImage (box.size.z, box.sizei.z, &diffz);
	//printf ("%d\t%d\t%f\t%f\n", ii,
	// ScalorType dr2;
	if (((diffx*diffx+diffy*diffy+diffz*diffz)) < rcut2 &&
	    targetIndexes[jj] != ii){
	  IndexType fidx(0);
	  fidx = AtomNBForceTable::calForceIndex (
	      const_nonBondedInteractionTable,
	      const_numAtomType[0],
	      reftype,
	      targettype[jj]);
	  // if (fidx != mdForceNULL) {
	  ScalorType fx, fy, fz, dp;
	  nbForcePoten (nonBondedInteractionType[fidx],
			&nonBondedInteractionParameter
			[nonBondedInteractionParameterPosition[fidx]],
			diffx, diffy, diffz,
			&fx, &fy, &fz, &dp);
	  myPoten += dp;
	  myVxx += fx * diffx;
	  myVyy += fy * diffy;
	  myVzz += fz * diffz;
	  fsumx += fx;
	  fsumy += fy;
	  fsumz += fz;
	  // }
	}
      }
    }
  }

  if (ii != MaxIndexValue){
    forcx[ii] += fsumx;
    forcy[ii] += fsumy;
    forcz[ii] += fsumz;
    statistic_nb_buff0[ii] = myPoten * 0.5f;
    statistic_nb_buff1[ii] = myVxx * 0.5f;
    statistic_nb_buff2[ii] = myVyy * 0.5f;
    statistic_nb_buff3[ii] = myVzz * 0.5f;
  }
}


__global__ void
calNonBondedInteraction (const IndexType		numAtom,
			 const CoordType *		coord,
			 ScalorType *			forcx,
			 ScalorType *			forcy, 
			 ScalorType *			forcz,
			 const TypeType *		type,
			 const RectangularBox		box,
			 DeviceCellList			clist,
			 const ScalorType		rcut,
			 DeviceNeighborList		nlist,
			 mdError_t *			ptr_de)
{
  // RectangularBoxGeometry::normalizeSystem (box, &ddata);
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType bidx, bidy, bidz;
  D1toD3 (clist.NCell, bid, bidx, bidy, bidz);

  IndexType Nneighbor = 0;
  // load index
  IndexType ii = getDeviceCellListData (clist, bid, tid);
  // load iith coordinate // use texturefetch instead
  CoordType ref;
  TypeType reftype;
  ScalorType fsumx (0.f), fsumy(0.f), fsumz(0.f);
  if (ii != MaxIndexValue){
#ifdef COMPILE_NO_TEX
    ref = coord[ii];
    reftype = type[ii];
#else
    ref = tex1Dfetch (global_texRef_interaction_coord, ii);
    reftype = tex1Dfetch(global_texRef_interaction_type, ii);
#endif
  }
  // ScalorType rlist = clist.rlist;

  // the target index and coordinates are shared

  extern __shared__ volatile char pub_sbuff[];
  
  volatile IndexType * targetIndexes =
      (volatile IndexType *) pub_sbuff;
  CoordType * target =
      (CoordType *) &targetIndexes[roundUp4(blockDim.x)];
  volatile TypeType * targettype =
      (volatile TypeType *) &target[roundUp4(blockDim.x)];

  ScalorType rlist2 = nlist.rlist * nlist.rlist;
  ScalorType rcut2  = rcut * rcut;
  bool oneCellX(false), oneCellY(false), oneCellZ(false);
  if (clist.NCell.x == 1) oneCellX = true;
  if (clist.NCell.y == 1) oneCellY = true;
  if (clist.NCell.z == 1) oneCellZ = true;

  for (IndexType i = 0; i < clist.numNeighborCell[bid]; ++i){
    __syncthreads();
    IndexType targetCellIdx = getNeighborCellIndex (clist, bid, i);
    CoordNoiType shiftNoi   = getNeighborCellShiftNoi (clist, bid, i);
    CoordType shift;
    shift.x = shiftNoi.x * box.size.x;
    shift.y = shiftNoi.y * box.size.y;
    shift.z = shiftNoi.z * box.size.z;
    targetIndexes[tid] = getDeviceCellListData(clist, targetCellIdx, tid);  
    if (targetIndexes[tid] != MaxIndexValue){
      target[tid] = tex1Dfetch(global_texRef_interaction_coord, targetIndexes[tid]);
      targettype[tid] = tex1Dfetch(global_texRef_interaction_type, targetIndexes[tid]);
    }

    __syncthreads();
    // find neighbor
    if (ii != MaxIndexValue){
      for (IndexType jj = 0; jj < clist.numbers[targetCellIdx]; ++jj){
	// if (targetIndexes[jj] == MaxIndexValue) break;
	ScalorType diffx = target[jj].x - shift.x - ref.x;
	ScalorType diffy = target[jj].y - shift.y - ref.y;
	ScalorType diffz = target[jj].z - shift.z - ref.z;
	if (oneCellX) shortestImage (box.size.x, box.sizei.x, &diffx);
	if (oneCellY) shortestImage (box.size.y, box.sizei.y, &diffy);
	if (oneCellZ) shortestImage (box.size.z, box.sizei.z, &diffz);
	ScalorType dr2 = (diffx*diffx+diffy*diffy+diffz*diffz);
	IndexType fidx(0);
	fidx = AtomNBForceTable::calForceIndex (
	    const_nonBondedInteractionTable,
	    const_numAtomType[0],
	    reftype,
	    targettype[jj]);
	if (dr2 < rcut2 &&
	    targetIndexes[jj] != ii){
	  // if (fidx != mdForceNULL) {
	  ScalorType fx, fy, fz, dp;
	  nbForcePoten (nonBondedInteractionType[fidx],
			&nonBondedInteractionParameter
			[nonBondedInteractionParameterPosition[fidx]],
			diffx, diffy, diffz,
			&fx, &fy, &fz, &dp);
	  fsumx += fx;
	  fsumy += fy;
	  fsumz += fz;
	  // }
	}
	if (dr2 < rlist2 &&
	    targetIndexes[jj] != ii){
	  IndexType listIdx = Nneighbor * nlist.stride + ii;
	  nlist.data[listIdx] = targetIndexes[jj];
	  nlist.forceIndex[listIdx] = fidx;
	  Nneighbor ++;
	}
      }
    }
  }

  if (ii != MaxIndexValue){
    forcx[ii] += fsumx;
    forcy[ii] += fsumy;
    forcz[ii] += fsumz;
    if (Nneighbor > nlist.listLength && ptr_de != NULL){
      *ptr_de = mdErrorShortNeighborList;
      return;
    }
    nlist.Nneighbor[ii] = Nneighbor;
  }
}


__global__ void
calNonBondedInteraction (const IndexType		numAtom,
			 const CoordType *		coord,
			 ScalorType *			forcx,
			 ScalorType *			forcy, 
			 ScalorType *			forcz,
			 const TypeType *		type,
			 const RectangularBox		box,
			 DeviceCellList			clist,
			 const ScalorType		rcut,
			 DeviceNeighborList		nlist,
			 ScalorType *			statistic_nb_buff0,
			 ScalorType *			statistic_nb_buff1,
			 ScalorType *			statistic_nb_buff2,
			 ScalorType *			statistic_nb_buff3,
			 mdError_t *			ptr_de)
{
  // RectangularBoxGeometry::normalizeSystem (box, &ddata);
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType bidx, bidy, bidz;
  D1toD3 (clist.NCell, bid, bidx, bidy, bidz);

  IndexType Nneighbor = 0;
  // load index
  IndexType ii = getDeviceCellListData (clist, bid, tid);
  // load iith coordinate // use texturefetch instead
  CoordType ref;
  TypeType reftype;
  ScalorType fsumx (0.f), fsumy(0.f), fsumz(0.f);
  ScalorType myPoten (0.0f), myVxx (0.0f), myVyy (0.0f), myVzz (0.0f);
  if (ii != MaxIndexValue){
#ifdef COMPILE_NO_TEX
    ref = coord[ii];
    reftype = type[ii];
#else
    ref = tex1Dfetch (global_texRef_interaction_coord, ii);
    reftype = tex1Dfetch(global_texRef_interaction_type, ii);
#endif
  }
  // ScalorType rlist = clist.rlist;

  // the target index and coordinates are shared

  extern __shared__ volatile char pub_sbuff[];
  
  volatile IndexType * targetIndexes =
      (volatile IndexType *) pub_sbuff;
  CoordType * target =
      (CoordType *) &targetIndexes[roundUp4(blockDim.x)];
  volatile TypeType * targettype =
      (volatile TypeType *) &target[roundUp4(blockDim.x)];

  ScalorType rlist2 = nlist.rlist * nlist.rlist;
  ScalorType rcut2  = rcut * rcut;
  bool oneCellX(false), oneCellY(false), oneCellZ(false);
  if (clist.NCell.x == 1) oneCellX = true;
  if (clist.NCell.y == 1) oneCellY = true;
  if (clist.NCell.z == 1) oneCellZ = true;

  for (IndexType i = 0; i < clist.numNeighborCell[bid]; ++i){
    __syncthreads();
    IndexType targetCellIdx = getNeighborCellIndex (clist, bid, i);
    CoordNoiType shiftNoi   = getNeighborCellShiftNoi (clist, bid, i);
    CoordType shift;
    shift.x = shiftNoi.x * box.size.x;
    shift.y = shiftNoi.y * box.size.y;
    shift.z = shiftNoi.z * box.size.z;
    targetIndexes[tid] = getDeviceCellListData(clist, targetCellIdx, tid);  
    if (targetIndexes[tid] != MaxIndexValue){
      target[tid] = tex1Dfetch(global_texRef_interaction_coord, targetIndexes[tid]);
      targettype[tid] = tex1Dfetch(global_texRef_interaction_type, targetIndexes[tid]);
    }

    __syncthreads();
    // find neighbor
    if (ii != MaxIndexValue){
      for (IndexType jj = 0; jj < clist.numbers[targetCellIdx]; ++jj){
	ScalorType diffx = target[jj].x - shift.x - ref.x;
	ScalorType diffy = target[jj].y - shift.y - ref.y;
	ScalorType diffz = target[jj].z - shift.z - ref.z;
	if (oneCellX) shortestImage (box.size.x, box.sizei.x, &diffx);
	if (oneCellY) shortestImage (box.size.y, box.sizei.y, &diffy);
	if (oneCellZ) shortestImage (box.size.z, box.sizei.z, &diffz);
	ScalorType dr2 = (diffx*diffx+diffy*diffy+diffz*diffz);
	IndexType fidx(0);
	fidx = AtomNBForceTable::calForceIndex (
	    const_nonBondedInteractionTable,
	    const_numAtomType[0],
	    reftype,
	    targettype[jj]);
	if (dr2 < rcut2 &&
	    targetIndexes[jj] != ii){
	  // if (fidx != mdForceNULL) {
	  ScalorType fx, fy, fz, dp;
	  nbForcePoten (nonBondedInteractionType[fidx],
			&nonBondedInteractionParameter
			[nonBondedInteractionParameterPosition[fidx]],
			diffx, diffy, diffz,
			&fx, &fy, &fz, &dp);
	  myPoten += dp;
	  myVxx += fx * diffx;
	  myVyy += fy * diffy;
	  myVzz += fz * diffz;
	  fsumx += fx;
	  fsumy += fy;
	  fsumz += fz;
	  // }
	}
	if (dr2 < rlist2 &&
	    targetIndexes[jj] != ii){
	  IndexType listIdx = Nneighbor * nlist.stride + ii;
	  nlist.data[listIdx] = targetIndexes[jj];
	  nlist.forceIndex[listIdx] = fidx;
	  Nneighbor ++;
	}
      }
    }
  }

  if (ii != MaxIndexValue){
    forcx[ii] += fsumx;
    forcy[ii] += fsumy;
    forcz[ii] += fsumz;
    statistic_nb_buff0[ii] = myPoten * 0.5f;
    statistic_nb_buff1[ii] = myVxx * 0.5f;
    statistic_nb_buff2[ii] = myVyy * 0.5f;
    statistic_nb_buff3[ii] = myVzz * 0.5f;
    if (Nneighbor > nlist.listLength && ptr_de != NULL){
      *ptr_de = mdErrorShortNeighborList;
      return;
    }
    nlist.Nneighbor[ii] = Nneighbor;
  }
}






__global__ void
calNonBondedInteraction_all (const IndexType		numAtom,
			     const CoordType *		coord,
			     ScalorType *		forcx,
			     ScalorType *		forcy, 
			     ScalorType *		forcz,
			     const TypeType *		type,
			     const RectangularBox	box,
			     const ScalorType		rcut,
			     mdError_t *		ptr_de)
{
  // RectangularBoxGeometry::normalizeSystem (box, &ddata);
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  
  IndexType numberAtom = numAtom;
  IndexType ii = tid + bid * blockDim.x;
  ScalorType fsumx (0.f), fsumy(0.f), fsumz(0.f);

  extern __shared__ volatile char pub_sbuff[];

  volatile CoordType * target =
      (volatile CoordType *) pub_sbuff;
  volatile TypeType * targettype =
      (volatile TypeType *) &target[roundUp4(blockDim.x)];

  __syncthreads();
  
  CoordType ref;
  TypeType reftype;
  if (ii < numberAtom){
    ref = coord[ii];
    reftype = type[ii];
  }
  ScalorType rcut2 = rcut * rcut;
  
  for (IndexType targetBlockId = 0;
       targetBlockId * blockDim.x < numberAtom; ++targetBlockId){
    IndexType jj = tid + targetBlockId * blockDim.x;
    __syncthreads();
    if (jj < numberAtom){
      target[tid].x = coord[jj].x;
      target[tid].y = coord[jj].y;
      target[tid].z = coord[jj].z;
      targettype[tid] = type[jj];
    }
    __syncthreads();
    if (ii < numberAtom){
      for (IndexType kk = 0; kk < blockDim.x; ++kk){
	if (kk + targetBlockId * blockDim.x >= numberAtom) break;
	ScalorType diffx = target[kk].x - ref.x;
	ScalorType diffy = target[kk].y - ref.y;
	ScalorType diffz = target[kk].z - ref.z;
	shortestImage (box, &diffx, &diffy, &diffz);
	ScalorType dr2;
	if ((dr2 = diffx*diffx+diffy*diffy+diffz*diffz) < rcut2 &&
	    kk + targetBlockId * blockDim.x != ii){
	  IndexType fidx = AtomNBForceTable::calForceIndex (
	      const_nonBondedInteractionTable,
	      const_numAtomType[0],
	      reftype,
	      targettype[kk]);
	  if (dr2 < rcut2 ) {
	    // if (fidx != mdForceNULL) {
	    ScalorType fx, fy, fz, dp;
	    nbForcePoten (nonBondedInteractionType[fidx],
			  &nonBondedInteractionParameter
			  [nonBondedInteractionParameterPosition[fidx]],
			  diffx, diffy, diffz,
			  &fx, &fy, &fz, &dp);
	    fsumx += fx;
	    fsumy += fy;
	    fsumz += fz;
	    // }
	  }
	}
      }
    }
  }
  if (ii < numberAtom){
    forcx[ii] += fsumx;
    forcy[ii] += fsumy;
    forcz[ii] += fsumz;
  }
}


__global__ void
calNonBondedInteraction_all  (const IndexType		numAtom,
			      const CoordType *		coord,
			      ScalorType *		forcx,
			      ScalorType *		forcy, 
			      ScalorType *		forcz,
			      const TypeType *		type,
			      const RectangularBox	box,
			      const ScalorType		rcut,
			      ScalorType *		statistic_nb_buff0,
			      ScalorType *		statistic_nb_buff1,
			      ScalorType *		statistic_nb_buff2,
			      ScalorType *		statistic_nb_buff3,
			      mdError_t *		ptr_de)
{
  // RectangularBoxGeometry::normalizeSystem (box, &ddata);
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  
  IndexType numberAtom = numAtom;
  IndexType ii = tid + bid * blockDim.x;
  ScalorType fsumx (0.f), fsumy(0.f), fsumz(0.f);
  ScalorType myPoten (0.0f), myVxx (0.0f), myVyy (0.0f), myVzz (0.0f);

  extern __shared__ volatile char pub_sbuff[];

  volatile CoordType * target =
      (volatile CoordType *) pub_sbuff;
  volatile TypeType * targettype =
      (volatile TypeType *) &target[roundUp4(blockDim.x)];

  __syncthreads();
  
  CoordType ref;
  TypeType reftype;
  if (ii < numberAtom){
    ref = coord[ii];
    reftype = type[ii];
  }
  ScalorType rcut2 = rcut * rcut;
  
  for (IndexType targetBlockId = 0;
       targetBlockId * blockDim.x < numberAtom; ++targetBlockId){
    IndexType jj = tid + targetBlockId * blockDim.x;
    __syncthreads();
    if (jj < numberAtom){
      target[tid].x = coord[jj].x;
      target[tid].y = coord[jj].y;
      target[tid].z = coord[jj].z;
      targettype[tid] = type[jj];
    }
    __syncthreads();
    if (ii < numberAtom){
      for (IndexType kk = 0; kk < blockDim.x; ++kk){
	if (kk + targetBlockId * blockDim.x >= numberAtom) break;
	ScalorType diffx = target[kk].x - ref.x;
	ScalorType diffy = target[kk].y - ref.y;
	ScalorType diffz = target[kk].z - ref.z;
	shortestImage (box, &diffx, &diffy, &diffz);
	ScalorType dr2;
	if ((dr2 = diffx*diffx+diffy*diffy+diffz*diffz) < rcut2 &&
	    kk + targetBlockId * blockDim.x != ii){
	  IndexType fidx = AtomNBForceTable::calForceIndex (
	      const_nonBondedInteractionTable,
	      const_numAtomType[0],
	      reftype,
	      targettype[kk]);
	  if (dr2 < rcut2 ) {
	    // if (fidx != mdForceNULL) {
	    ScalorType fx, fy, fz, dp;
	    nbForcePoten (nonBondedInteractionType[fidx],
			  &nonBondedInteractionParameter
			  [nonBondedInteractionParameterPosition[fidx]],
			  diffx, diffy, diffz,
			  &fx, &fy, &fz, &dp);
	    myPoten += dp;
	    myVxx += fx * diffx;
	    myVyy += fy * diffy;
	    myVzz += fz * diffz;
	    fsumx += fx;
	    fsumy += fy;
	    fsumz += fz;
	    // }
	  }
	}
      }
    }
  }
  if (ii < numberAtom){
    forcx[ii] += fsumx;
    forcy[ii] += fsumy;
    forcz[ii] += fsumz;
    statistic_nb_buff0[ii] = myPoten * 0.5f;
    statistic_nb_buff1[ii] = myVxx * 0.5f;
    statistic_nb_buff2[ii] = myVyy * 0.5f;
    statistic_nb_buff3[ii] = myVzz * 0.5f;
  }
}



__global__ void
calTwinRangeCorrection_cell (const IndexType		numAtom,
			     const CoordType *		coord,
			     ScalorType *		forcx,
			     ScalorType *		forcy, 
			     ScalorType *		forcz,
			     const TypeType *		type,
			     const RectangularBox	box,
			     const DeviceCellList	clist,
			     const ScalorType		rcut1,
			     const ScalorType		rcut2,
			     ScalorType *		statistic_nb_buff0,
			     ScalorType *		statistic_nb_buff1,
			     ScalorType *		statistic_nb_buff2,
			     ScalorType *		statistic_nb_buff3,
			     mdError_t *		ptr_de)
{
  // RectangularBoxGeometry::normalizeSystem (box, &ddata);
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType bidx, bidy, bidz;
  D1toD3 (clist.NCell, bid, bidx, bidy, bidz);
  
  // load index
  IndexType ii = getDeviceCellListData (clist, bid, tid);
  // load iith coordinate // use texturefetch instead
  CoordType ref;
  TypeType reftype;
  ScalorType fsumx (0.f), fsumy(0.f), fsumz(0.f);
  ScalorType myPoten (0.0f), myVxx (0.0f), myVyy (0.0f), myVzz (0.0f);
  if (ii != MaxIndexValue){
#ifdef COMPILE_NO_TEX
    ref = coord[ii];
    reftype = type[ii];
#else
    ref = tex1Dfetch (global_texRef_interaction_coord, ii);
    reftype = tex1Dfetch(global_texRef_interaction_type, ii);
#endif
  }
  
  // the target index and coordinates are shared

  extern __shared__ volatile char pub_sbuff[];
  
  volatile IndexType * targetIndexes =
      (volatile IndexType *) pub_sbuff;
  CoordType * target =
      (CoordType *) &targetIndexes[roundUp4(blockDim.x)];
  volatile TypeType * targettype =
      (volatile TypeType *) &target[roundUp4(blockDim.x)];

  ScalorType rcut12 = rcut1 * rcut1;
  ScalorType rcut22 = rcut2 * rcut2;
  for (IndexType i = 0; i < clist.numNeighborCell[bid]; ++i){
    __syncthreads();
    IndexType targetCellIdx = getNeighborCellIndex (clist, bid, i);
    CoordNoiType shiftNoi   = getNeighborCellShiftNoi (clist, bid, i);
    CoordType shift;
    shift.x = shiftNoi.x * box.size.x;
    shift.y = shiftNoi.y * box.size.y;
    shift.z = shiftNoi.z * box.size.z;
    targetIndexes[tid] = getDeviceCellListData(clist, targetCellIdx, tid);  
    if (targetIndexes[tid] != MaxIndexValue){
      target[tid] = tex1Dfetch(global_texRef_interaction_coord, targetIndexes[tid]);
      targettype[tid] = tex1Dfetch(global_texRef_interaction_type, targetIndexes[tid]);
    }
    bool oneCellX(false), oneCellY(false), oneCellZ(false);
    if (clist.NCell.x == 1) oneCellX = true;
    if (clist.NCell.y == 1) oneCellY = true;
    if (clist.NCell.z == 1) oneCellZ = true;

    __syncthreads();
    // find neighbor
    if (ii != MaxIndexValue){
      for (IndexType jj = 0; jj < clist.numbers[targetCellIdx]; ++jj){
	// if (targetIndexes[jj] == MaxIndexValue) break;
	ScalorType diffx = target[jj].x - shift.x - ref.x;
	ScalorType diffy = target[jj].y - shift.y - ref.y;
	ScalorType diffz = target[jj].z - shift.z - ref.z;
	if (oneCellX) shortestImage (box.size.x, box.sizei.x, &diffx);
	if (oneCellY) shortestImage (box.size.y, box.sizei.y, &diffy);
	if (oneCellZ) shortestImage (box.size.z, box.sizei.z, &diffz);
	ScalorType dr2 = (diffx*diffx+diffy*diffy+diffz*diffz);
	if (dr2 < rcut22 && dr2 >= rcut12 &&
	    targetIndexes[jj] != ii){
	  IndexType fidx(0);
	  fidx = AtomNBForceTable::calForceIndex (
	      const_nonBondedInteractionTable,
	      const_numAtomType[0],
	      reftype,
	      targettype[jj]);
	  // if (fidx != mdForceNULL) {
	  ScalorType fx, fy, fz, dp;
	  nbForcePoten (nonBondedInteractionType[fidx],
			&nonBondedInteractionParameter
			[nonBondedInteractionParameterPosition[fidx]],
			diffx, diffy, diffz,
			&fx, &fy, &fz, &dp);
	  // printf ("# %d\t%d\t%f\t%f\t%f\n",
	  // 	  ii, targetIndexes[jj],
	  // 	  ref.z, target[jj].z, fz);
	  myPoten += dp;
	  myVxx += fx * diffx;
	  myVyy += fy * diffy;
	  myVzz += fz * diffz;
	  fsumx += fx;
	  fsumy += fy;
	  fsumz += fz;
	  // }
	}
      }
    }
  }

  if (ii != MaxIndexValue){
    forcx[ii] = fsumx;
    forcy[ii] = fsumy;
    forcz[ii] = fsumz;
    statistic_nb_buff0[ii] = myPoten * 0.5f;
    statistic_nb_buff1[ii] = myVxx * 0.5f;
    statistic_nb_buff2[ii] = myVyy * 0.5f;
    statistic_nb_buff3[ii] = myVzz * 0.5f;
  }
}


__global__ void
calTwinRangeCorrection_all   (const IndexType		numAtom,
			      const CoordType *		coord,
			      ScalorType *		forcx,
			      ScalorType *		forcy, 
			      ScalorType *		forcz,
			      const TypeType *		type,
			      const RectangularBox	box,
			      const ScalorType		rcut1,
			      const ScalorType		rcut2,
			      ScalorType *		statistic_nb_buff0,
			      ScalorType *		statistic_nb_buff1,
			      ScalorType *		statistic_nb_buff2,
			      ScalorType *		statistic_nb_buff3,
			      mdError_t *		ptr_de)
{
  // RectangularBoxGeometry::normalizeSystem (box, &ddata);
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  
  IndexType numberAtom = numAtom;
  IndexType ii = tid + bid * blockDim.x;
  ScalorType fsumx (0.f), fsumy(0.f), fsumz(0.f);
  ScalorType myPoten (0.0f), myVxx (0.0f), myVyy (0.0f), myVzz (0.0f);

  extern __shared__ volatile char pub_sbuff[];

  volatile CoordType * target =
      (volatile CoordType *) pub_sbuff;
  volatile TypeType * targettype =
      (volatile TypeType *) &target[roundUp4(blockDim.x)];

  __syncthreads();
  
  CoordType ref;
  TypeType reftype;
  if (ii < numberAtom){
    ref = coord[ii];
    reftype = type[ii];
  }
  ScalorType rcut12 = rcut1 * rcut1;
  ScalorType rcut22 = rcut2 * rcut2;
  
  for (IndexType targetBlockId = 0;
       targetBlockId * blockDim.x < numberAtom; ++targetBlockId){
    IndexType jj = tid + targetBlockId * blockDim.x;
    __syncthreads();
    if (jj < numberAtom){
      target[tid].x = coord[jj].x;
      target[tid].y = coord[jj].y;
      target[tid].z = coord[jj].z;
      targettype[tid] = type[jj];
    }
    __syncthreads();
    if (ii < numberAtom){
      for (IndexType kk = 0; kk < blockDim.x; ++kk){
	if (kk + targetBlockId * blockDim.x >= numberAtom) break;
	ScalorType diffx = target[kk].x - ref.x;
	ScalorType diffy = target[kk].y - ref.y;
	ScalorType diffz = target[kk].z - ref.z;
	shortestImage (box, &diffx, &diffy, &diffz);
	ScalorType dr2 = diffx*diffx+diffy*diffy+diffz*diffz;
	if (dr2 < rcut22 && dr2 >= rcut12 &&
	    kk + targetBlockId * blockDim.x != ii){
	  IndexType fidx = AtomNBForceTable::calForceIndex (
	      const_nonBondedInteractionTable,
	      const_numAtomType[0],
	      reftype,
	      targettype[kk]);
	  // if (fidx != mdForceNULL) {
	  ScalorType fx, fy, fz, dp;
	  nbForcePoten (nonBondedInteractionType[fidx],
			&nonBondedInteractionParameter
			[nonBondedInteractionParameterPosition[fidx]],
			diffx, diffy, diffz,
			&fx, &fy, &fz, &dp);
	  myPoten += dp;
	  myVxx += fx * diffx;
	  myVyy += fy * diffy;
	  myVzz += fz * diffz;
	  fsumx += fx;
	  fsumy += fy;
	  fsumz += fz;
	  // }
	}
      }
    }
  }
  if (ii < numberAtom){
    forcx[ii] = fsumx;
    forcy[ii] = fsumy;
    forcz[ii] = fsumz;
    statistic_nb_buff0[ii] = myPoten * 0.5f;
    statistic_nb_buff1[ii] = myVxx * 0.5f;
    statistic_nb_buff2[ii] = myVyy * 0.5f;
    statistic_nb_buff3[ii] = myVzz * 0.5f;
  }
}



__global__ void
buildNeighborListCalTwinRangeCorr_cell (const IndexType		numAtom,
					const CoordType *	coord,
					ScalorType *		forcx,
					ScalorType *		forcy, 
					ScalorType *		forcz,
					const TypeType *	type,
					const RectangularBox	box,
					const DeviceCellList	clist,
					const ScalorType	rcut1,
					const ScalorType	rcut2,
					DeviceNeighborList	nlist,
					ScalorType *		statistic_nb_buff0,
					ScalorType *		statistic_nb_buff1,
					ScalorType *		statistic_nb_buff2,
					ScalorType *		statistic_nb_buff3,
					mdError_t *		ptr_de)
{
  // RectangularBoxGeometry::normalizeSystem (box, &ddata);
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType bidx, bidy, bidz;
  D1toD3 (clist.NCell, bid, bidx, bidy, bidz);
  
  // set number of neighbor to 0
  IndexType Nneighbor = 0;
  // load index
  IndexType ii = getDeviceCellListData (clist, bid, tid);
  // load iith coordinate // use texturefetch instead
  CoordType ref;
  TypeType reftype;
  ScalorType fsumx (0.f), fsumy(0.f), fsumz(0.f);
  ScalorType myPoten (0.0f), myVxx (0.0f), myVyy (0.0f), myVzz (0.0f);
  if (ii != MaxIndexValue){
#ifdef COMPILE_NO_TEX
    ref = coord[ii];
    reftype = type[ii];
#else
    ref = tex1Dfetch (global_texRef_interaction_coord, ii);
    reftype = tex1Dfetch(global_texRef_interaction_type, ii);
#endif
  }
  
  // the target index and coordinates are shared

  extern __shared__ volatile char pub_sbuff[];
  
  volatile IndexType * targetIndexes =
      (volatile IndexType *) pub_sbuff;
  CoordType * target =
      (CoordType *) &targetIndexes[roundUp4(blockDim.x)];
  volatile TypeType * targettype =
      (volatile TypeType *) &target[roundUp4(blockDim.x)];

  ScalorType rcut12 = rcut1 * rcut1;
  ScalorType rcut22 = rcut2 * rcut2;
  for (IndexType i = 0; i < clist.numNeighborCell[bid]; ++i){
    __syncthreads();
    IndexType targetCellIdx = getNeighborCellIndex (clist, bid, i);
    CoordNoiType shiftNoi   = getNeighborCellShiftNoi (clist, bid, i);
    CoordType shift;
    shift.x = shiftNoi.x * box.size.x;
    shift.y = shiftNoi.y * box.size.y;
    shift.z = shiftNoi.z * box.size.z;
    targetIndexes[tid] = getDeviceCellListData(clist, targetCellIdx, tid);  
    if (targetIndexes[tid] != MaxIndexValue){
      target[tid] = tex1Dfetch(global_texRef_interaction_coord, targetIndexes[tid]);
      targettype[tid] = tex1Dfetch(global_texRef_interaction_type, targetIndexes[tid]);
    }
    bool oneCellX(false), oneCellY(false), oneCellZ(false);
    if (clist.NCell.x == 1) oneCellX = true;
    if (clist.NCell.y == 1) oneCellY = true;
    if (clist.NCell.z == 1) oneCellZ = true;

    __syncthreads();
    // find neighbor
    if (ii != MaxIndexValue){
      for (IndexType jj = 0; jj < clist.numbers[targetCellIdx]; ++jj){
	// if (targetIndexes[jj] == MaxIndexValue) break;
	ScalorType diffx = target[jj].x - shift.x - ref.x;
	ScalorType diffy = target[jj].y - shift.y - ref.y;
	ScalorType diffz = target[jj].z - shift.z - ref.z;
	if (oneCellX) shortestImage (box.size.x, box.sizei.x, &diffx);
	if (oneCellY) shortestImage (box.size.y, box.sizei.y, &diffy);
	if (oneCellZ) shortestImage (box.size.z, box.sizei.z, &diffz);
	ScalorType dr2 = (diffx*diffx+diffy*diffy+diffz*diffz);
	if (targetIndexes[jj] != ii){
	  if (dr2 < rcut22 && dr2 >= rcut12 ){
	    IndexType fidx(0);
	    fidx = AtomNBForceTable::calForceIndex (
		const_nonBondedInteractionTable,
		const_numAtomType[0],
		reftype,
		targettype[jj]);
	    // if (fidx != mdForceNULL) {
	    ScalorType fx, fy, fz, dp;
	    nbForcePoten (nonBondedInteractionType[fidx],
			  &nonBondedInteractionParameter
			  [nonBondedInteractionParameterPosition[fidx]],
			  diffx, diffy, diffz,
			  &fx, &fy, &fz, &dp);
	    // printf ("# %d\t%d\t%f\t%f\t%f\n",
	    // 	  ii, targetIndexes[jj],
	    // 	  ref.z, target[jj].z, fz);
	    myPoten += dp;
	    myVxx += fx * diffx;
	    myVyy += fy * diffy;
	    myVzz += fz * diffz;
	    fsumx += fx;
	    fsumy += fy;
	    fsumz += fz;
	    // }
	  }
	  else if (dr2 < rcut12){
	    IndexType fidx(0);
	    fidx = AtomNBForceTable::calForceIndex (
		const_nonBondedInteractionTable,
		const_numAtomType[0],
		reftype,
		targettype[jj]);
	    IndexType listIdx = Nneighbor * nlist.stride + ii;
	    nlist.data[listIdx] = targetIndexes[jj];
	    nlist.forceIndex[listIdx] = fidx;
	    Nneighbor ++;
	  }
	}  
      }
    }
  }

  if (ii != MaxIndexValue){
    forcx[ii] = fsumx;
    forcy[ii] = fsumy;
    forcz[ii] = fsumz;
    statistic_nb_buff0[ii] = myPoten * 0.5f;
    statistic_nb_buff1[ii] = myVxx * 0.5f;
    statistic_nb_buff2[ii] = myVyy * 0.5f;
    statistic_nb_buff3[ii] = myVzz * 0.5f;
    if (Nneighbor > nlist.listLength && ptr_de != NULL){
      *ptr_de = mdErrorShortNeighborList;
    }
    nlist.Nneighbor[ii] = Nneighbor;
  }
}


__global__ void
buildNeighborListCalTwinRangeCorr_all (const IndexType		numAtom,
				       const CoordType *	coord,
				       ScalorType *		forcx,
				       ScalorType *		forcy, 
				       ScalorType *		forcz,
				       const TypeType *		type,
				       const RectangularBox	box,
				       const ScalorType		rcut1,
				       const ScalorType		rcut2,
				       DeviceNeighborList	nlist,
				       ScalorType *		statistic_nb_buff0,
				       ScalorType *		statistic_nb_buff1,
				       ScalorType *		statistic_nb_buff2,
				       ScalorType *		statistic_nb_buff3,
				       mdError_t *		ptr_de)
{
  // RectangularBoxGeometry::normalizeSystem (box, &ddata);
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  
  IndexType numberAtom = numAtom;
  IndexType Nneighbor = 0;
  IndexType ii = tid + bid * blockDim.x;
  ScalorType fsumx (0.f), fsumy(0.f), fsumz(0.f);
  ScalorType myPoten (0.0f), myVxx (0.0f), myVyy (0.0f), myVzz (0.0f);

  extern __shared__ volatile char pub_sbuff[];

  volatile CoordType * target =
      (volatile CoordType *) pub_sbuff;
  volatile TypeType * targettype =
      (volatile TypeType *) &target[roundUp4(blockDim.x)];

  __syncthreads();
  
  CoordType ref;
  TypeType reftype;
  if (ii < numberAtom){
    ref = coord[ii];
    reftype = type[ii];
  }
  ScalorType rcut12 = rcut1 * rcut1;
  ScalorType rcut22 = rcut2 * rcut2;
  
  for (IndexType targetBlockId = 0;
       targetBlockId * blockDim.x < numberAtom; ++targetBlockId){
    IndexType jj = tid + targetBlockId * blockDim.x;
    __syncthreads();
    if (jj < numberAtom){
      target[tid].x = coord[jj].x;
      target[tid].y = coord[jj].y;
      target[tid].z = coord[jj].z;
      targettype[tid] = type[jj];
    }
    __syncthreads();
    if (ii < numberAtom){
      for (IndexType kk = 0; kk < blockDim.x; ++kk){
	if (kk + targetBlockId * blockDim.x >= numberAtom) break;
	ScalorType diffx = target[kk].x - ref.x;
	ScalorType diffy = target[kk].y - ref.y;
	ScalorType diffz = target[kk].z - ref.z;
	shortestImage (box, &diffx, &diffy, &diffz);
	ScalorType dr2 = diffx*diffx+diffy*diffy+diffz*diffz;
	if (kk + targetBlockId * blockDim.x != ii){
	  if (dr2 < rcut22 && dr2 >= rcut12 ){   
	    IndexType fidx = AtomNBForceTable::calForceIndex (
		const_nonBondedInteractionTable,
		const_numAtomType[0],
		reftype,
		targettype[kk]);
	    // if (fidx != mdForceNULL) {
	    ScalorType fx, fy, fz, dp;
	    nbForcePoten (nonBondedInteractionType[fidx],
			  &nonBondedInteractionParameter
			  [nonBondedInteractionParameterPosition[fidx]],
			  diffx, diffy, diffz,
			  &fx, &fy, &fz, &dp);
	    myPoten += dp;
	    myVxx += fx * diffx;
	    myVyy += fy * diffy;
	    myVzz += fz * diffz;
	    fsumx += fx;
	    fsumy += fy;
	    fsumz += fz;
	    // }
	  }
	  else if (dr2 < rcut12){
	    IndexType fidx = AtomNBForceTable::calForceIndex (
		const_nonBondedInteractionTable,
		const_numAtomType[0],
		reftype,
		targettype[kk]);
	    IndexType listIdx = Nneighbor * nlist.stride + ii;
	    nlist.data[listIdx] = kk + targetBlockId * blockDim.x;
	    nlist.forceIndex[listIdx] = fidx;
	    Nneighbor ++;
	  }
	}
      }
    }
  }
  if (ii < numberAtom){
    forcx[ii] = fsumx;
    forcy[ii] = fsumy;
    forcz[ii] = fsumz;
    statistic_nb_buff0[ii] = myPoten * 0.5f;
    statistic_nb_buff1[ii] = myVxx * 0.5f;
    statistic_nb_buff2[ii] = myVyy * 0.5f;
    statistic_nb_buff3[ii] = myVzz * 0.5f;
    if (Nneighbor > nlist.listLength && ptr_de != NULL){
      *ptr_de = mdErrorShortNeighborList;
    }
    nlist.Nneighbor[ii] = Nneighbor;
  }
}






__global__ void
widomDeltaPoten_NVT (const IndexType		numTestParticle,
		     const CoordType *		coordTestParticle,
		     const TypeType *		typeTestParticle,
		     const IndexType		numAtom,
		     const CoordType *		coord,
		     const TypeType *		type,
		     const RectangularBox	box,
		     DeviceCellList		clist,
		     ScalorType *		statistic_nb_buff0,
		     mdError_t *		ptr_de)
{
  // RectangularBoxGeometry::normalizeSystem (box, &ddata);
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
//  IndexType ii = tid + bid * blockDim.x;
  if (bid >= numTestParticle) return;

  // extern __shared__ volatile char pub_sbuff_widom[];
  // volatile ScalorType * sumbuff = (volatile ScalorType *) pub_sbuff_widom;
  extern __shared__ volatile ScalorType sumbuff [];

  CoordType refCoord = coordTestParticle[bid];
  TypeType  refType = typeTestParticle[bid];
  ScalorType myPoten (0.0f);

  IndexType refCelli, refCellj, refCellk;
  refCelli = IndexType (refCoord.x * box.sizei.x * ScalorType(clist.NCell.x));
  refCellj = IndexType (refCoord.y * box.sizei.y * ScalorType(clist.NCell.y));
  refCellk = IndexType (refCoord.z * box.sizei.z * ScalorType(clist.NCell.z));
  
  if (refCelli == clist.NCell.x){
    refCelli -= clist.NCell.x;
    refCoord.x -= box.size.x;
  }
  if (refCellj == clist.NCell.y){
    refCellj -= clist.NCell.y;
    refCoord.y -= box.size.y;
  }
  if (refCellk == clist.NCell.z){
    refCellk -= clist.NCell.z;
    refCoord.z -= box.size.z;
  }

  IndexType refCellIndex = D3toD1 (clist.NCell, refCelli, refCellj, refCellk);
  for (IndexType i = 0; i < clist.numNeighborCell[refCellIndex]; ++i){
    __syncthreads ();
    IndexType targetCellIdx = getNeighborCellIndex    (clist, refCellIndex, i);
    CoordNoiType shiftNoi   = getNeighborCellShiftNoi (clist, refCellIndex, i);
    CoordType shift;
    shift.x = shiftNoi.x * box.size.x;
    shift.y = shiftNoi.y * box.size.y;
    shift.z = shiftNoi.z * box.size.z;
    IndexType targetIndex = getDeviceCellListData(clist, targetCellIdx, tid);  
    if (targetIndex != MaxIndexValue){
      TypeType targettype = tex1Dfetch(global_texRef_interaction_type, targetIndex);
      if (refType == targettype){
	CoordType targetCoord = tex1Dfetch(global_texRef_interaction_coord, targetIndex);
	ScalorType diffx = targetCoord.x - shift.x - refCoord.x;
	ScalorType diffy = targetCoord.y - shift.y - refCoord.y;
	ScalorType diffz = targetCoord.z - shift.z - refCoord.z;
	ScalorType dr2 = ((diffx*diffx+diffy*diffy+diffz*diffz));
	if (dr2 < clist.rlist*clist.rlist && dr2 > 1e-4){
	  IndexType fidx(0);
	  ScalorType dp;
	  fidx = AtomNBForceTable::
	      calForceIndex (const_nonBondedInteractionTable,
			     const_numAtomType[0],
			     refType,
			     refType);
	  nbPoten (nonBondedInteractionType[fidx],
		   &nonBondedInteractionParameter
		   [nonBondedInteractionParameterPosition[fidx]],
		   diffx, diffy, diffz, &dp);
	  myPoten += dp;
	  // printf ("dp: %f,  %f %f %f\n", dp, diffx, diffy, diffz);
	}
      }
    }
  }

  sumbuff[tid] = myPoten;
  __syncthreads();
  sumVectorBlockBuffer_2 (sumbuff);
  __syncthreads();
  if (tid == 0){
    statistic_nb_buff0[bid] = sumbuff[0];
  }
}
  
//   if (tid == 0){
//     // printf ("### du is %f\n", sumbuff[0]);
//     statistic_nb_buff0[bid] = expf(- (sumbuff[0] + energyCorrection) / temperature);
//   }
// }




__global__ void
widomDeltaPoten_allPair_NVT (const IndexType		numTestParticle,
			     const CoordType *		coordTestParticle,
			     const TypeType *		typeTestParticle,
			     const IndexType		numAtom,
			     const CoordType *		coord,
			     const TypeType *		type,
			     const RectangularBox	box,
			     const ScalorType		rlist,
			     ScalorType *		statistic_nb_buff0,
			     mdError_t *		ptr_de)
{
  // RectangularBoxGeometry::normalizeSystem (box, &ddata);
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  // IndexType ii = tid + bid * blockDim.x;

  if (bid >= numTestParticle) return;
  
  CoordType refCoord = coordTestParticle[bid];
  TypeType refType = typeTestParticle[bid];

  ScalorType myPoten = 0.;
  extern __shared__ volatile ScalorType sumbuff [];
  
  for (IndexType start = 0; start < numAtom; start += blockDim.x){
    IndexType targetIndex = start + tid;
    if (targetIndex >= numAtom) break;
    TypeType targetType = type[targetIndex];
    if (targetType != refType) continue;
    CoordType targetCoord = coord[targetIndex];
    ScalorType diffx = targetCoord.x - refCoord.x;
    ScalorType diffy = targetCoord.y - refCoord.y;
    ScalorType diffz = targetCoord.z - refCoord.z;
    RectangularBoxGeometry::shortestImage (box, &diffx, &diffy, &diffz);
    ScalorType dr2 = (diffx*diffx+diffy*diffy+diffz*diffz);
    if (dr2 < rlist * rlist && dr2 > 1e-4 ){
      IndexType fidx(0);
      ScalorType dp;
      fidx = AtomNBForceTable::
	  calForceIndex (const_nonBondedInteractionTable,
			 const_numAtomType[0],
			 refType,
			 refType);
      nbPoten (nonBondedInteractionType[fidx],
	       &nonBondedInteractionParameter
	       [nonBondedInteractionParameterPosition[fidx]],
	       diffx, diffy, diffz, &dp);
      myPoten += dp;
    }
  }

  sumbuff[tid] = myPoten;
  __syncthreads();
  sumVectorBlockBuffer_2 (sumbuff);
  __syncthreads();
  if (tid == 0){
    statistic_nb_buff0[bid] = sumbuff[0];
  }
}

//   if (tid == 0){
//     // printf ("### du is %f\n", sumbuff[0]);
//     statistic_nb_buff0[bid] = expf(- (sumbuff[0] + energyCorrection) / temperature);
//   }  
// }

