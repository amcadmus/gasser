#ifndef __Parallel_MDData_h_wanghan__
#define __Parallel_MDData_h_wanghan__

#include "common.h"
#include "Topology.h"
#include "BoxGeometry.h"
#include "Parallel_DataTransferBlock.h"

extern "C"{
#include "xdrfile/xdrfile.h"
#include "xdrfile/xdrfile_xtc.h"
}

using namespace RectangularBoxGeometry;

namespace Parallel{
  class MDExcptInconsistentMemorySizeOnHostMDData : public MDException {};
  class MDExcptNumAtomMoreThanMemSize : public MDException {};
  class MDExcptInvalidTopology : public MDException {};

  class HostMDData ;
  class GlobalHostMDData;
  class DeviceMDData;

  enum MDDataItemShift
  {
    MDDataItemShift_Coordinate		= 0,
    MDDataItemShift_CoordinateNoi	= 1,
    MDDataItemShift_Velocity		= 2,
    MDDataItemShift_Force		= 3,
    MDDataItemShift_GlobalIndex		= 4,
    MDDataItemShift_Type		= 5,
    MDDataItemShift_Mass		= 6,
    MDDataItemShift_Charge		= 7
  } ;
  enum MDDataItemMask
  {
    MDDataItemMask_All			= MaxIndexValue,
    MDDataItemMask_None			= 0,
    MDDataItemMask_Coordinate		= (1<<MDDataItemShift_Coordinate),
    MDDataItemMask_CoordinateNoi	= (1<<MDDataItemShift_CoordinateNoi),
    MDDataItemMask_Velocity		= (1<<MDDataItemShift_Velocity),
    MDDataItemMask_Force		= (1<<MDDataItemShift_Force),
    MDDataItemMask_AllExceptForce	= (MDDataItemMask_All^MDDataItemMask_Force),
    MDDataItemMask_GlobalIndex		= (1<<MDDataItemShift_GlobalIndex),
    MDDataItemMask_Type			= (1<<MDDataItemShift_Type),
    MDDataItemMask_Mass			= (1<<MDDataItemShift_Mass),
    MDDataItemMask_Charge		= (1<<MDDataItemShift_Charge)
  };
  typedef IndexType MDDataItemMask_t;
  
  class HostMDData 
  {
    friend class DeviceMDData ;
protected:
    IndexType numData_;
    IndexType memSize_;
    HostCoordType * coord;
    HostCoordNoiType * coordNoi;
    ScalorType * velox;
    ScalorType * veloy;
    ScalorType * veloz;
    ScalorType * forcx;
    ScalorType * forcy;
    ScalorType * forcz;
    IndexType  * globalIndex;
    RectangularBox globalBox;
    // RectangularBox localBox;
protected:
    // topology related
    TypeType   * type;
    ScalorType * mass;
    ScalorType * charge;
protected:
    void reallocCoordinate  (const IndexType & memSize);
    void reallocVelocity    (const IndexType & memSize);
    void reallocForce       (const IndexType & memSize);
    void reallocGlobalIndex (const IndexType & memSize);
    void reallocTopProperty (const IndexType & memSize);
public:
    HostMDData ();
    HostMDData (const HostMDData & hdata);
    ~HostMDData ();
public:
    const IndexType & numData () const {return numData_;}
    const IndexType & memSize () const {return memSize_;}
    IndexType & numData ()  {return numData_;}
    // IndexType & memSize ()  {return memSize_;}
    void clear ();
    void clearData () {numData_ = 0;}
    void easyRealloc (const IndexType & memSize);
    void copy (const HostMDData & hdata,
	       const MDDataItemMask_t mask = MDDataItemMask_All);
public:
    const RectangularBox & getGlobalBox	() const {return globalBox;}
    void setGlobalBox (const ScalorType & bx,
		       const ScalorType & by,
		       const ScalorType & bz)
	{ setBoxSize (bx, by, bz, &globalBox); }
    void setGlobalBox (const RectangularBox & box) { globalBox = box; }
    void pushBackAtom (const HostCoordType & coord,
		       const HostCoordNoiType & coordNoi,
		       const ScalorType & velox,
		       const ScalorType & veloy,
		       const ScalorType & veloz,
		       const IndexType & globalIndex,
		       const TypeType & type,
		       const ScalorType & mass,
		       const ScalorType & charge);
    void writeData_SimpleFile (const char * filename);
// public:
//     void formDataTransferBlock (const IndexType & startIndex,
// 				const IndexType & num,
// 				DataTransferBlock & block);
public:
    HostCoordType * cptr_coordinate		() {return coord;}
    HostCoordNoiType * cptr_coordinateNoi	() {return coordNoi;}
    ScalorType * cptr_velocityX			() {return velox;}
    ScalorType * cptr_velocityY			() {return veloy;}
    ScalorType * cptr_velocityZ			() {return veloz;}
    ScalorType * cptr_forceX			() {return forcx;}
    ScalorType * cptr_forceY			() {return forcy;}
    ScalorType * cptr_forceZ			() {return forcz;}
    IndexType  * cptr_globalIndex		() {return globalIndex;}
    TypeType   * cptr_type			() {return type;}
    ScalorType * cptr_mass			() {return mass;}
    ScalorType * cptr_charge			() {return charge;}
    const HostCoordType * cptr_coordinate	() const {return coord;}
    const HostCoordNoiType * cptr_coordinateNoi	() const {return coordNoi;}
    const ScalorType * cptr_velocityX		() const {return velox;}
    const ScalorType * cptr_velocityY		() const {return veloy;}
    const ScalorType * cptr_velocityZ		() const {return veloz;}
    const ScalorType * cptr_forceX		() const {return forcx;}
    const ScalorType * cptr_forceY		() const {return forcy;}
    const ScalorType * cptr_forceZ		() const {return forcz;}
    const IndexType  * cptr_globalIndex		() const {return globalIndex;}
    const TypeType   * cptr_type		() const {return type;}
    const ScalorType * cptr_mass		() const {return mass;}
    const ScalorType * cptr_charge		() const {return charge;}
  };


  class GlobalHostMDData : public HostMDData
  {
    IndexType findMolIndex (const Topology::System & sysTop,
			    const IndexType & globalIndex);
    XDRFILE *xdfile;
    matrix xdbox;
    rvec *xdx;
    float xdprec;
public:
    IndexType numAtomInGroFile (const char * filename);
    void initConf_GroFile (const char * filename,
			   char * atomName, IndexType * atomIndex,
			   char * resdName, IndexType * resdIndex);
    void initTopology (const Topology::System & sysTop);
    void writeData_GroFile (const char * filename,
			    const char * atomName, const IndexType * atomIndex,
			    const char * resdName, const IndexType * resdIndex);
    void initWriteData_xtcFile (const char * filename, float prec=1000.f);
    void writeData_xtcFile (int step, float time);
    void endWriteData_xtcFile ();
  };


  void distributeGlobalMDData (const GlobalHostMDData & gdata,
			       HostMDData & ldata);

  
#ifdef DEVICE_CODE
  class DeviceMDData 
  {
public:
    IndexType numData_;
    IndexType memSize_;
    CoordType * coord;
    CoordNoiType * coordNoi;
    ScalorType * velox;
    ScalorType * veloy;
    ScalorType * veloz;
    ScalorType * forcx;
    ScalorType * forcy;
    ScalorType * forcz;
    IndexType  * globalIndex;
    RectangularBox globalBox;
    // RectangularBox localBox;
protected:
    // topology related
    TypeType   * type;
    ScalorType * mass;
    ScalorType * charge;
protected:
    bool malloced;
public:
    DeviceMDData ();
    DeviceMDData (const DeviceMDData & ddata);
    ~DeviceMDData ();
public:
    const IndexType & numData () const {return numData_;}
    const IndexType & memSize () const {return memSize_;}
    IndexType & numData ()  {return numData_;}

    const RectangularBox & getGlobalBox     () const {return globalBox;}
    const HostVectorType & getGlobalBoxSize () const {return globalBox.size;}
    void setGlobalBox (const RectangularBox & box)
	{ globalBox = box; }
    void setGlobalBox (const ScalorType & bx,
		       const ScalorType & by,
		       const ScalorType & bz)
	{ setBoxSize (bx, by, bz, &globalBox); }
    void easyMalloc (const IndexType &memSize);
    void initZero ();
    void clear ();
    void clearData () {numData_ = 0;}
    void copyFromHost (const HostMDData & hdata,
		       const MDDataItemMask_t mask = MDDataItemMask_All);
    void copyToHost   (HostMDData & hdata,
		       const MDDataItemMask_t mask = MDDataItemMask_All) const;
    void copyFromDevice (const DeviceMDData & ddata,
			 const MDDataItemMask_t mask = MDDataItemMask_All);
public:
    CoordType * dptr_coordinate			() {return coord;}
    CoordNoiType * dptr_coordinateNoi		() {return coordNoi;}
    ScalorType * dptr_velocityX			() {return velox;}
    ScalorType * dptr_velocityY			() {return veloy;}
    ScalorType * dptr_velocityZ			() {return veloz;}
    ScalorType * dptr_forceX			() {return forcx;}
    ScalorType * dptr_forceY			() {return forcy;}
    ScalorType * dptr_forceZ			() {return forcz;}
    IndexType  * dptr_globalIndex		() {return globalIndex;}
    TypeType   * dptr_type			() {return type;}
    ScalorType * dptr_mass			() {return mass;}
    ScalorType * dptr_charge			() {return charge;}
    const CoordType * dptr_coordinate		() const {return coord;}
    const CoordNoiType * dptr_coordinateNoi	() const {return coordNoi;}
    const ScalorType * dptr_velocityX		() const {return velox;}
    const ScalorType * dptr_velocityY		() const {return veloy;}
    const ScalorType * dptr_velocityZ		() const {return veloz;}
    const ScalorType * dptr_forceX		() const {return forcx;}
    const ScalorType * dptr_forceY		() const {return forcy;}
    const ScalorType * dptr_forceZ		() const {return forcz;}
    const IndexType  * dptr_globalIndex		() const {return globalIndex;}
    const TypeType   * dptr_type		() const {return type;}
    const ScalorType * dptr_mass		() const {return mass;}
    const ScalorType * dptr_charge		() const {return charge;}
  };
#endif
  // void cpyHostMDDataToDevice (const HostMDData & hdata,
  // 			      DeviceMDData

}

#endif

