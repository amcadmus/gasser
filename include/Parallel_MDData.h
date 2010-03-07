#ifndef __Parallel_MDData_h_wanghan__
#define __Parallel_MDData_h_wanghan__

#include "common.h"
#include "Topology.h"
#include "BoxGeometry.h"

using namespace RectangularBoxGeometry;

namespace Parallel{
  class MDExcptInconsistentMemorySizeOnHostMDData : public MDException {};
  class MDExcptNumAtomMoreThanMemSize : public MDException {};
  class MDExcptInvalidTopology : public MDException {};

  class HostMDData ;
  class GlobalHostMDData;
  class DeviceMDData;
  
  class HostMDData 
  {
    friend class DeviceMDData ;
protected:
    IndexType numAtom_;
    IndexType memSize_;
    HostCoordType * coord;
    IntScalorType * coordNoix;
    IntScalorType * coordNoiy;
    IntScalorType * coordNoiz;
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
    ~HostMDData ();
public:
    const IndexType & numAtom () const {return numAtom_;}
    const IndexType & memSize () const {return memSize_;}
    IndexType & numAtom ()  {return numAtom_;}
    // IndexType & memSize ()  {return memSize_;}
    void clearAll ();
    void clearData () {numAtom_ = 0;}
    void reallocAll (const IndexType & memSize);
public:
    const RectangularBox & getGlobalBox	() const {return globalBox;}
    void setGlobalBox (const ScalorType & bx,
		       const ScalorType & by,
		       const ScalorType & bz)
	{ setBoxSize (bx, by, bz, &globalBox); }
    void setGlobalBox (const RectangularBox & box) { globalBox = box; }
    void pushBackAtom (const HostCoordType & coord,
		       const IntScalorType & coordNoix,
		       const IntScalorType & coordNoiy,
		       const IntScalorType & coordNoiz,
		       const ScalorType & velox,
		       const ScalorType & veloy,
		       const ScalorType & veloz,
		       const IndexType & globalIndex,
		       const TypeType & type,
		       const ScalorType & mass,
		       const ScalorType & charge);
    void writeData_SimpleFile (const char * filename);
public:
    HostCoordType * cptr_coordinate		() {return coord;}
    IntScalorType * cptr_coordinateNoiX		() {return coordNoix;}
    IntScalorType * cptr_coordinateNoiY		() {return coordNoiy;}
    IntScalorType * cptr_coordinateNoiZ		() {return coordNoiz;}
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
    const IntScalorType * cptr_coordinateNoiX	() const {return coordNoix;}
    const IntScalorType * cptr_coordinateNoiY	() const {return coordNoiy;}
    const IntScalorType * cptr_coordinateNoiZ	() const {return coordNoiz;}
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
public:
    IndexType numAtomInGroFile (const char * filename);
    void initConf_GroFile (const char * filename,
			   char * atomName, IndexType * atomIndex,
			   char * resdName, IndexType * resdIndex);
    void initTopology (const Topology::System & sysTop);
  };


  void distributeGlobalMDData (const GlobalHostMDData & gdata,
			       HostMDData & ldata);

  
#ifdef DEVICE_CODE
  class DeviceMDData 
  {
protected:
    IndexType numAtom_;
    IndexType memSize_;
    CoordType * coord;
    IntScalorType * coordNoix;
    IntScalorType * coordNoiy;
    IntScalorType * coordNoiz;
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
    ~DeviceMDData ();
public:
    const IndexType & numAtom () const {return numAtom_;}
    const IndexType & memSize () const {return memSize_;}
    IndexType & numAtom ()  {return numAtom_;}

    void setGlobalBox (const RectangularBox & box)
	{ globalBox = box; }
    void setGlobalBox (const ScalorType & bx,
		       const ScalorType & by,
		       const ScalorType & bz)
	{ setBoxSize (bx, by, bz, &globalBox); }

    void mallocAll (const IndexType &memSize);
    void clearAll ();
    void clearData () {numAtom_ = 0;}
    void copyFromHost (const HostMDData & hdata);
    void copyToHost   (HostMDData & hdata) const;
public:
    CoordType * dptr_coordinate			() {return coord;}
    IntScalorType * dptr_coordinateNoiX		() {return coordNoix;}
    IntScalorType * dptr_coordinateNoiY		() {return coordNoiy;}
    IntScalorType * dptr_coordinateNoiZ		() {return coordNoiz;}
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
    const IntScalorType * dptr_coordinateNoiX	() const {return coordNoix;}
    const IntScalorType * dptr_coordinateNoiY	() const {return coordNoiy;}
    const IntScalorType * dptr_coordinateNoiZ	() const {return coordNoiz;}
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

