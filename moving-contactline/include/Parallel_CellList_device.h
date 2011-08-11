#ifndef __Parallel_CellList_device_wanghan__
#define __Parallel_CellList_device_wanghan__

#define DEVICE_CODE

namespace Parallel{
  namespace CudaGlobal {
    __global__ void
    initZeroCell (const IntVectorType numCell,
		  IndexType * numAtomInCell);
    
    __global__ void
    formCellStructure (const VectorType frameLow,
		       const VectorType frameUp,
		       const IntVectorType numCell,
		       IndexType * numAtomInCell,
		       const IndexType numAtom,
		       const CoordType  * bk_coord,
		       const CoordNoiType * bk_coordNoix,
		       const ScalorType * bk_velox,
		       const ScalorType * bk_veloy,
		       const ScalorType * bk_veloz,
		       const IndexType  * bk_globalIndex,
		       const TypeType   * bk_type,
		       const ScalorType * bk_mass,
		       const ScalorType * bk_charge,
		       CoordType  * coord,
		       CoordNoiType * coordNoix,
		       ScalorType * velox,
		       ScalorType * veloy,
		       ScalorType * veloz,
		       IndexType  * globalIndex,
		       TypeType   * type,
		       ScalorType * mass,
		       ScalorType * charge,
		       mdError_t * ptr_de);
    __global__ void
    reinitCellStructure_calForwardMap (const IndexType * bk_numAtomInCell,
				       const CoordType * bk_coord,
				       const VectorType frameLow,
				       const VectorType frameUp,
				       const IntVectorType numCell,
				       IndexType * numAtomInCell,
				       IndexType * forwardMap,
				       mdError_t * ptr_de);
    __global__ void
    reinitCellStructure_step1 (const IndexType * bk_numAtomInCell,
			       const CoordType  * bk_coord,
			       const CoordNoiType * bk_coordNoi,
			       const ScalorType * bk_velox,
			       const ScalorType * bk_veloy,
			       const ScalorType * bk_veloz,
			       const ScalorType * bk_forcx,
			       const ScalorType * bk_forcy,
			       const ScalorType * bk_forcz,
			       const IndexType  * bk_globalIndex,
			       const TypeType   * bk_type,
			       const ScalorType * bk_mass,
			       const ScalorType * bk_charge,
			       const IndexType * forwardMap,
			       CoordType  * coord,
			       CoordNoiType * coordNoi,
			       ScalorType * velox,
			       ScalorType * veloy,
			       ScalorType * veloz,
			       ScalorType * forcx,
			       ScalorType * forcy,
			       ScalorType * forcz,
			       IndexType  * globalIndex,
			       TypeType   * type,
			       ScalorType * mass,
			       ScalorType * charge);
    __global__ void
    reinitCellStructure_step2 (const IndexType * bk_numAtomInCell,
			       const IndexType bk_maxNumBond,
			       const IndexType * bk_numBond,
			       const IndexType * bk_bondIndex,
			       const IndexType * bk_bondNeighbor_globalIndex,
			       const IndexType bk_maxNumAngle,
			       const IndexType * bk_numAngle,
			       const IndexType * bk_angleIndex,
			       const IndexType * bk_anglePosi,			   
			       const IndexType * bk_angleNeighbor_globalIndex,
			       const IndexType bk_maxNumDihedral,
			       const IndexType * bk_numDihedral,
			       const IndexType * bk_dihedralIndex,
			       const IndexType * bk_dihedralPosi,
			       const IndexType * bk_dihedralNeighbor_globalIndex,
			       const IndexType bk_bondTopStride,
			       const IndexType * forwardMap,
			       IndexType * numBond,
			       IndexType * bondIndex,
			       IndexType * bondNeighbor_globalIndex,
			       IndexType * numAngle,
			       IndexType * angleIndex,
			       IndexType * anglePosi,
			       IndexType * angleNeighbor_globalIndex,
			       IndexType * numDihedral,
			       IndexType * dihedralIndex,
			       IndexType * dihedralPosi,
			       IndexType * dihedralNeighbor_globalIndex,
			       const IndexType bondTopStride);

    __global__ void
    rebuildCellList_step1 (const VectorType frameLow,
			   const VectorType frameUp,
			   const IntVectorType numCell,
			   const IndexType * bk_numAtomInCell,
			   IndexType * numAtomInCell,
			   CoordType  * coord,
			   CoordNoiType * coordNoi,
			   ScalorType * velox,
			   ScalorType * veloy,
			   ScalorType * veloz,
			   ScalorType * forcx,
			   ScalorType * forcy,
			   ScalorType * forcz,
			   IndexType  * globalIndex,
			   TypeType   * type,
			   ScalorType * mass,
			   ScalorType * charge,
			   mdError_t * ptr_de);			 
    __global__ void
    rebuildCellList_step2 (IndexType * numAtomInCell,
			   CoordType  * coord,
			   CoordNoiType * coordNoi,
			   ScalorType * velox,
			   ScalorType * veloy,
			   ScalorType * veloz,
			   ScalorType * forcx,
			   ScalorType * forcy,
			   ScalorType * forcz,
			   IndexType  * globalIndex,
			   TypeType   * type,
			   ScalorType * mass,
			   ScalorType * charge,
			   mdError_t * ptr_de);
    __global__ void 
    rebuildCellList_step1 (const VectorType frameLow,
			   const VectorType frameUp,
			   const IntVectorType numCell,
			   const IndexType * bk_numAtomInCell,
			   IndexType * numAtomInCell,
			   CoordType * coord,
			   CoordNoiType * coordNoi,
			   ScalorType * velox,
			   ScalorType * veloy,
			   ScalorType * veloz,
			   ScalorType * forcx,
			   ScalorType * forcy,
			   ScalorType * forcz,
			   IndexType  * globalIndex,
			   TypeType   * type,
			   ScalorType * mass,
			   ScalorType * charge,
			   IndexType * forwardMap,
			   mdError_t * ptr_de,
			   IndexType * erridx,
			   ScalorType * errsrc);
    __global__ void 
    rebuildCellList_step2 (IndexType * numAtomInCell,
			   CoordType  * coord,
			   CoordNoiType * coordNoi,
			   ScalorType * velox,
			   ScalorType * veloy,
			   ScalorType * veloz,
			   ScalorType * forcx,
			   ScalorType * forcy,
			   ScalorType * forcz,
			   IndexType  * globalIndex,
			   TypeType   * type,
			   ScalorType * mass,
			   ScalorType * charge,
			   IndexType * forwardMap,
			   mdError_t * ptr_de);
    // __global__ void 
    // rebuildCellList_step1 (const VectorType frameLow,
    // 			   const VectorType frameUp,
    // 			   const IntVectorType numCell,
    // 			   const IndexType * bk_numAtomInCell,
    // 			   IndexType * numAtomInCell,
    // 			   CoordType * coord,
    // 			   CoordNoiType * coordNoi,
    // 			   ScalorType * velox,
    // 			   ScalorType * veloy,
    // 			   ScalorType * veloz,
    // 			   ScalorType * forcx,
    // 			   ScalorType * forcy,
    // 			   ScalorType * forcz,
    // 			   IndexType  * globalIndex,
    // 			   TypeType   * type,
    // 			   ScalorType * mass,
    // 			   ScalorType * charge,
    // 			   IndexType * global_numBond,
    // 			   IndexType * global_neighborIndex,
    // 			   IndexType * global_bondIndex,
    // 			   IndexType * neighborIndex,
    // 			   IndexType bondListStride,
    // 			   mdError_t * ptr_de);
    // __global__ void 
    // rebuildCellList_step2 (IndexType * numAtomInCell,
    // 			   CoordType  * coord,
    // 			   CoordNoiType * coordNoi,
    // 			   ScalorType * velox,
    // 			   ScalorType * veloy,
    // 			   ScalorType * veloz,
    // 			   ScalorType * forcx,
    // 			   ScalorType * forcy,
    // 			   ScalorType * forcz,
    // 			   IndexType  * globalIndex,
    // 			   TypeType   * type,
    // 			   ScalorType * mass,
    // 			   ScalorType * charge,
    // 			   IndexType * global_numBond,
    // 			   IndexType * global_neighborIndex,
    // 			   IndexType * global_bondIndex,
    // 			   IndexType * neighborIndex,
    // 			   IndexType bondListStride,
    // 			   IndexType maxNumBond,
    // 			   mdError_t * ptr_de);
    __global__ void 
    rebuildCellList_step1_mapBondTop (const IndexType * forwardMap,
				      const IndexType * bk_numAtomInCell,
				      IndexType maxNumBond,
				      IndexType * numBond,
				      IndexType * bondIndex,
				      IndexType * bondNeighbor_globalIndex,
				      IndexType * bondNeighbor_localIndex,
				      IndexType maxNumAngle,
				      IndexType * numAngle,
				      IndexType * angleIndex,
				      IndexType * anglePosi,
				      IndexType * angleNeighbor_globalIndex,
				      IndexType * angleNeighbor_localIndex,
				      IndexType maxNumDihedral,
				      IndexType * numDihedral,
				      IndexType * dihedralIndex,
				      IndexType * dihedralPosi,
				      IndexType * dihedralNeighbor_globalIndex,
				      IndexType * dihedralNeighbor_localIndex,
				      IndexType bondTopStride);
    __global__ void 
    rebuildCellList_step2_mapBondTop (const IndexType * forwardMap,
				      const IndexType * bk_numAtomInCell,
				      IndexType maxNumBond,
				      IndexType * numBond,
				      IndexType * bondIndex,
				      IndexType * bondNeighbor_globalIndex,
				      IndexType * bondNeighbor_localIndex,
				      IndexType maxNumAngle,
				      IndexType * numAngle,
				      IndexType * angleIndex,
				      IndexType * anglePosi,
				      IndexType * angleNeighbor_globalIndex,
				      IndexType * angleNeighbor_localIndex,
				      IndexType maxNumDihedral,
				      IndexType * numDihedral,
				      IndexType * dihedralIndex,
				      IndexType * dihedralPosi,
				      IndexType * dihedralNeighbor_globalIndex,
				      IndexType * dihedralNeighbor_localIndex,
				      IndexType bondTopStride);
    __global__ void 
    rebuildCellList_step1_mapBondTop (const IndexType * forwardMap,
				      const IndexType * bk_numAtomInCell,
				      IndexType maxNumBond,
				      IndexType * numBond,
				      IndexType * bondIndex,
				      IndexType * bondNeighbor_globalIndex,
				      IndexType maxNumAngle,
				      IndexType * numAngle,
				      IndexType * angleIndex,
				      IndexType * anglePosi,
				      IndexType * angleNeighbor_globalIndex,
				      IndexType maxNumDihedral,
				      IndexType * numDihedral,
				      IndexType * dihedralIndex,
				      IndexType * dihedralPosi,
				      IndexType * dihedralNeighbor_globalIndex,
				      IndexType bondTopStride);
    __global__ void 
    rebuildCellList_step2_mapBondTop (const IndexType * forwardMap,
				      const IndexType * bk_numAtomInCell,
				      IndexType maxNumBond,
				      IndexType * numBond,
				      IndexType * bondIndex,
				      IndexType * bondNeighbor_globalIndex,
				      IndexType maxNumAngle,
				      IndexType * numAngle,
				      IndexType * angleIndex,
				      IndexType * anglePosi,
				      IndexType * angleNeighbor_globalIndex,
				      IndexType maxNumDihedral,
				      IndexType * numDihedral,
				      IndexType * dihedralIndex,
				      IndexType * dihedralPosi,
				      IndexType * dihedralNeighbor_globalIndex,
				      IndexType bondTopStride);
    __global__ void
    packDeviceMDData (const IndexType * cellIndex,
		      const IndexType * numAtomInCell,
		      const IndexType * cellStartIndex,
		      const MDDataItemMask_t mask,
		      const CoordType  * source_coord,
		      const CoordNoiType * source_coordNoi,
		      const ScalorType * source_velox,
		      const ScalorType * source_veloy,
		      const ScalorType * source_veloz,
		      const ScalorType * source_forcx,
		      const ScalorType * source_forcy,
		      const ScalorType * source_forcz,
		      const IndexType  * source_globalIndex,
		      const TypeType   * source_type,
		      const ScalorType * source_mass,
		      const ScalorType * source_charge,
		      CoordType  * coord,
		      CoordNoiType * coordNoi,
		      ScalorType * velox,
		      ScalorType * veloy,
		      ScalorType * veloz,
		      ScalorType * forcx,
		      ScalorType * forcy,
		      ScalorType * forcz,
		      IndexType  * globalIndex,
		      TypeType   * type,
		      ScalorType * mass,
		      ScalorType * charge); 
    __global__ void 
    packDeviceMDData_bondTop (const IndexType * cellIndex,
			      const IndexType * numAtomInCell,
			      const IndexType * cellStartIndex,
			      const MDDataItemMask_t mask,
			      const IndexType * source_numBond,
			      const IndexType * source_bondIndex,
			      const IndexType * source_bondNeighbor_globalIndex,
			      const IndexType * source_numAngle,
			      const IndexType * source_angleIndex,
			      const IndexType * source_anglePosi ,
			      const IndexType * source_angleNeighbor_globalIndex,
			      const IndexType * source_numDihedral,
			      const IndexType * source_dihedralIndex,
			      const IndexType * source_dihedralPosi ,
			      const IndexType * source_dihedralNeighbor_globalIndex,
			      const IndexType   source_bondTopStride,
			      const IndexType   maxNumBond,
			      const IndexType   maxNumAngle,
			      const IndexType   maxNumDihedral,
			      IndexType * numBond,
			      IndexType * bondIndex,
			      IndexType * bondNeighbor_globalIndex,
			      IndexType * numAngle,
			      IndexType * angleIndex,
			      IndexType * anglePosi ,
			      IndexType * angleNeighbor_globalIndex,
			      IndexType * numDihedral,
			      IndexType * dihedralIndex,
			      IndexType * dihedralPosi ,
			      IndexType * dihedralNeighbor_globalIndex,
			      const IndexType bondTopStiide);
    
    __global__ void
    unpackDeviceMDData_replace (const IndexType * cellIndex,
				const IndexType * cellStartIndex,
				const MDDataItemMask_t mask,
				const CoordType  * source_coord,
				const CoordNoiType * source_coordNoi,
				const ScalorType * source_velox,
				const ScalorType * source_veloy,
				const ScalorType * source_veloz,
				const ScalorType * source_forcx,
				const ScalorType * source_forcy,
				const ScalorType * source_forcz,
				const IndexType  * source_globalIndex,
				const TypeType   * source_type,
				const ScalorType * source_mass,
				const ScalorType * source_charge,
				IndexType * numAtomInCell,
				CoordType  * coord,
				CoordNoiType * coordNoi,
				ScalorType * velox,
				ScalorType * veloy,
				ScalorType * veloz,
				ScalorType * forcx,
				ScalorType * forcy,
				ScalorType * forcz,
				IndexType  * globalIndex,
				TypeType   * type,
				ScalorType * mass,
				ScalorType * charge);
    
    __global__ void 
    unpackDeviceMDData_bondTop_replace (const IndexType * cellIndex,
					const IndexType * cellStartIndex,
					const MDDataItemMask_t mask,
					const IndexType * source_numBond,
					const IndexType * source_bondIndex,
					const IndexType * source_bondNeighbor_globalIndex,
					const IndexType * source_numAngle,
					const IndexType * source_angleIndex,
					const IndexType * source_anglePosi ,
					const IndexType * source_angleNeighbor_globalIndex,
					const IndexType * source_numDihedral,
					const IndexType * source_dihedralIndex,
					const IndexType * source_dihedralPosi ,
					const IndexType * source_dihedralNeighbor_globalIndex,
					const IndexType   source_bondTopStride,
					const IndexType   maxNumBond,
					const IndexType   maxNumAngle,
					const IndexType   maxNumDihedral,
					IndexType * numBond,
					IndexType * bondIndex,
					IndexType * bondNeighbor_globalIndex,
					IndexType * numAngle,
					IndexType * angleIndex,
					IndexType * anglePosi ,
					IndexType * angleNeighbor_globalIndex,
					IndexType * numDihedral,
					IndexType * dihedralIndex,
					IndexType * dihedralPosi ,
					IndexType * dihedralNeighbor_globalIndex,
					const IndexType bondTopStiide);
    __global__ void
    unpackDeviceMDData_add (const IndexType * cellIndex,
			    const IndexType * cellStartIndex,
			    const MDDataItemMask_t mask,
			    const CoordType  * source_coord,
			    const CoordNoiType * source_coordNoi,
			    const ScalorType * source_velox,
			    const ScalorType * source_veloy,
			    const ScalorType * source_veloz,
			    const ScalorType * source_forcx,
			    const ScalorType * source_forcy,
			    const ScalorType * source_forcz,
			    const IndexType  * source_globalIndex,
			    const TypeType   * source_type,
			    const ScalorType * source_mass,
			    const ScalorType * source_charge,
			    IndexType * numAtomInCell,
			    CoordType  * coord,
			    CoordNoiType * coordNoi,
			    ScalorType * velox,
			    ScalorType * veloy,
			    ScalorType * veloz,
			    ScalorType * forcx,
			    ScalorType * forcy,
			    ScalorType * forcz,
			    IndexType  * globalIndex,
			    TypeType   * type,
			    ScalorType * mass,
			    ScalorType * charge,
			    mdError_t * ptr_de);
    __global__ void 
    unpackDeviceMDData_bondTop_add (const IndexType * cellIndex,
				    const IndexType * cellStartIndex,
				    const IndexType * oldNumAtomInCell,
				    const MDDataItemMask_t mask,
				    const IndexType * source_numBond,
				    const IndexType * source_bondIndex,
				    const IndexType * source_bondNeighbor_globalIndex,
				    const IndexType * source_numAngle,
				    const IndexType * source_angleIndex,
				    const IndexType * source_anglePosi ,
				    const IndexType * source_angleNeighbor_globalIndex,
				    const IndexType * source_numDihedral,
				    const IndexType * source_dihedralIndex,
				    const IndexType * source_dihedralPosi ,
				    const IndexType * source_dihedralNeighbor_globalIndex,
				    const IndexType   source_bondTopStride,
				    const IndexType   maxNumBond,
				    const IndexType   maxNumAngle,
				    const IndexType   maxNumDihedral,
				    IndexType * numBond,
				    IndexType * bondIndex,
				    IndexType * bondNeighbor_globalIndex,
				    IndexType * numAngle,
				    IndexType * angleIndex,
				    IndexType * anglePosi ,
				    IndexType * angleNeighbor_globalIndex,
				    IndexType * numDihedral,
				    IndexType * dihedralIndex,
				    IndexType * dihedralPosi ,
				    IndexType * dihedralNeighbor_globalIndex,
				    const IndexType bondTopStiide,
				    mdError_t * ptr_de);
    __global__ void
    clearCellListData (const IndexType * deviceList,
		       IndexType num,
		       IndexType * numAtomInCell);
    __global__ void
    normalizeSystem_CellListed (RectangularBox box,
				const IndexType * numAtomInCell,
				CoordType * coord,
				CoordNoiType * coordNoi);
    __global__ void
    normalizeSystemOnGhost_CellListed (RectangularBox box,
				       const IndexType * numAtomInCell,
				       const IntVectorType numCell,
				       const IndexType divideLevel,
				       const int nx,
				       const int ny,
				       const int nz,
				       const int rankx,
				       const int ranky,
				       const int rankz,
				       CoordType * coord,
				       CoordNoiType * coordNoi);
    __global__ void
    buildCellNeighborhood (const IntVectorType numCell,
			   const IndexType devideLevel,
			   const ScalorType rlist,
			   const HostVectorType boxSize,
			   IndexType * numNeighbor,
			   IndexType * neighborCellIndex,
			   const IndexType stride);
    __global__ void
    buildCellNeighborhood (const IntVectorType numCell,
			   const IndexType devideLevel,
			   const ScalorType rlist,
			   const HostVectorType globalBoxSize,
			   const int rankx,
			   const int ranky,
			   const int rankz,
			   const int nProcDimx,
			   const int nProcDimy,
			   const int nProcDimz,
			   IndexType * numNeighbor,
			   IndexType * neighborCellIndex,
			   CoordNoiType * neighborShiftNoi,
			   const IndexType stride);
    __global__ void
    buildCellNeighborhood (const IntVectorType numCell,
			   const IndexType devideLevel,
			   const ScalorType rlist,
			   const HostVectorType globalBoxSize,
			   const int rankx,
			   const int ranky,
			   const int rankz,
			   const int nProcDimx,
			   const int nProcDimy,
			   const int nProcDimz,
			   const IndexType * subList0,
			   const IndexType length0,
			   const IndexType * subList1,
			   const IndexType length1,
			   IndexType * numNeighbor,
			   IndexType * neighborCellIndex,
			   CoordNoiType * neighborShiftNoi,
			   const IndexType stride);
    __global__ void 
    rescaleCoordinate (const IndexType * numAtomInCell,
		       const HostVectorType scale,
		       CoordType * coord);
    __global__ void 
    rescaleVelocity   (const IndexType * numAtomInCell,
		       const HostVectorType scale,
		       ScalorType * velox,
		       ScalorType * veloy,
		       ScalorType * veloz);
  }
  namespace CudaDevice{
    template <typename VEC, typename T>
    __device__ T
    D3toD1 (const VEC & NCell,
	    const T &ix,
	    const T &iy,
	    const T &iz)
    {
      return iz +
	  NCell.z * iy +
	  NCell.z * NCell.y * ix;
    }
    template <typename VEC, typename T>
    __device__ void
    D1toD3 (const VEC & NCell,
	    const T &i, 
	    T &x,
	    T &y,
	    T &z)
    {
      T tmp = i;
      z = tmp % (NCell.z);
      tmp = (tmp - z) / NCell.z;
      y = tmp % (NCell.y);
      x = (tmp - y) / NCell.y;
    }
  }
}


#endif