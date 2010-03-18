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
		       const IntScalorType * bk_coordNoix,
		       const IntScalorType * bk_coordNoiy,
		       const IntScalorType * bk_coordNoiz,
		       const ScalorType * bk_velox,
		       const ScalorType * bk_veloy,
		       const ScalorType * bk_veloz,
		       const IndexType  * bk_globalIndex,
		       const TypeType   * bk_type,
		       const ScalorType * bk_mass,
		       const ScalorType * bk_charge,
		       CoordType  * coord,
		       IntScalorType * coordNoix,
		       IntScalorType * coordNoiy,
		       IntScalorType * coordNoiz,
		       ScalorType * velox,
		       ScalorType * veloy,
		       ScalorType * veloz,
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
			   CoordType  * coord,
			   IntScalorType * coordNoix,
			   IntScalorType * coordNoiy,
			   IntScalorType * coordNoiz,
			   ScalorType * velox,
			   ScalorType * veloy,
			   ScalorType * veloz,
			   IndexType  * globalIndex,
			   TypeType   * type,
			   ScalorType * mass,
			   ScalorType * charge,
			   mdError_t * ptr_de);			 
    __global__ void
    rebuildCellList_step2 (IndexType * numAtomInCell,
			   CoordType  * coord,
			   IntScalorType * coordNoix,
			   IntScalorType * coordNoiy,
			   IntScalorType * coordNoiz,
			   ScalorType * velox,
			   ScalorType * veloy,
			   ScalorType * veloz,
			   IndexType  * globalIndex,
			   TypeType   * type,
			   ScalorType * mass,
			   ScalorType * charge,
			   mdError_t * ptr_de);
    __global__ void
    packDeviceMDData (const IndexType * cellIndex,
		      const IndexType * numAtomInCell,
		      const IndexType * cellStartIndex,
		      const MDDataItemMask_t mask,
		      const CoordType  * source_coord,
		      const IntScalorType * source_coordNoix,
		      const IntScalorType * source_coordNoiy,
		      const IntScalorType * source_coordNoiz,
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
		      IntScalorType * coordNoix,
		      IntScalorType * coordNoiy,
		      IntScalorType * coordNoiz,
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
    
    void __global__
    unpackDeviceMDData_replace (const IndexType * cellIndex,
				const IndexType * cellStartIndex,
				const MDDataItemMask_t mask,
				const CoordType  * source_coord,
				const IntScalorType * source_coordNoix,
				const IntScalorType * source_coordNoiy,
				const IntScalorType * source_coordNoiz,
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
				IntScalorType * coordNoix,
				IntScalorType * coordNoiy,
				IntScalorType * coordNoiz,
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
    
    void __global__
    unpackDeviceMDData_add (const IndexType * cellIndex,
			    const IndexType * cellStartIndex,
			    const MDDataItemMask_t mask,
			    const CoordType  * source_coord,
			    const IntScalorType * source_coordNoix,
			    const IntScalorType * source_coordNoiy,
			    const IntScalorType * source_coordNoiz,
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
			    IntScalorType * coordNoix,
			    IntScalorType * coordNoiy,
			    IntScalorType * coordNoiz,
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
    void __global__
    clearCellListData (const IndexType * deviceList,
		       IndexType num,
		       IndexType * numAtomInCell);
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