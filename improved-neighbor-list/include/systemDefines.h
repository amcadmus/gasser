#ifndef __systemDefines_h_wanghan__
#define __systemDefines_h_wanghan__

#include <limits.h>

typedef float		ScalorType;
typedef int		IntScalorType;
typedef unsigned	IndexType;
// typedef unsigned	ForceIndexType;
typedef int		TypeType;

// #ifdef HOST_CODE
// __align(64) struct float4
// {
//   float x, y, z, w;
// };
// struct float3
// {
//   float x, y, z;
// };
// #endif // hostcode

#define MaxThreadsPerBlock			512
#define NThreadForSum				64
#define DefaultNThreadPerBlock			64
#define MinThreadsPerBlock			16
#define MaxIndexValue				UINT_MAX
#define MaxForceIndexValue			UINT_MAX
#define NUintBit				32
#define StringSize				8
#define MaxExceptionMsgLength			512
#define MaxNBForceTableBuffSize			512
#define MaxNumberParamPerForce			16	
#define MaxNumberNonBondedInteraction		32
#define MaxNumberNonBondedInteractionParameter	(32*5)
#define MaxNumberBondedInteraction		32
#define MaxNumberBondedInteractionParamemter	(32*3)
#define MaxLengthNonBondedInteractionTable	512

#define MemAllocExtension			2

#define SystemSharedBuffSize		16384
#define GlobalFunctionParamSizeLimit	256

#define M_PIF				3.14159265f

#ifdef DEVICE_CODE
typedef float4		CoordType;
typedef int3		CoordNoiType;
typedef float3		VectorType;
typedef int3		IntVectorType;
typedef uint3		IndexVectorType;
#endif


struct floatH4 
{
  float x, y, z, w;
  floatH4 () : x(0.f), y(0.f), z(0.f), w(0.f) {}
  // floatH4 & operator = (const floatH4 & a) 
  //     {x = a.x; y = a.y; z = a.z; w = a.w;
};
struct floatH3
{
  float x, y, z;
  floatH3 () : x(0.f), y(0.f), z(0.f) {}
};
struct intH3
{
  int x, y, z;
  intH3 () : x(0), y(0), z(0) {}
};

typedef floatH4		HostCoordType;
typedef intH3		HostCoordNoiType;
typedef floatH3		HostVectorType;
typedef intH3		HostIntVectorType;
// typedef floatH4		CoordType;

#endif
