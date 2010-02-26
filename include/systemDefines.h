#ifndef __systemDefines_h_wanghan__
#define __systemDefines_h_wanghan__

//#include <limit.h>

typedef float		ScalorType;
typedef int		IntScalorType;
typedef unsigned	IndexType;
typedef unsigned	ForceIndexType;
typedef int		TypeType;

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
#define MaxNumberNonBondedInteraction		512
#define MaxNumberNonBondedInteractionParameter	(512*5)
#define MaxNumberBondedInteraction		256
#define MaxNumberBondedInteractionParamemter	(256*3)
#define MaxLengthNonBondedInteractionTable	128

#define SystemSharedBuffSize		16384
#define GlobalFunctionParamSizeLimit	256

#define M_PIF				3.14159265f

#ifndef CPP_FILE
typedef float4		CoordType;
typedef float3		VectorType;
typedef int3		IntVectorType;
typedef uint3		IndexVectorType;
#endif

#endif
