#include "Integrator.h"
#include "Auxiliary.h"
#include <stdio.h>
#include "Statistic_interface.h"

__global__ void leapFrog1Step (const IndexType numAtom,
			       const ScalorType * massi,
			       ScalorType * coordx,
			       ScalorType * coordy, 
			       ScalorType * coordz,
			       ScalorType * velox,
			       ScalorType * veloy, 
			       ScalorType * veloz,
			       const ScalorType * forcx,
			       const ScalorType * forcy, 
			       const ScalorType * forcz,
			       const ScalorType dt)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  if (ii >= numAtom) return;
  ScalorType mi = massi[ii];
  ScalorType vx, vy, vz;
  vx = (velox[ii] += dt * forcx[ii] * mi);
  coordx[ii] += dt * vx;
  vy = (veloy[ii] += dt * forcy[ii] * mi);
  coordy[ii] += dt * vy;
  vz = (veloz[ii] += dt * forcz[ii] * mi);
  coordz[ii] += dt * vz;
}

__global__ void leapFrog1Step (const IndexType numAtom,
			       const ScalorType * mass,
			       const ScalorType * massi,
			       ScalorType * coordx,
			       ScalorType * coordy, 
			       ScalorType * coordz,
			       ScalorType * velox,
			       ScalorType * veloy, 
			       ScalorType * veloz,
			       const ScalorType * forcx,
			       const ScalorType * forcy, 
			       const ScalorType * forcz,
			       const ScalorType dt,
			       ScalorType * statistic_buffxx,
			       ScalorType * statistic_buffyy,
			       ScalorType * statistic_buffzz)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  ScalorType vx, vy, vz;
  
  if (ii < numAtom) {
    ScalorType mi = massi[ii];
    vx = (velox[ii] += dt * forcx[ii] * mi);
    coordx[ii] += dt * vx;
    vy = (veloy[ii] += dt * forcy[ii] * mi);
    coordy[ii] += dt * vy;
    vz = (veloz[ii] += dt * forcz[ii] * mi);
    coordz[ii] += dt * vz;
  }

  extern __shared__ volatile ScalorType buff [];

  ScalorType scalor;
  if (ii < numAtom) scalor = 0.5f * mass[ii];
  else scalor = 0.f;
  
  if (ii < numAtom){
    buff[threadIdx.x] = scalor * vx * vx;
  }
  else {
    buff[threadIdx.x] = 0.f;
  }
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) statistic_buffxx[bid] = buff[0];
  __syncthreads();
  if (ii < numAtom){
    buff[threadIdx.x] = scalor * vy * vy;
  }
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) statistic_buffyy[bid] = buff[0];
  __syncthreads();
  if (ii < numAtom){
    buff[threadIdx.x] = scalor * vz * vz;
  }
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) statistic_buffzz[bid] = buff[0];
}

__global__ void leapFrogStepX (const IndexType numAtom,
			       const ScalorType * massi,
			       ScalorType * coordx,
			       ScalorType * coordy, 
			       ScalorType * coordz,
			       const ScalorType * velox,
			       const ScalorType * veloy, 
			       const ScalorType * veloz,
			       const ScalorType dt)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  if (ii >= numAtom) return;
  coordx[ii] += dt * velox[ii];
  coordy[ii] += dt * veloy[ii];
  coordz[ii] += dt * veloz[ii];
}

__global__ void leapFrogStepV (const IndexType numAtom,
			       const ScalorType * massi,
			       ScalorType * velox,
			       ScalorType * veloy, 
			       ScalorType * veloz,
			       const ScalorType * forcx,
			       const ScalorType * forcy, 
			       const ScalorType * forcz,
			       const ScalorType dt)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  
  if (ii < numAtom) {
    ScalorType mi = massi[ii];
    (velox[ii] += dt * forcx[ii] * mi);
    (veloy[ii] += dt * forcy[ii] * mi);
    (veloz[ii] += dt * forcz[ii] * mi);
  }
}

__global__ void leapFrogStepV (const IndexType numAtom,
			       const ScalorType * mass,
			       const ScalorType * massi,
			       ScalorType * velox,
			       ScalorType * veloy, 
			       ScalorType * veloz,
			       const ScalorType * forcx,
			       const ScalorType * forcy, 
			       const ScalorType * forcz,
			       const ScalorType dt,
			       ScalorType * statistic_buffxx,
			       ScalorType * statistic_buffyy,
			       ScalorType * statistic_buffzz)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  ScalorType vx, vy, vz;
  
  if (ii < numAtom) {
    ScalorType mi = massi[ii];
    vx = (velox[ii] += dt * forcx[ii] * mi);
    vy = (veloy[ii] += dt * forcy[ii] * mi);
    vz = (veloz[ii] += dt * forcz[ii] * mi);
  }

  extern __shared__ volatile ScalorType buff [];

  ScalorType scalor;
  if (ii < numAtom) scalor = 0.5f * mass[ii];
  else scalor = 0.f;
  
  if (ii < numAtom){
    buff[threadIdx.x] = scalor * vx * vx;
  }
  else {
    buff[threadIdx.x] = 0.f;
  }
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) statistic_buffxx[bid] = buff[0];
  __syncthreads();
  if (ii < numAtom){
    buff[threadIdx.x] = scalor * vy * vy;
  }
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) statistic_buffyy[bid] = buff[0];
  __syncthreads();
  if (ii < numAtom){
    buff[threadIdx.x] = scalor * vz * vz;
  }
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) statistic_buffzz[bid] = buff[0];
}




__device__ IndexType integrator_counter_prepare_x = 0;
__device__ IndexType integrator_counter_prepare_y = 0;
__device__ IndexType integrator_counter_prepare_z = 0;
__global__ void initRemoveTranslationalFreedom ()
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  if (tid + bid == 0){
    integrator_counter_prepare_x = 0;
    integrator_counter_prepare_y = 0;
    integrator_counter_prepare_z = 0;
  }
}

__global__ void prepareRemoveTranslationalFreedom (IndexType numAtom,
						   ScalorType * mass,
						   ScalorType * velox,
						   ScalorType * veloy,
						   ScalorType * veloz,
						   ScalorType * buffx,
						   ScalorType * buffy,
						   ScalorType * buffz,
						   ScalorType * sums)
{
  sumVectorMomentum(mass, velox, numAtom,
		    buffx, &integrator_counter_prepare_x,
		    &sums[0]);
  sumVectorMomentum(mass, veloy, numAtom,
		    buffy, &integrator_counter_prepare_y,
		    &sums[1]);
  sumVectorMomentum(mass, veloz, numAtom,
		    buffz, &integrator_counter_prepare_z,
		    &sums[2]);
}

__global__ void removeFreedom (IndexType numAtom,
			       ScalorType * velox, 
			       ScalorType * veloy,
			       ScalorType * veloz,
			       ScalorType totalMassi,
			       ScalorType * sums)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;
  if (ii >= numAtom) return;
  velox[ii] -= sums[0] * totalMassi;
  veloy[ii] -= sums[1] * totalMassi;
  veloz[ii] -= sums[2] * totalMassi;
}



__global__ void velocityVerlet_part1 (const IndexType numAtom,
				      const ScalorType * massi,
				      ScalorType * coordx,
				      ScalorType * coordy, 
				      ScalorType * coordz,
				      ScalorType * velox,
				      ScalorType * veloy, 
				      ScalorType * veloz,
				      const ScalorType * forcx,
				      const ScalorType * forcy, 
				      const ScalorType * forcz,
				      const ScalorType dt)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;
  if (ii >= numAtom) return;
  
  ScalorType hdtmi = 0.5f*dt*massi[ii];
  velox[ii]  += hdtmi * forcx[ii];
  coordx[ii] += dt * velox[ii];
  veloy[ii]  += hdtmi * forcy[ii];
  coordy[ii] += dt * veloy[ii];
  veloz[ii]  += hdtmi * forcz[ii];
  coordz[ii] += dt * veloz[ii];
}


__global__ void velocityVerlet_part2 (const IndexType numAtom,
				      const ScalorType * massi,
				      ScalorType * velox,
				      ScalorType * veloy, 
				      ScalorType * veloz,
				      const ScalorType * forcx,
				      const ScalorType * forcy, 
				      const ScalorType * forcz,
				      const ScalorType dt)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;
  if (ii >= numAtom) return;

  ScalorType hdtmi = 0.5f*dt*massi[ii];
  velox[ii] += hdtmi * forcx[ii];
  veloy[ii] += hdtmi * forcy[ii];
  veloz[ii] += hdtmi * forcz[ii];
}


__global__ void velocityVerlet_part2 (const IndexType numAtom,
				      const ScalorType * mass,
				      const ScalorType * massi,
				      ScalorType * velox,
				      ScalorType * veloy, 
				      ScalorType * veloz,
				      const ScalorType * forcx,
				      const ScalorType * forcy, 
				      const ScalorType * forcz,
				      const ScalorType dt,
				      ScalorType * statistic_buffxx,
				      ScalorType * statistic_buffyy,
				      ScalorType * statistic_buffzz)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  extern __shared__ volatile ScalorType buff [];
ScalorType vx(0.f), vy(0.f), vz(0.f);
  if (ii < numAtom) {
    ScalorType hdtmi = 0.5f*dt*massi[ii];
    vx = (velox[ii] += hdtmi * forcx[ii]);
    vy = (veloy[ii] += hdtmi * forcy[ii]);
    vz = (veloz[ii] += hdtmi * forcz[ii]);
  }

  ScalorType scalor;
  if (ii < numAtom) scalor = 0.5f * mass[ii];
  else scalor = 0.f;
  
  if (ii < numAtom){
    buff[threadIdx.x] = scalor * vx * vx;
  }
  else {
    buff[threadIdx.x] = 0.f;
  }
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) statistic_buffxx[bid] = buff[0];
  __syncthreads();
  if (ii < numAtom){
    buff[threadIdx.x] = scalor * vy * vy;
  }
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) statistic_buffyy[bid] = buff[0];
  __syncthreads();
  if (ii < numAtom){
    buff[threadIdx.x] = scalor * vz * vz;
  }
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) statistic_buffzz[bid] = buff[0];
  
  
  // __syncthreads();
  // IndexType num = ((bid+1) * blockDim.x > numAtom) ? 
  //     (blockDim.x - ((bid+1) * blockDim.x) + numAtom) : blockDim.x;
  // __shared__ volatile ScalorType buff [MaxThreadsPerBlock * 2];
  // buff[tid] = 0.0f;
  // buff[tid + blockDim.x] = 0.0f;
  // if (tid < num)
  //   buff[tid] = mass[ii] * (vx*vx + vy*vy + vz*vz);
  // __syncthreads();
  
  // __shared__ volatile  bool isLastBlockDone_for_velocityVerlet;
  // ScalorType partialSum = sumVectorBlockBuffer (buff, num);
  // if (tid == 0){
  //   statistic_buff[bid] = partialSum;
  //   IndexType value = atomicInc(&integrator_counter_for_velocityVerlet, gridDim.x*gridDim.y);
  //   isLastBlockDone_for_velocityVerlet = (value == (gridDim.x*gridDim.y - 1));
  //   // printf ("bid %d, sum %f, mark %d\n", bid, partialSum, isLastBlockDone_for_velocityVerlet);
  // }
  // __threadfence();
  // __syncthreads();
  
  // if (isLastBlockDone_for_velocityVerlet){
  //   IndexType p = 0;
  //   ScalorType tmpsum = 0.0f;
  //   while (p < gridDim.x*gridDim.y){
  //     IndexType tmp = gridDim.x*gridDim.y - p;
  //     IndexType n = ((tmp < blockDim.x) ? tmp : blockDim.x);
  //     // printf ("%d\n", blockDim.x);
  //     tmpsum += sumVectorBlock (statistic_buff, p, n);
  //     p += blockDim.x;
  //   }
  //   if (tid == 0){
  //     integrator_counter_for_velocityVerlet = 0;
  //     stddata[mdStatisticKineticEnergy] = tmpsum * 0.5f;
  //   }
  // }  
  // __threadfence();
}

__global__ void velocityVerlet_part2a (const IndexType numAtom,
				       const ScalorType * mass,
				       const ScalorType * massi,
				       ScalorType * velox,
				       ScalorType * veloy, 
				       ScalorType * veloz,
				       const ScalorType * forcx,
				       const ScalorType * forcy, 
				       const ScalorType * forcz,
				       const ScalorType dt,
				       ScalorType * statistic_buff)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;

  ScalorType vx, vy, vz;
  if (ii < numAtom) { 
    ScalorType hdtmi = 0.5f*dt*massi[ii];
    vx = (velox[ii] += hdtmi * forcx[ii]);
    vy = (veloy[ii] += hdtmi * forcy[ii]);
    vz = (veloz[ii] += hdtmi * forcz[ii]);
  }

  extern __shared__ volatile ScalorType buff [];
  if (ii < numAtom)
    buff[tid] = 0.5 * mass[ii] * (vx*vx + vy*vy + vz*vz);
  else
    buff[tid] = 0.f;
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) statistic_buff[bid] = buff[0];
  
  // __syncthreads();
  // IndexType num = ((bid+1) * blockDim.x > numAtom) ? 
  //     (blockDim.x - ((bid+1) * blockDim.x) + numAtom) : blockDim.x;
  // __shared__ volatile ScalorType buff [MaxThreadsPerBlock * 2];
  // buff[tid] = 0.0f;
  // buff[tid + blockDim.x] = 0.0f;
  // if (tid < num)
  //   buff[tid] = mass[ii] * (vx*vx + vy*vy + vz*vz);
  // __syncthreads();
  
  // __shared__ volatile bool isLastBlockDone_for_velocityVerlet;
  // ScalorType partialSum = sumVectorBlockBuffer (buff, num);
  // if (tid == 0){
  //   statistic_buff[bid] = partialSum;
  //   IndexType value = atomicInc(&integrator_counter_for_velocityVerlet, gridDim.x*gridDim.y);
  //   isLastBlockDone_for_velocityVerlet = (value == (gridDim.x*gridDim.y - 1));
  //   // printf ("bid %d, sum %f, mark %d\n", bid, partialSum, isLastBlockDone_for_velocityVerlet);
  // }
  // __threadfence();
  // __syncthreads();
  
  // if (isLastBlockDone_for_velocityVerlet){
  //   IndexType p = 0;
  //   ScalorType tmpsum = 0.0f;
  //   while (p < gridDim.x*gridDim.y){
  //     IndexType tmp = gridDim.x*gridDim.y - p;
  //     IndexType n = ((tmp < blockDim.x) ? tmp : blockDim.x);
  //     // printf ("%d\n", blockDim.x);
  //     tmpsum += sumVectorBlock (statistic_buff, p, n);
  //     p += blockDim.x;
  //   }
  //   if (tid == 0){
  //     integrator_counter_for_velocityVerlet = 0;
  //     *kineticE = tmpsum * 0.5f;
  //   }
  // }  
}





__global__ void velocityRescale_rescale (const IndexType numAtom,
					 ScalorType * velox,
					 ScalorType * veloy, 
					 ScalorType * veloz,
					 const ScalorType alpha)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;
  
  if (ii < numAtom) {
    (velox[ii] *= alpha);
    (veloy[ii] *= alpha);
    (veloz[ii] *= alpha);
  }
}

__global__ void velocityRescale_rescale (const IndexType numAtom,
					 const ScalorType * mass,
					 ScalorType * velox,
					 ScalorType * veloy, 
					 ScalorType * veloz,
					 const ScalorType alpha,
					 ScalorType * statistic_buffxx,
					 ScalorType * statistic_buffyy,
					 ScalorType * statistic_buffzz)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;
  ScalorType vx, vy, vz;
  
  if (ii < numAtom) {
    vx = (velox[ii] *= alpha);
    vy = (veloy[ii] *= alpha);
    vz = (veloz[ii] *= alpha);
  }

  extern __shared__ volatile ScalorType buff [];

  ScalorType scalor;
  if (ii < numAtom) scalor = 0.5 * mass[ii];
  else scalor = 0.f;
  
  if (ii < numAtom){
    buff[threadIdx.x] = scalor * vx * vx;
  }
  else {
    buff[threadIdx.x] = 0.f;
  }
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) statistic_buffxx[bid] = buff[0];
  __syncthreads();
  if (ii < numAtom){
    buff[threadIdx.x] = scalor * vy * vy;
  }
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) statistic_buffyy[bid] = buff[0];
  __syncthreads();
  if (ii < numAtom){
    buff[threadIdx.x] = scalor * vz * vz;
  }
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) statistic_buffzz[bid] = buff[0];

  // __syncthreads();
  // IndexType num = ((bid+1) * blockDim.x > numAtom) ? 
  //     (blockDim.x - ((bid+1) * blockDim.x) + numAtom) : blockDim.x;
  // __shared__ volatile ScalorType buff [MaxThreadsPerBlock * 2];
  // buff[tid] = 0.0f;
  // buff[tid + blockDim.x] = 0.0f;
  // if (tid < num)
  //   buff[tid] = mass[ii] * (vx*vx +vy*vy + vz*vz);
  // __syncthreads();
  
  // __shared__ volatile bool isLastBlockDone_for_velocityVerlet;
  // ScalorType partialSum = sumVectorBlockBuffer (buff, num);
  // if (tid == 0){
  //   statistic_buff[bid] = partialSum;
  //   IndexType value = atomicInc(&integrator_counter_for_velocityVerlet, gridDim.x*gridDim.y);
  //   isLastBlockDone_for_velocityVerlet = (value == (gridDim.x*gridDim.y - 1));
  //   // printf ("bid %d, sum %f, mark %d\n", bid, partialSum, isLastBlockDone_for_velocityVerlet);
  // }
  // __threadfence();
  // __syncthreads();
  
  // if (isLastBlockDone_for_velocityVerlet){
  //   IndexType p = 0;
  //   ScalorType tmpsum = 0;
  //   while (p < gridDim.x*gridDim.y){
  //     IndexType tmp = gridDim.x*gridDim.y - p;
  //     IndexType n = ((tmp < blockDim.x) ? tmp : blockDim.x);
  //     tmpsum += sumVectorBlock (statistic_buff, p, n);
  //     p += blockDim.x;
  //   }
  //   if (tid == 0){
  //     integrator_counter_for_velocityVerlet = 0;
  //     stddata[mdStatisticKineticEnergy] = tmpsum * 0.5f;
  //   }
  // }  
}



__global__ void rescaleData (const IndexType numAtom,
			     ScalorType * data,
			     ScalorType alpha)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  
  if (ii < numAtom) data[ii] *= alpha;
}

