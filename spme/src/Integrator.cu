#include "MDSystem_interface.h"
#include "common.h"
#include "Integrator.h"
#include "Auxiliary.h"
#include <stdio.h>
#include "Statistic_interface.h"


__global__ void leapFrog1Step (const IndexType numAtom,
			       const ScalorType * massi,
			       CoordType * coord,
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
  coord[ii].x += dt * vx;
  vy = (veloy[ii] += dt * forcy[ii] * mi);
  coord[ii].y += dt * vy;
  vz = (veloz[ii] += dt * forcz[ii] * mi);
  coord[ii].z += dt * vz;
}

__global__ void leapFrog1Step (const IndexType numAtom,
			       const ScalorType * mass,
			       const ScalorType * massi,
			       CoordType * coord,
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
    coord[ii].x += dt * vx;
    vy = (veloy[ii] += dt * forcy[ii] * mi);
    coord[ii].y += dt * vy;
    vz = (veloz[ii] += dt * forcz[ii] * mi);
    coord[ii].z += dt * vz;
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

#ifndef COORD_IN_ONE_VEC
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
#else
__global__ void leapFrogStepX (const IndexType numAtom,
			       const ScalorType * massi,
			       CoordType * coord,
			       const ScalorType * velox,
			       const ScalorType * veloy, 
			       const ScalorType * veloz,
			       const ScalorType dt)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  if (ii >= numAtom) return;
  coord[ii].x += dt * velox[ii];
  coord[ii].y += dt * veloy[ii];
  coord[ii].z += dt * veloz[ii];
}
#endif


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
						   ScalorType * buffz)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = threadIdx.x + bid * blockDim.x;

  extern __shared__ volatile ScalorType buff[];
  buff[tid] = 0.f;
  buff[tid+blockDim.x] = 0.f;
  __syncthreads();
  if (ii < numAtom){
    buff[tid] = mass[ii] * velox[ii];
  }
  __syncthreads();
  sumVectorBlockBuffer (buff, blockDim.x);
  if (tid == 0) buffx[bid] = buff[0];
  __syncthreads();
  if (ii < numAtom){
    buff[tid] = mass[ii] * veloy[ii];
  }
  __syncthreads();
  sumVectorBlockBuffer (buff, blockDim.x);
  if (tid == 0) buffy[bid] = buff[0];
  __syncthreads();
  if (ii < numAtom){
    buff[tid] = mass[ii] * veloz[ii];
  }
  __syncthreads();
  sumVectorBlockBuffer (buff, blockDim.x);
  if (tid == 0) buffz[bid] = buff[0];
}

__global__ void removeFreedom (IndexType numAtom,
			       ScalorType * velox, 
			       ScalorType * veloy,
			       ScalorType * veloz,
			       ScalorType totalMassi,
			       ScalorType * sums)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  if (ii < numAtom) {
    velox[ii] -= sums[0] * totalMassi;
    veloy[ii] -= sums[1] * totalMassi;
    veloz[ii] -= sums[2] * totalMassi;
  }
}


#ifndef COORD_IN_ONE_VEC
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
#else
__global__ void velocityVerlet_part1 (const IndexType numAtom,
				      const ScalorType * massi,
				      CoordType * coord,
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
  coord[ii].x += dt * velox[ii];
  veloy[ii]  += hdtmi * forcy[ii];
  coord[ii].y += dt * veloy[ii];
  veloz[ii]  += hdtmi * forcz[ii];
  coord[ii].z += dt * veloz[ii];
}
#endif


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
}



__global__ void rescaleData (const IndexType numAtom,
			     ScalorType * data,
			     ScalorType alpha)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  
  if (ii < numAtom) data[ii] *= alpha;
}



__global__ void
leapFrogStepV_VCouple (const IndexType numAtom,
		       const ScalorType * massi,
		       ScalorType * velox,
		       ScalorType * veloy, 
		       ScalorType * veloz,
		       const ScalorType * forcx,
		       const ScalorType * forcy, 
		       const ScalorType * forcz,
		       const ScalorType lambda0,
		       const ScalorType lambda1,
		       const ScalorType lambda2,
		       const ScalorType dt)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  
  if (ii < numAtom) {
    ScalorType mi = massi[ii];
    ScalorType hdt = 0.5f * dt;
    ScalorType tmp = 1.f - hdt * lambda0;
    ScalorType tmp1= 1.f + hdt * lambda0;
    velox[ii] = (tmp * velox[ii] + dt * mi * forcx[ii]) / tmp1;
    tmp = 1.f - hdt * lambda1;
    tmp1= 1.f + hdt * lambda1;
    veloy[ii] = (tmp * veloy[ii] + dt * mi * forcy[ii]) / tmp1;
    tmp = 1.f - hdt * lambda2;
    tmp1= 1.f + hdt * lambda2;
    veloz[ii] = (tmp * veloz[ii] + dt * mi * forcz[ii]) / tmp1;
  }
}

__global__ void
leapFrogStepV_VCouple (const IndexType numAtom,
		       const ScalorType * mass,
		       const ScalorType * massi,
		       ScalorType * velox,
		       ScalorType * veloy, 
		       ScalorType * veloz,
		       const ScalorType * forcx,
		       const ScalorType * forcy, 
		       const ScalorType * forcz,
		       const ScalorType lambda0,
		       const ScalorType lambda1,
		       const ScalorType lambda2,
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
    ScalorType hdt = 0.5f * dt;
    ScalorType tmp = 1.f - hdt * lambda0;
    ScalorType tmp1= 1.f + hdt * lambda0 ;
    vx = (velox[ii] = ((tmp * velox[ii] + dt * mi * forcx[ii]) / tmp1));
    tmp = 1.f - hdt * lambda1;
    tmp1= 1.f + hdt * lambda1;
    vy = (veloy[ii] = ((tmp * veloy[ii] + dt * mi * forcy[ii]) / tmp1));
    tmp = 1.f - hdt * lambda2;
    tmp1= 1.f + hdt * lambda2;
    vz = (veloz[ii] = ((tmp * veloz[ii] + dt * mi * forcz[ii]) / tmp1));
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


