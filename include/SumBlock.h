#ifndef __SumBlock_h_wanghan__
#define __SumBlock_h_wanghan__

namespace SumBlock {
    template <typename TYPENAME>
    __device__ TYPENAME sum32_1bsize (TYPENAME * buff);
    template <typename TYPENAME>
    __device__ TYPENAME sum32_2bsize (TYPENAME * buff);
    template <typename TYPENAME>
    __device__ TYPENAME sum64_1bsize (TYPENAME * buff);
    template <typename TYPENAME>
    __device__ TYPENAME sum128_1bsize (TYPENAME * buff);
};


template <typename TYPENAME> __device__ TYPENAME
SumBlock::sum32_1bsize (TYPENAME * buff)
{
  __syncthreads();
  TYPENAME tmp;

  // // +1 cycle
  if (threadIdx.x < (32-1)){
    tmp = buff[threadIdx.x] + buff[threadIdx.x+1];
  }
  else {
    tmp = buff[threadIdx.x];
  }
  // __syncthreads();
  buff[threadIdx.x] = tmp;
  // __syncthreads();
  // // +2 cycle
  if (threadIdx.x < (32-2)){
    tmp = buff[threadIdx.x] + buff[threadIdx.x+2];
  }
  else {
    tmp = buff[threadIdx.x];
  }
  // __syncthreads();
  buff[threadIdx.x] = tmp;
  // __syncthreads();
  // // +4 cycle
  if (threadIdx.x < (32-4)){
    tmp = buff[threadIdx.x] + buff[threadIdx.x+4];
  }
  else {
    tmp = buff[threadIdx.x];
  }
  // __syncthreads();
  buff[threadIdx.x] = tmp;
  // __syncthreads();
  // // +8 cycle
  if (threadIdx.x < (32-8)){
    tmp = buff[threadIdx.x] + buff[threadIdx.x+8];
  }
  else {
    tmp = buff[threadIdx.x];
  }
  // __syncthreads();
  buff[threadIdx.x] = tmp;
  // __syncthreads();
  // // +16 cycle
  if (threadIdx.x < (32-16)){
    tmp = buff[threadIdx.x] + buff[threadIdx.x+16];
  }
  else {
    tmp = buff[threadIdx.x];
  }
  // __syncthreads();
  buff[threadIdx.x] = tmp;
  // __syncthreads();

  return buff[0];
}


template <typename TYPENAME> __device__ TYPENAME
SumBlock::sum64_1bsize (TYPENAME * buff)
{
  __syncthreads();
  TYPENAME tmp;

  // +1 cycle
  if (threadIdx.x < (64-1)){
    tmp = buff[threadIdx.x] + buff[threadIdx.x+1];
  }
  else {
    tmp = buff[threadIdx.x];
  }
  __syncthreads();
  buff[threadIdx.x] = tmp;
  __syncthreads();
  // +2 cycle
  if (threadIdx.x < (64-2)){
    tmp = buff[threadIdx.x] + buff[threadIdx.x+2];
  }
  else {
    tmp = buff[threadIdx.x];
  }
  __syncthreads();
  buff[threadIdx.x] = tmp;
  __syncthreads();
  // +4 cycle
  if (threadIdx.x < (64-4)){
    tmp = buff[threadIdx.x] + buff[threadIdx.x+4];
  }
  else {
    tmp = buff[threadIdx.x];
  }
  __syncthreads();
  buff[threadIdx.x] = tmp;
  __syncthreads();
  // +8 cycle
  if (threadIdx.x < (64-8)){
    tmp = buff[threadIdx.x] + buff[threadIdx.x+8];
  }
  else {
    tmp = buff[threadIdx.x];
  }
  __syncthreads();
  buff[threadIdx.x] = tmp;
  __syncthreads();
  // +16 cycle
  if (threadIdx.x < (64-16)){
    tmp = buff[threadIdx.x] + buff[threadIdx.x+16];
  }
  else {
    tmp = buff[threadIdx.x];
  }
  __syncthreads();
  buff[threadIdx.x] = tmp;
  __syncthreads();
  // +64 cycle
  if (threadIdx.x < (64-32)){
    tmp = buff[threadIdx.x] + buff[threadIdx.x+32];
  }
  else {
    tmp = buff[threadIdx.x];
  }
  __syncthreads();
  buff[threadIdx.x] = tmp;
  __syncthreads();

  return buff[0];
}


template <typename TYPENAME> __device__ TYPENAME
SumBlock::sum128_1bsize (TYPENAME * buff)
{
  __syncthreads();
  TYPENAME tmp;

  // +1 cycle
  if (threadIdx.x < (128-1)){
    tmp = buff[threadIdx.x] + buff[threadIdx.x+1];
  }
  else {
    tmp = buff[threadIdx.x];
  }
  __syncthreads();
  buff[threadIdx.x] = tmp;
  __syncthreads();
  // +2 cycle
  if (threadIdx.x < (128-2)){
    tmp = buff[threadIdx.x] + buff[threadIdx.x+2];
  }
  else {
    tmp = buff[threadIdx.x];
  }
  __syncthreads();
  buff[threadIdx.x] = tmp;
  __syncthreads();
  // +4 cycle
  if (threadIdx.x < (128-4)){
    tmp = buff[threadIdx.x] + buff[threadIdx.x+4];
  }
  else {
    tmp = buff[threadIdx.x];
  }
  __syncthreads();
  buff[threadIdx.x] = tmp;
  __syncthreads();
  // +8 cycle
  if (threadIdx.x < (128-8)){
    tmp = buff[threadIdx.x] + buff[threadIdx.x+8];
  }
  else {
    tmp = buff[threadIdx.x];
  }
  __syncthreads();
  buff[threadIdx.x] = tmp;
  __syncthreads();
  // +16 cycle
  if (threadIdx.x < (128-16)){
    tmp = buff[threadIdx.x] + buff[threadIdx.x+16];
  }
  else {
    tmp = buff[threadIdx.x];
  }
  __syncthreads();
  buff[threadIdx.x] = tmp;
  __syncthreads();
  // +32 cycle
  if (threadIdx.x < (128-32)){
    tmp = buff[threadIdx.x] + buff[threadIdx.x+32];
  }
  else {
    tmp = buff[threadIdx.x];
  }
  __syncthreads();
  buff[threadIdx.x] = tmp;
  __syncthreads();
  // +64 cycle
  if (threadIdx.x < (128-64)){
    tmp = buff[threadIdx.x] + buff[threadIdx.x+64];
  }
  else {
    tmp = buff[threadIdx.x];
  }
  __syncthreads();
  buff[threadIdx.x] = tmp;
  __syncthreads();

  return buff[0];
}




template <typename TYPENAME> __device__ TYPENAME
SumBlock::sum32_2bsize (TYPENAME * buff)
{
  buff[threadIdx.x + blockDim.x] = TYPENAME(0);
  __syncthreads();
  TYPENAME tmp;

  // // +1 cycle
  if (threadIdx.x < (32-1)){
    tmp = buff[threadIdx.x] + buff[threadIdx.x+1];
  }
  else {
    tmp = buff[threadIdx.x];
  }
  buff[threadIdx.x] = tmp;
  // // +2 cycle
  if (threadIdx.x < (32-2)){
    tmp = buff[threadIdx.x] + buff[threadIdx.x+2];
  }
  else {
    tmp = buff[threadIdx.x];
  }
  buff[threadIdx.x] = tmp;
  // // +4 cycle
  if (threadIdx.x < (32-4)){
    tmp = buff[threadIdx.x] + buff[threadIdx.x+4];
  }
  else {
    tmp = buff[threadIdx.x];
  }
  buff[threadIdx.x] = tmp;
  // // +8 cycle
  if (threadIdx.x < (32-8)){
    tmp = buff[threadIdx.x] + buff[threadIdx.x+8];
  }
  else {
    tmp = buff[threadIdx.x];
  }
  buff[threadIdx.x] = tmp;
  // // +16 cycle
  if (threadIdx.x < (32-16)){
    tmp = buff[threadIdx.x] + buff[threadIdx.x+16];
  }
  else {
    tmp = buff[threadIdx.x];
  }
  buff[threadIdx.x] = tmp;

  return buff[0];
}




#endif
