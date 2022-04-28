/** cuda_kernels.cu
 *
 * Copyright (C) 2022 Pierre BLANCHART
 * pierre.blanchart@gmail.com
 * CEA/LIST/DM2I/SID/LI3A
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 **/


#include "cuda_kernels.cuh"


/*****************************************************************************************************************/
// sum of elements inside vector A (1 x len)
template <typename T>
__global__ void kernel_vecSum(const T * __restrict__ A, T *res, const int len) {

  // blockDim.x = MAX_THREADS_PER_BLOCK
  int stride = gridDim.x * MAX_THREADS_PER_BLOCK;
  int tx = threadIdx.x;
  int bx = blockIdx.x;
  int tid = MAX_THREADS_PER_BLOCK*bx + tx;

  T sum = 0;
  for (int i = tid; i < len; i += stride) {
    sum += A[i];
  }

  __shared__ T AS[MAX_THREADS_PER_BLOCK];
  AS[tx] = sum;

  // ensure all shared loaded
  __syncthreads();

  // Final summation
  sum = 0;
  int n = len < blockDim.x ? len : blockDim.x;
  while (n > 1) {
    sum += (n & 1) ? AS[n - 1] : 0;
    n >>= 1; // n /= 2
    if (tx < n) {
      AS[tx] += AS[n + tx];
    }
    // ensure all shared summed
    __syncthreads();
  }

  if (tx==0) {
    res[bx] = sum+AS[0];
  }
}

template __global__ void kernel_vecSum<float>(const float * __restrict__ A, float *res, const int len);
template __global__ void kernel_vecSum<int>(const int * __restrict__ A, int *res, const int len);
template __global__ void kernel_vecSum<unsigned int>(const unsigned int * __restrict__ A, unsigned int *res, const int len);

template <typename T>
T fun_kernel_vecSum(const T *A, T *&temp_sum, const int len) {
  if (!temp_sum) cudaMalloc(&temp_sum, GRIDSIZE*sizeof(T));

  dim3 dimBlock(MAX_THREADS_PER_BLOCK);
  dim3 dimGrid(GRIDSIZE);
  kernel_vecSum<T><<<dimGrid, dimBlock>>>(A, temp_sum, len);
  kernel_vecSum<T><<<1, dimBlock>>>(temp_sum, temp_sum, GRIDSIZE);
  cudaDeviceSynchronize();

  T sum;
  cudaMemcpy(&sum, temp_sum, sizeof(T), cudaMemcpyDeviceToHost);
  return sum;
}

// instantiation
template float fun_kernel_vecSum<float>(const float *A, float *&temp_sum, const int len);
template int fun_kernel_vecSum<int>(const int *A, int *&temp_sum, const int len);
template unsigned int fun_kernel_vecSum<unsigned int>(const unsigned int *A, unsigned int *&temp_sum, const int len);


/*****************************************************************************************************************/
__global__ void kernel_vecCumsum_int(int *v1, int *v2, const int offset, const int len) {
  int stride_grid = gridDim.x * blockDim.x; // total number of threads
  int tid = blockDim.x * blockIdx.x + threadIdx.x;

  int temp;
  for (int i=tid; i < len; i+=stride_grid) {
    if (i < offset) {
      v2[i] = v1[i];
    } else {
      v2[i] = v1[i-offset] + v1[i];
    }
  }
  __syncthreads();

  for (int i=tid; i < len; i+=stride_grid) {
    // swap(v1, v2);
    if (i >= offset) {
      temp = v1[i];
      v1[i] = v2[i];
      v2[i] = temp;
    }
  }
  __syncthreads();

}

/**
__global__ void scan(float *g_odata, float *g_idata, int n) {
  extern __shared__ float temp[]; // allocated on invocation
  int thid = threadIdx.x;
  int pout = 0, pin = 1; // Load input into shared memory. This is exclusive scan, so shift right by one, and set first element to 0
  temp[pout*n + thid] = (thid > 0) ? g_idata[thid-1] : 0;
  __syncthreads();
  for (int offset = 1; offset < n; offset *= 2)   {
    pout = 1 - pout; // swap double buffer indices
    pin = 1 - pout;
    if (thid >= offset) temp[pout*n+thid] += temp[pin*n+thid - offset]; else temp[pout*n+thid] = temp[pin*n+thid];
    __syncthreads();
  }
  g_odata[thid] = temp[pout*n+thid]; // write output
}
**/

void fun_kernel_vecCumsum_int(int *v, const int len) {
  int *vbuff;
  cudaMalloc(&vbuff, len*sizeof(int));

  int offset = 1;
  while (offset < len) {
    kernel_vecCumsum_int<<<GRIDSIZE,MAX_THREADS_PER_BLOCK>>>(v, vbuff, offset, len);
    offset *= 2;
  }
  cudaFree(vbuff);
}


/*****************************************************************************************************************/
// sum per col of colwise order matrix A (nrow x ncol)
__global__ void kernel_colSum(const float * __restrict__ A, float *res, const int nrow) {

  __shared__ float AS[REDUCE_SIZE];

  int bx = blockIdx.x; // from 0 to ncol - 1
  int tx = threadIdx.x; // from 0 to REDUCE_SIZE - 1
  int offs = bx * nrow + tx;

  float sum = 0;

  for (int i = 0; i < nrow / REDUCE_SIZE; i++, offs += REDUCE_SIZE) {
    sum += A[offs];
  }
  // Sum the remaining part
  if ((nrow % REDUCE_SIZE) != 0) {
    if (tx < nrow % REDUCE_SIZE) {
      sum += A[offs];
    }
  }

  AS[tx] = sum;
  // ensure all shared loaded
  __syncthreads();

  // Final summation
  sum = 0;
  int n = nrow < REDUCE_SIZE ? nrow : REDUCE_SIZE;
  while (n > 1) {
    sum += (n & 1) ? AS[n - 1] : 0;
    n >>= 1; // n /= 2
    if (tx < n) {
      AS[tx] += AS[n + tx];
    }
    // ensure all shared summed
    __syncthreads();
  }

  if (tx==0) {
    res[bx] = sum+AS[0];
  }
}


void fun_kernel_colSum(const float *A, float *&res, const int nrow, const int ncol) {
  if (!res) cudaMalloc(&res, ncol*sizeof(float));
  if (ncol > MAX_GRIDSIZE) {
    int ind_cur=0;
    int Nchunks = ncol/MAX_GRIDSIZE;
    if (Nchunks*MAX_GRIDSIZE < ncol) Nchunks++;
    for (int c=0; c < Nchunks; c++) {
      int dim_chunk = min(MAX_GRIDSIZE, ncol-ind_cur);
      printf("Computing chunk %d [%d]\n", c, dim_chunk);
      dim3 dimGrid(dim_chunk);
      dim3 dimBlock(REDUCE_SIZE);
      kernel_colSum<<<dimGrid, dimBlock>>>(&A[ind_cur*nrow], &res[ind_cur], nrow);
      ind_cur += MAX_GRIDSIZE;
    }
  } else {
    dim3 dimGrid(ncol);
    dim3 dimBlock(REDUCE_SIZE);
    kernel_colSum<<<dimGrid, dimBlock>>>(A, res, nrow);
  }
}


/*************************************************************************************************************************/
// sum per col of colwise order matrix A (nrow x ncol)
__global__ void kernel_colMean_call(const float * __restrict__ A, const float * __restrict__ S_target,
                                    float *res, const int nrow, const float E) {

  __shared__ float AS[REDUCE_SIZE];

  int bx = blockIdx.x; // from 0 to ncol - 1
  int tx = threadIdx.x; // from 0 to REDUCE_SIZE - 1
  int offs = bx * nrow + tx;

  float sum = 0;

  for (int i = 0; i < nrow / REDUCE_SIZE; i++, offs += REDUCE_SIZE) {
    sum += max(S_target[bx]*exp(A[offs]) - E, 0.f);
  }
  // Sum the remaining part
  if ((nrow % REDUCE_SIZE) != 0) {
    if (tx < nrow % REDUCE_SIZE) {
      sum += max(S_target[bx]*exp(A[offs]) - E, 0.f);
    }
  }

  AS[tx] = sum;
  // ensure all shared loaded
  __syncthreads();

  // Final summation
  sum = 0;
  int n = nrow < REDUCE_SIZE ? nrow : REDUCE_SIZE;
  while (n > 1) {
    sum += (n & 1) ? AS[n - 1] : 0;
    n >>= 1; // n /= 2
    if (tx < n) {
      AS[tx] += AS[n + tx];
    }
    // ensure all shared summed
    __syncthreads();
  }

  if (tx==0) {
    res[bx] = (sum+AS[0])/float(nrow);
  }
}


void fun_kernel_colMean_call(const float *A, const float *S_target, float *&res, const int nrow, const int ncol, const float E) {
  if (!res) cudaMalloc(&res, ncol*sizeof(float));
  if (ncol > MAX_GRIDSIZE) {
    int ind_cur=0;
    int Nchunks = ncol/MAX_GRIDSIZE;
    if (Nchunks*MAX_GRIDSIZE < ncol) Nchunks++;
    for (int c=0; c < Nchunks; c++) {
      int dim_chunk = min(MAX_GRIDSIZE, ncol-ind_cur);
      printf("Computing chunk %d [%d]\n", c, dim_chunk);
      dim3 dimGrid(dim_chunk);
      dim3 dimBlock(REDUCE_SIZE);
      kernel_colMean_call<<<dimGrid, dimBlock>>>(&A[ind_cur*nrow], &S_target[ind_cur], &res[ind_cur], nrow, E);
      ind_cur += MAX_GRIDSIZE;
    }
  } else {
    dim3 dimGrid(ncol);
    dim3 dimBlock(REDUCE_SIZE);
    kernel_colMean_call<<<dimGrid, dimBlock>>>(A, S_target, res, nrow, E);
  }
}


/*****************************************************************************************************************/
__global__ void kernel_add_with_index(float *A, const int * __restrict__ index, const float * __restrict__ val, const int len_index) {
  int stride = gridDim.x * blockDim.x; // total number of threads
  int tid = blockDim.x * blockIdx.x + threadIdx.x;

  for (int i=tid; i < len_index; i+=stride) {
    A[index[i]] += val[i];
  }
}

void fun_kernel_add_with_index(float *A, const int *index, const float *val, const int len_index) {
  kernel_add_with_index<<<GRIDSIZE,MAX_THREADS_PER_BLOCK>>>(A, index, val, len_index);
}

