/** cuda_kernels.cuh
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


#ifndef MY_CUDA_KERNELS_H
#define MY_CUDA_KERNELS_H

#include "utils_GPU.h"


#define MAX_THREADS_PER_BLOCK 1024
#define GRIDSIZE 48
#define MAX_GRIDSIZE 65520
#define REDUCE_SIZE 64


// sum reduce
template <typename T>
__global__ void kernel_vecSum(const T * __restrict__ A, T *res, const int len);

template <typename T>
T fun_kernel_vecSum(const T *A, T *&temp_sum, const int len);

// utils
void fun_kernel_vecCumsum_int(int *v, const int len);

// colsum reduce
void fun_kernel_colSum(const float *A, float *&res, const int nrow, const int ncol);
void fun_kernel_colMean_call(const float *A, const float *S_target, float *&res, const int nrow, const int ncol, const float E);

void fun_kernel_add_with_index(float *A, const int *index, const float *val, const int len_index);

#endif
