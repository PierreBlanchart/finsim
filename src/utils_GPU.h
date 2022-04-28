/** utils_GPU.h
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


#ifndef DEF_UTILS_GPU_H
#define DEF_UTILS_GPU_H

#include "def_utils.h"

#include <cuda_runtime.h>
#include <curand.h>


inline void cudaCheck(cudaError_t code) {
  if (code != cudaSuccess) {
    printf("GPUassert: %s\n", cudaGetErrorString(code));
  }
}

/*****************************************************************************************************/
// load data from host to GPU - and vice versa
template <typename T>
inline T *host2device(const T *data, const int sz) {
  if (sz <= 0) return NULL;
  T *data_GPU;
  cudaMalloc(&data_GPU, sz*sizeof(T));
  cudaMemcpy(data_GPU, data, sz*sizeof(T), cudaMemcpyHostToDevice);
  return data_GPU;
}

template <typename T>
inline void host2device(T *&data_GPU, const T *data, const int sz) {
  if (sz <= 0) return;
  if (!data_GPU) cudaMalloc(&data_GPU, sz*sizeof(T));
  cudaMemcpy(data_GPU, data, sz*sizeof(T), cudaMemcpyHostToDevice);
}


// load data from GPU to host
template <typename T>
inline T *device2host(const T *data, const int sz) {
  if (sz <= 0) return NULL;
  T *data_host = new T[sz];
  cudaMemcpy(data_host, data, sz*sizeof(T), cudaMemcpyDeviceToHost);
  return data_host;
}

template <typename T>
inline void device2host(T *&data_host, const T *data, const int sz) {
  if (sz <= 0) return;
  if (!data_host) data_host = new T[sz];
  cudaMemcpy(data_host, data, sz*sizeof(T), cudaMemcpyDeviceToHost);
}

#endif
