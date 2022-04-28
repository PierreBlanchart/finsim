/** def_utils.h
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


#ifndef DEF_UTILS_H
#define DEF_UTILS_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <string>
#include <math.h>
#include <chrono>
#include <time.h>
#include <limits>
#include <ctype.h>
#include <float.h>
#include <stdarg.h>
#include <vector>
#include <random>
#include <numeric>
#include <algorithm>
#include <functional>
#include <stack>
#include <queue>
#include <cstring>
using namespace std;

// color display in terminal
#define KNRM  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"

typedef std::chrono::high_resolution_clock myclock;
typedef std::chrono::high_resolution_clock::time_point timepoint;
using namespace std::chrono;

// open mp
#include <omp.h>
#define NTHREADS_OMP 16 // open MP

// simd
#include <xmmintrin.h>

typedef unsigned char uchar;


inline int even_number(int N) {
  return N + (N%2);
}


inline void printVec(const float *vec, const int N, const int n_digits=2) {
  stringstream temp;
  temp << "%." << n_digits << "f ";
  for (size_t i=0; i < N; i++) printf(temp.str().c_str(), vec[i]);
  printf("\n");
}


#endif
