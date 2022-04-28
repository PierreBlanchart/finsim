/** interface.h
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


#ifndef INTERFACE_H
#define INTERFACE_H

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

#include "cuda_kernels.cuh"


// interface helper functions
inline float *arma2GPU(const arma::mat &A) {
  arma::fmat fA = arma::conv_to<arma::fmat>::from(A);
  float *res = host2device<float>(fA.memptr(), fA.n_rows*fA.n_cols);
  return res;
}


inline arma::vec GPU2arma_vec(float *data, const int N) {
  float *data_host = device2host<float>(data, N);
  arma::fvec fres(data_host, N);
  arma::vec res = arma::conv_to<arma::vec>::from(fres);
  return res;
}


#endif
