/** test_simu.cpp
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


#include "interface.h"


// [[Rcpp::export]]
arma::vec cuda_sim_BS_call(const arma::vec &S, const float H, const float E,
                           const float r, const float sigma,
                           const int Nsample) {

  int Npredict = S.n_elem;

  curandGenerator_t gen;
  curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);
  curandSetPseudoRandomGeneratorSeed(gen, 1234ULL); // set seed

  float *temp_gen;
  cudaMalloc(&temp_gen, Nsample*Npredict*sizeof(float));

  const float mu_H = (r - (sigma*sigma / 2.f))*H;
  const float sd_H = sigma*sqrt(H);
  curandGenerateNormal(gen, temp_gen, Nsample*Npredict, mu_H, sd_H);

  float *S_GPU = arma2GPU(S);
  float *res=NULL;
  fun_kernel_colMean_call(temp_gen, S_GPU, res, Nsample, Npredict, E);

  cudaFree(temp_gen);
  cudaFree(S_GPU);

  return GPU2arma_vec(res, Npredict);
}


// [[Rcpp::export]]
arma::vec cuda_sim_JDP_GPU(const float H,
                           const float r, const float sigma, // BM params
                           const double lambda_jump, const float mu_jump, const float sigma_jump, // Poisson jumps
                           const float dt, const int Nsample) {

  curandGenerator_t gen;
  curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);
  curandSetPseudoRandomGeneratorSeed(gen, 1234ULL); // set seed

  int Ntime = round(H/dt);
  printf("Simulating %d trajectories with %d time steps\n", Nsample, Ntime);

  // BM
  const float mu_H = (r - (sigma*sigma / 2.f))*dt, sd_H = sigma*sqrt(dt);
  float *temp_gen_BM;
  cudaMalloc(&temp_gen_BM, even_number(Ntime*Nsample)*sizeof(float)); // Ntime x Nsample
  curandGenerateNormal(gen, temp_gen_BM, even_number(Ntime*Nsample), mu_H, sd_H);

  // Poisson jumps
  unsigned int *NJ_per_traj;
  cudaMalloc(&NJ_per_traj, Nsample*sizeof(unsigned int));
  curandGeneratePoisson(gen, NJ_per_traj, Nsample, lambda_jump*H);

  unsigned int *temp_sum=NULL;
  unsigned int NJ_tot = fun_kernel_vecSum<unsigned int>(NJ_per_traj, temp_sum, Nsample);
  printf("%d jumps for %d trajectories\n", NJ_tot, Nsample);

  // positions of jumps
  float *temp_gen_pos;
  cudaMalloc(&temp_gen_pos, NJ_tot*sizeof(float));
  curandGenerateUniform(gen, temp_gen_pos, NJ_tot);

  // amplitude of jumps
  float *temp_val_JP;
  cudaMalloc(&temp_val_JP, even_number(NJ_tot)*sizeof(float));
  curandGenerateNormal(gen, temp_val_JP, even_number(NJ_tot), mu_jump, sigma_jump);
  // float *temp_val_JP_host = new float[NJ_tot];
  // cudaMemcpy(temp_val_JP_host, temp_val_JP, NJ_tot*sizeof(float), cudaMemcpyDeviceToHost);
  // printVec(temp_val_JP_host, NJ_tot, 3);

  // compute indexes of Poisson jumps - CPU version for now
  int *indexes = new int[NJ_tot];
  unsigned int *NJ_per_traj_host = new unsigned int[Nsample];
  cudaMemcpy(NJ_per_traj_host, NJ_per_traj, Nsample*sizeof(unsigned int), cudaMemcpyDeviceToHost);
  float *temp_gen_pos_host = new float[NJ_tot];
  cudaMemcpy(temp_gen_pos_host, temp_gen_pos, NJ_tot*sizeof(float), cudaMemcpyDeviceToHost);
  // printf("Sampling %d jumps for traj 1\n", NJ_per_traj_host[0]);
  // printVec(temp_gen_pos_host, NJ_per_traj_host[0], 3);

  int ind_cur = 0;
  for (int s=0; s < Nsample; s++) {
    for (int n=0; n < NJ_per_traj_host[s]; n++) {
      indexes[ind_cur] = int(temp_gen_pos_host[ind_cur]*float(Ntime)) + s*Ntime;
      ind_cur++;
    }
  }

  int *indexes_GPU;
  cudaMalloc(&indexes_GPU, NJ_tot*sizeof(int));
  cudaMemcpy(indexes_GPU, indexes, NJ_tot*sizeof(int), cudaMemcpyHostToDevice);
  delete[] indexes;
  delete[] NJ_per_traj_host;
  delete[] temp_gen_pos_host;

  // superimpose processes
  fun_kernel_add_with_index(temp_gen_BM, indexes_GPU, temp_val_JP, NJ_tot);

  cudaFree(indexes_GPU);
  cudaFree(NJ_per_traj);
  cudaFree(temp_val_JP);

  // sum each trajectory
  float *res=NULL;
  fun_kernel_colSum(temp_gen_BM, res, Ntime, Nsample);
  return GPU2arma_vec(res, Nsample);

  // return GPU2arma_vec(temp_gen_BM, Ntime);
}


