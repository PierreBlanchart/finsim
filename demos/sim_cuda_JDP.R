library(finsim)

S0 <- 110
H <- 2
dt <- 1e-2
r <- 0.02
sigma <- 0.2

lambda.jump <- 1
mu.jump <- 0.05
sigma.jump <- 0.2

Nsample <- 1e5
res <- as.numeric(cuda_sim_JDP_GPU(H=H, r=r, sigma=sigma, lambda_jump=lambda.jump, mu_jump=mu.jump, sigma_jump=sigma.jump, dt=dt, Nsample=Nsample))

hist(S0*exp(res), breaks=100, freq=FALSE, main="")
grid()


