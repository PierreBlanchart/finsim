library(finsim)

S0 <- 110
H <- 4
dt <- 1e-3
E <- 120
r <- 0.02
sigma <- 0.2

lambda.jump <- 0
mu.jump <- 0
sigma.jump <- 0

Ntime <- round(H/dt)

obj.traj <- sim_JDP(S=S0, r=r, sigma=sigma,
                    lambda.jump=lambda.jump, mu.jump=mu.jump, sigma.jump=sigma.jump,
                    dt=dt, Ntime=Ntime,
                    Nsample=10)

plot_fincurve(obj.traj, ind2plot = 1:10, (0:(Ntime-1))*dt)

