library(finsim)

S0 <- 110
H <- 4
dt <- 1e-4
E <- 120
r <- 0.02
sigma <- 0.2

lambda.jump <- 1
mu.jump <- 0.05
sigma.jump <- 0.1

Ntime <- round(H/dt)

obj.traj <- sim_JDP(S=S0, r=r, sigma=sigma,
                    lambda.jump=lambda.jump, mu.jump=mu.jump, sigma.jump=sigma.jump,
                    dt=dt, Ntime=Ntime,
                    Nsample=1)

plot_fincurve(obj.traj, ind2plot = 1, (0:(Ntime-1))*dt)

