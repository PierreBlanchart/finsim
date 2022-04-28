library(finsim)

S0 <- 110
H <- 2
dt <- 1e-3
E <- 120
r <- 0.02
sigma <- 0.2

lambda.jump <- 1
mu.jump <- 0.05
sigma.jump <- 0.1

Ntime <- round(H/dt)
Nsample <- 5e4

obj.traj <- sim_JDP(S=S0, r=r, sigma=sigma,
                    lambda.jump=lambda.jump, mu.jump=mu.jump, sigma.jump=sigma.jump,
                    dt=dt, Ntime=Ntime,
                    Nsample=Nsample)

# distrib at t=H
par(
  mfrow=c(1, 2),
  mar=c(2, 3, 1, 1),
  xaxt='n'
)
hist(obj.traj$traj[, Ntime], breaks=80, freq=FALSE, main="", xlim=c(0, 300), ylim=c(0, 1.3e-2))
grid()
bw <- mean(abs(diff(obj.traj$traj[, Ntime]))) / 4
obj.density <- density(obj.traj$traj[, Ntime], bw=bw, kernel="gaussian")
plot(obj.density$x, obj.density$y, xlim=c(0, 300), ylim=c(0, 1.3e-2), type="l", bty="n")
grid()

