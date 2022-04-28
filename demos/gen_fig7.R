library(finsim)

S0 <- 110
H <- 2
dt <- 1e-3
E <- 120
r <- 0.02
sigma <- 0.2

lambda.jump <- 0.2
mu.jump <- 0
sigma.jump <- 0.1

Ntime <- round(H/dt)
Nsample <- 1e4

# simulate Nsample trajectories
obj.traj <- sim_JDP(S=S0, r=r, sigma=sigma,
                    lambda.jump=lambda.jump, mu.jump=mu.jump, sigma.jump=sigma.jump,
                    dt=dt, Ntime=Ntime,
                    Nsample=Nsample)

# plot
# plot_fincurve(obj.traj, ind2plot = 1:10, (0:(Ntime-1))*dt)

# compute value of barrier options for different barrier values using the previously simulated trajectories
S.range <- seq(50, 150, by=1)
N.eval <- length(S.range)
traj <- obj.traj$traj/S0
B <- c(130, 170, 210) # barrier values
N.B <- length(B)

R.barrier <- matrix(NA, N.B, N.eval)
for (b in 1:N.B) {
  pb <- txtProgressBar(min = 0, max = N.eval, style = 3)
  for (i in 1:N.eval) {
    R.barrier[b, i] <- barrier_call_European(S.range[i]*traj, E = E, B = B[b], r=r, H=H)
    setTxtProgressBar(pb, i)
  }
  close(pb)
}

# plot C(S, t) as a function of S for each barrier value
colors <- glasbey()
plot(S.range, R.barrier[1, ], type="l", lwd=2, col=colors[1], ylim=c(0, max(R.barrier)), ylab="Rbarrier", cex.lab=1.5)
for (b in 2:N.B) {
  lines(S.range, R.barrier[b, ], type="l", lwd=2, col=colors[b], cex.lab=1.5)
}
legend("topleft", inset=-1e-2, legend=paste0("B=", B), col=colors, lty=rep(1, N.B), cex=1.3, lwd=3)
grid()


