library(finsim)

# comparison between the analytical call option return as computed using the BS analytical formula, and the estimated return obtained by GPU simulation
H <- 2
E <- 120
r <- 0.02
sigma <- 0.2

S.range <- c(0, 200)
S.val <- as.numeric(seq(S.range[1], S.range[2], 1e-1))

Nsample <- 1e5
res <- cuda_sim_BS_call(S = S.val, H = H, E = E, r = r, sigma = sigma, Nsample = Nsample) # GPU simulation
res.comp <- BS_call(S = S.val, H = H, E = E, r = r, sigma = sigma) # analytical BS
par(
  mfrow=c(1, 2),
  mar=c(2, 3, 1, 1)
)
plot(S.val, res.comp, type="l", main="analytical")
grid()
plot(S.val, exp(-r*H)*res, type="l", col="red", main="simulated")
grid()
