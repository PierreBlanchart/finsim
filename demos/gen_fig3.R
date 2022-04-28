library(finsim)

H.max <- 3
E <- 180
r <- 0.02
sigma <- 0.2

t1 <- 2
D1 <- 0.05

step <- 1e-2
S.range <- seq(50, 800, by=0.01)

res.H.max <- BS_call.dividend(S = S.range, H = H.max, E = E, r = r, sigma = sigma, D=D1, PLOT_ = FALSE)
res.t1 <- BS_call.dividend(S = S.range, H = t1, E = E, r = r, sigma = sigma, D=D1, PLOT_ = FALSE)*(1+D1)*exp(r*(H.max-t1))

# compute optimal exercise price
plot(S.range, res.H.max, ylim = c(0, max(c(res.H.max, res.t1))), type='l')
lines(S.range, res.t1, type='l', col='red')
grid()

diff.A <- res.t1-res.H.max
S.crossing <- S.range[round(which.min(abs(diff.A)))]
plot(S.range, diff.A, type='l', lwd=2, xlab="S", ylab = "Exercise - No_exercise", cex.lab=1.5)
abline(v=S.crossing, col="blue", lty=2, lwd=2)
text(S.crossing, 34, paste0("Optimal exercise price:\nS_opt = ", round(S.crossing, 1)), pos=4, col="blue", cex=1.3)
grid()

