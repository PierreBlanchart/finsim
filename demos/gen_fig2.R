library(finsim)

H <- 2
r <- 0.02
sigma <- 0.1

S.range <- c(30, 180)
S.val <- as.numeric(seq(S.range[1], S.range[2], 1e0))

E.range <- c(80, 150)
E.val <- as.numeric(seq(E.range[1], E.range[2], 1e0))

surf.BS <- matrix(NA, length(E.val), length(S.val))
for (i in 1:length(E.val)) {
  surf.BS[i, ] <- BS_call(S = S.val, H = H, E = E.val[i], r = r, sigma = sigma, PLOT_ = FALSE)
}

library(plotly)
fig <- plot_ly(x = S.val, y = E.val, z = surf.BS, type = "surface") %>%
  layout(
    title="European Call",
    scene= list(
      xaxis=list(title="S.val"),
      yaxis=list(title="E.val"),
      zaxis=list(title="C.val")
    )
  )
fig


