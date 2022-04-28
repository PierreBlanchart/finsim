# /** utils_fin.R
# *
# * Copyright (C) 2022 Pierre BLANCHART
# * pierre.blanchart@gmail.com
# * CEA/LIST/DM2I/SID/LI3A
# * This program is free software: you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation, either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program.  If not, see <https://www.gnu.org/licenses/>.
# **/


#' @export BS_call
#' BS: european call
BS_call <- function(S, H, E, r, sigma, PLOT_=FALSE) {
  log_S_E <- log(S/E)
  sigma_H <- sigma*sqrt(H)
  sigma_sq <- sigma^2

  d1 <- (log_S_E + (r + sigma_sq/2)*H) / sigma_H
  d2 <- d1 - sigma_H # (log_S_E + (r - sigma_sq/2)*H) / sigma_H

  N_d1 <- pnorm(d1)
  N_d2 <- pnorm(d2)
  C_S_H <- N_d1*S - N_d2*E*exp(-r*H)

  if (PLOT_) {
    plot(S, C_S_H, type="l")
    grid()
  }

  return(C_S_H)
}


#' @export BS_call.dividend
#' BS: european call with one discrete dividend
#' 0 <= D < 1
BS_call.dividend <- function(S, H, E, r, sigma, D, PLOT_=FALSE) {
  E.D <- E/(1-D)
  log_S_E <- log(S/E.D)
  sigma_H <- sigma*sqrt(H)
  sigma_sq <- sigma^2

  d1 <- (log_S_E + (r + sigma_sq/2)*H) / sigma_H
  d2 <- d1 - sigma_H # (log_S_E + (r - sigma_sq/2)*H) / sigma_H

  N_d1 <- pnorm(d1)
  N_d2 <- pnorm(d2)
  C_S_H <- (N_d1*S - N_d2*E.D*exp(-r*H))*(1-D)

  if (PLOT_) {
    plot(S, C_S_H, type="l")
    grid()
  }

  return(C_S_H)
}


#' @export S_distrib
#' log-normal distribution of S at time H
S_distrib <- function(S0, S.range, H, r, sigma, PLOT_=FALSE) {

  # underlying normal law parameters (distribution of log(S(H)))
  mu.N <- log(S0) + (r - (sigma^2)/2)*H
  sigma.N <- sigma*sqrt(H)

  S.density <- (1/(sigma.N*sqrt(2*pi)*S.range)) * exp(-(log(S.range) - mu.N)^2 / (2*sigma.N^2))

  if (PLOT_) {
    plot(S.range, S.density, type="l")
    grid()
  }

  return(
    list(
      density=S.density,
      mu=exp(mu.N + (sigma.N^2)/2),
      sigma=exp(sigma.N^2 - 1)*exp(2*mu.N + sigma.N^2)
    )
  )
}


colCumsum <- function(A) {
  Nrow <- nrow(A)
  Ncol <- ncol(A)
  temp.cum <- cumsum(A)
  if (Ncol==1) return(temp.cum)
  temp.res <- temp.cum[seq(Nrow, Nrow*Ncol, by=Nrow)]
  return(matrix(temp.cum - rep(c(0, temp.res[1:(Ncol-1)]), each=Nrow), Nrow))
}


myrep <- function(val, each) {
  Nval <- length(val)
  max.each <- max(each)+1
  temp.cum <- matrix(0, max.each, Nval)
  temp.cum[(0:(Nval-1))*max.each + each+1] <- 1
  ind.rep <- rep(1:Nval, each=max.each)[(1 - colCumsum(temp.cum)) > 0]
  return(val[ind.rep])
}


rowprod <- function(A) {
  return(apply(A, MARGIN = c(1), FUN = function(row) prod(row)))
}


rowCumprod <- function(A) {
  return(t(apply(A, MARGIN = c(1), FUN = function(row) cumprod(row))))
}


#' @export sim_JDP
sim_JDP <- function(S, r, sigma, lambda.jump, mu.jump, sigma.jump, dt, Ntime, Nsample) {

  # BM
  mu.BM <- (r - (sigma*sigma)/2)*dt
  sigma.BM <- sigma*sqrt(dt)
  val.BM <- matrix(exp(rnorm(Nsample*Ntime, mean=mu.BM, sd=sigma.BM)), Nsample) # Nsample x Ntime
  val.BM[, 1] <- S

  # poisson jumps
  H <- Ntime*dt
  NH <- rpois(Nsample, lambda.jump*H) # number of jumps per sampled trajectory

  N.jumps <- sum(NH)
  ind.sample <- myrep(1:Nsample, each=NH) # index of sample trajectory to which belong a given jump

  ind.jumps <- pmax(round(runif(N.jumps)*Ntime), 2) # index of jumps between 2 and Ntime for each sampled trajectory
  val.jumps <- exp(rnorm(N.jumps, mean=mu.jump, sd=sigma.jump))

  val.J <- matrix(1, Nsample, Ntime)
  val.J[(ind.jumps-1)*Nsample + ind.sample] <- val.jumps

  # traj.sim <- rowprod(val.BM*val.J)
  # return(traj.sim) # plot histo of this to get distribution

  traj.sim <- rowCumprod(val.BM*val.J)
  return(
    list(
      traj=traj.sim,
      jump2traj=ind.sample,
      pos.jumps=ind.jumps-1
    )
  )
}


#' @export barrier_call_European
#' traj: Nsample x Ntime
barrier_call_European <- function(traj, E, B, r, H) {
  Ntime <- ncol(traj)
  max.per.traj <- apply(traj, MARGIN = c(1), FUN=function(row) max(row))
  is.passing <- max.per.traj > B
  return(exp(-r*H)*mean(is.passing * pmax(traj[, Ntime]-E, 0)))
}


glasbey <- function() {
  return(
    c("#0000FF", "#FF0000", "#00FF00", "#000033", "#FF00B6", "#005300", "#FFD300", "#009FFF", "#9A4D42", "#00FFBE", "#783FC1",
      "#1F9698", "#FFACFD", "#B1CC71", "#F1085C",   "#FE8F42", "#DD00FF", "#201A01", "#720055", "#766C95", "#02AD24", "#C8FF00",
      "#886C00", "#FFB79F", "#858567", "#A10300", "#14F9FF", "#00479E", "#DC5E93", "#93D4FF", "#004CFF", "#F2F318")
  )
}


#' @export plot_with_interrupt
plot_with_interrupt <- function(ind2plot, obj.traj, seq.time, col, ymin=NA, ymax=NA) {
  Ntime <- length(seq.time)
  H <- seq.time[Ntime]

  ind.jumps <- which(obj.traj$jump2traj == ind2plot)
  Njumps <- length(ind.jumps) # number of jumps in this trajectory
  pos.jumps <- sort(c(obj.traj$pos.jumps[ind.jumps], Ntime))
  traj2plot <- obj.traj$traj[ind2plot, ]

  ind.cur <- 1
  for (i in 1:(Njumps+1)) {
    if (i==1 && !is.na(ymin)) {
      plot(
        seq.time[ind.cur:pos.jumps[i]], traj2plot[ind.cur:pos.jumps[i]],
        xlim=c(0, H), ylim=c(ymin, ymax),
        type="l",
        col=col,
        lwd=2,
        xlab="time",
        ylab="S",
        cex.lab=1.5
      )
    } else {
      lines(
        seq.time[ind.cur:pos.jumps[i]], traj2plot[ind.cur:pos.jumps[i]],
        type="l",
        col=col,
        lwd=2
      )
    }
    ind.cur <- pos.jumps[i]+1
  }

}


#' @export plot_fincurve
#' traj: Ntraj x Ntime
plot_fincurve <- function(obj.traj, ind2plot, seq.time=NULL, add.legend=FALSE) {
  Ntraj <- length(ind2plot)
  Ntime <- ncol(obj.traj$traj)

  col.choice <- glasbey()
  colors <- col.choice[sample(length(col.choice), Ntraj, replace=TRUE)]

  if (is.null(seq.time)) seq.time <- seq(0, 1, length.out=Ntime)

  plot_with_interrupt(ind2plot[1], obj.traj, seq.time, col=colors[1],
                      ymin=min(obj.traj$traj[ind2plot, ]),
                      ymax=max(obj.traj$traj[ind2plot, ]))
  grid(NA, NULL, col = "lightgray", lty = "dotted", lwd = 1) # horizontal grid

  if (Ntraj > 1) {
    for (t in 2:Ntraj) {
      plot_with_interrupt(ind2plot[t], obj.traj, seq.time, col=colors[t])
    }
  }

  if (add.legend) {
    legend("topright", inset=-1e-2, legend=paste0("traj_", ind2plot), col=colors, lty=rep(1, Ntraj), cex=1, lwd=3)
  }
}




