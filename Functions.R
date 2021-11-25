### PACKAGES

library(tidyverse)
library(deSolve)
library(gridExtra)
library(RColorBrewer)
library(reshape2)
library(Spbsampling)
library(viridis)
library(Matrix)

### FUNCTIONS

# Building matrices

BuildCircMat <-function(x) {
  n <- length(x)
  suppressWarnings(
    matrix(x[matrix(1:n,n+1,n+1,byrow=T)[c(1,n:2),1:n]],n,n))
}

BuildConstSumMat <- function(S, mu, bC) {
  mat <- matrix(runif(S^2, mu - bC, mu + bC), S, S)
  mat <- 0.5 * (mat + t(mat))
  diag(mat) <- 0
  ret <- stsum(mat, rep((S-1) * mu, times = S))$mat
  return(ret)
}

BuildLowerTriMat <- function(S, mu, bC) {
  mat <- matrix(runif(S^2, mu - bC, mu + bC), S, S)
  mat[upper.tri(mat)] <- 0
  return(mat)
}

BuildBandMat <- function(S, mu, bC) {
  mat <- matrix(runif(S^2, mu - bC, mu + bC), S, S)
  mat <- as.matrix(band(mat, 1, 1)) + as.matrix(band(mat, -1, -1))
  mat <- diag((1 / rowSums(mat)), S, S) %*% mat
  mat <- mat * mu
  return(mat)
}

BuildSparseMat <- function(S, mu, bC) {
  mat <- matrix(0, S, S)
  for(i in 1:S) {
    rand_size <- sample(x = 1:floor(S/3), size = 1)
    inds <- sample(x = setdiff(1:S, i), replace = FALSE, size = rand_size)
    mat[i,inds] <- 1
  }
  #adj_mu <- mu * (S-1) * S / sum(mat != 0)
  mat <- mat * matrix(runif(S^2, mu - bC, mu + bC), S, S)
  return(mat)
}

BuildCorrC <- function(S, mu, bC, p = 0) {
  mat <- matrix(0, S, S)
  for(i in 1:S) {
    for(j in i:S) {
      cur_sam <- p * runif(1, min = mu - bC, max = mu + bC)
      mat[i, j] <- cur_sam + (1 - p) * runif(1, min = mu - bC, max = mu + bC)
      mat[j, i] <- cur_sam + (1 - p) * runif(1, min = mu - bC, max = mu + bC)
    }
  }
  return(mat)
}

BuildCFromP <- function(S, Cd, C0, bC, P) {
  
  
}

BuildC <- function(S, MatType, Cd, C0, bC, p = 0) {
  
  if(MatType == "Random") {
    
    C <- matrix(runif(S^2, C0 - bC, C0 + bC), S, S)
    diag(C) <- Cd
    
  } else if(MatType == "Symmetric") {
    
    C <- matrix(runif(S^2, C0 - bC, C0 + bC), S, S)
    C <- 0.5 * (C + t(C))
    diag(C) <- Cd
    
    
  } else if(MatType == "ConstCirc") {
    
    c_vals <- seq(0, 1, length.out = S) * 1:S
    c_vals <- c_vals / max(c_vals) * C0
    C <- BuildCircMat(c_vals)
    diag(C) <- Cd
    C <- 0.5 * (C + t(C))
    
  } else if(MatType == "Circulant") {
    
    C <- BuildCircMat(runif(S, C0 - bC, C0 + bC))
    diag(C) <- Cd
    C <- 0.5 * (C + t(C))
    
  } else if(MatType == "Tradeoff") {
    
    C <- BuildConstSumMat(S, mu = C0, bC = bC)
    diag(C) <- Cd
    
  } else if(MatType == "Lower Triangular") {
    
    C <- BuildLowerTriMat(S, mu = C0, bC = bC)
    diag(C) <- Cd
    
  } else if(MatType == "Banded") {
    
    C <- BuildBandMat(S, mu = C0, bC = bC)
    diag(C) <- Cd
    
  } else if(MatType == "Sparse") {
    C <- BuildSparseMat(S, mu = C0, bC = bC)
    diag(C) <- Cd
    
  } else if(MatType == "Non-Sym Tradeoff") {
    C <- matrix(runif(S^2, min = C0 - bC, max = C0 + bC), S, S)
    C <- diag((C0 * (S-1) / rowSums(C)), S, S) %*% C
    diag(C) <- Cd
  } else if(MatType == "Correlated") {
    C <- BuildCorrC(S, mu = C0, bC = bC, p = p)
    diag(C) <- Cd
  }
  
  return(C)
}

BuildP <- function(S, MatType) {
  
  if(MatType == "Random") {
    
    P <- matrix(runif(S^2, 0, 1), S, S)
    diag(P) <- 0
    P <- diag((1 / rowSums(P)), S, S) %*% P
    
  } else if(MatType == "Symmetric") {
    
    P <- BuildConstSumMat(S, mu = 1 / (S-1), 1 / (S-1))
    
  } else if(MatType == "ConstCirc") {
    
    p_vals <- seq(0, 1, length.out = S) * 1:S
    P <- BuildCircMat(p_vals)
    diag(P) <- 0
    P <- 0.5 * (P + t(P))
    P <- diag((1 / rowSums(P)), S, S) %*% P
    
  } else if(MatType == "Circulant") {
    
    P <- BuildCircMat(runif(S, 0, 1))
    diag(P) <- 0
    P <- 0.5 * (P + t(P))
    P <- diag((1 / rowSums(P)), S, S) %*% P
    
  } else if(MatType == "Constant") {
    P <- matrix(1 / (S-1), S, S)
    diag(P) <- 0
  } else if(MatType == "Upper Triangular") {
    P <- t(BuildLowerTriMat(S, mu = 0.5, bC = 0.5))
    P <- diag((1 / rowSums(P)), S, S) %*% P
  } else if(MatType == "No Cross-Feeding") {
    P <- matrix(0, S, S)
  }
  return(P)
}

# Building parameter sets

BuildPars <- function(S, Ctype, Cd, C0, bC, p = 0, Ptype, eps_val, n, r) {
  
  C <- BuildC(S = S, MatType = Ctype, Cd = Cd, C0 = C0, bC = bC, p = p)
  P <- BuildP(S = S, MatType = Ptype)
  eps <- matrix(1, S, S)
  diag(eps) <- eps_val
  
  N <- rep(n, times = S)
  R <- rep(r, times = S)
  
  pars <- list(S = S, C = C, P = P, eps = eps, N = N, R = R)
  return(pars)
}

# Computing stability and feasibility

BuildJ <- function(pars) {
  with(pars, {
    
    Ctilde <- C
    diag(Ctilde) <- diag(C) * (1 - diag(eps))
    
    J <- matrix(0, 2*S, 2*S)
    J[1:S, 1:S] <- - diag(as.vector(C %*% N), S, S) + t(P) %*% diag(as.vector(Ctilde %*% N), S, S)
    J[1:S, (S+1):(2*S)] <- - diag(R, S, S) %*% C + t(P) %*% diag(R, S, S) %*% Ctilde
    J[(S+1):(2*S), 1:S] <- diag(diag(C) * N * diag(eps), S, S)
    
    return(J)
  })
}

GetEig <- function(J) {
  eigs <- eigen(J, only.values = TRUE)$values
  re_eigs <- Re(eigs)
  max_eig <- max(re_eigs)
  return(max_eig)
}

GetStab <- function(Cd, pars) {
  diag(pars$C) <- Cd
  J <- BuildJ(pars)
  eig <- GetEig(J)
  return(eig)
}

GetBound <- function(Cd, pars) {
  diag(pars$C) <- Cd
  bound <- with(pars, {
    Ctilde <- C
    diag(Ctilde) <- diag(C) * (1 - diag(eps))
    B <- - C + t(P) %*% Ctilde
    diagB <- diag(B)
    B <- abs(B)
    diag(B) <- diagB
    bound <- max(rowSums(B))
    return(bound)})
  return(bound)
}

GetFeas <- function(Cd, pars) {
  diag(pars$C) <- Cd
  feasibility <- with(pars, {
    Ctilde <- C
    diag(Ctilde) <- diag(C) * (1 - diag(eps))
    rho <- (diag(R, S, S) %*% C - t(P) %*% diag(R, S, S) %*% Ctilde) %*% N
    feasibility <- min(rho)
    return(feasibility)})
  return(feasibility)
}

StabB <- function(Cd, pars) {
  eig <- with(pars, {
    diag(C) <- Cd
    Ctilde <- C
    diag(Ctilde) <- diag(C) * (1 - diag(eps))
    B <- - C + t(P) %*% Ctilde
    Bprime <- B %*% diag(diag(eps) * n * Cd)
    return(GetEig(B))
  })
  return(eig)
}

FindCds <- function(pars, maxCd_mult) {
  Cd <- 1
  max_val <- maxCd_mult * max(pars$C)
  
  max_bound <- GetBound(max_val, pars)
  if(max_bound < 0) {
    bound <- uniroot(GetBound, interval = c(0, max_val), pars)$root
    Bpred <- uniroot(StabB, interval = c(0, max_val), pars)$root
  } else {
    bound <- NA
    Bpred <- NA
    print("Gershgorin bound or B always unstable")
  }
  
  min_feas <- GetFeas(0, pars)
  max_feas <- GetFeas(max_val, pars)
  
  if(min_feas * max_feas < 0) {
    feas <- uniroot(GetFeas, interval = c(0, max_val), pars)$root
  } else if(max_feas < 0) {
    feas <- NA
    print("Never feasible")
  } else if(min_feas >= 0) {
    feas <- 0
  }
  
  min_Cd <- GetStab(1e-10, pars)
  max_Cd <- GetStab(max_val, pars)
  
  if(min_Cd * max_Cd < 0) {
    Cd <- uniroot(GetStab, interval = c(1e-10, max_val), pars)$root
  } else if(max_Cd > 0) {
    Cd <- NA
    print("Never stable")
  }
  
  Cds <- data.frame(EigCd = Cd, Bound = bound, Bpred = Bpred, Feasibility = feas, Cd = max(Cd, feas))
  return(Cds)
}

GetR <- function(eta, pars) {
  R <- with(pars, eta / diag(C) / diag(eps))
  return(R)
}

GetN <- function(rho, pars) {
  N <- with(pars, {
    Ctilde <- C
    diag(Ctilde) <- diag(C) * (1 - diag(eps))
    N <- solve(diag(R, S, S) %*% C - t(P) %*% diag(R, S, S) %*% Ctilde) %*% rho
    return(N)}) 
  N <- as.vector(N)
  return(N)
}

GetFeasAndStab <- function(rho, pars) {
  N <- GetN(rho, pars)
  #print(N)
  feas <- prod(N > 0)
  pars$N <- N
  #print(feas)
  if(feas) {
    stab <- GetEig(BuildJ(pars))
    #print(stab)
    stab <- (stab < 0)
    #print(stab)
  } else {
    stab <- 0
  }
  S <- length(rho)
  
  grs <- with(pars, diag(R, S, S) %*% (eps * C))
  grs[grs == 0] <- NA
  sel_grs <- apply(grs, 2, FUN = min, na.rm = TRUE)
  inds <- which(grs == t(matrix(sel_grs, S, S)), arr.ind = TRUE)
  LLfeas <- prod(inds[,1] == 1:S)
  
  upper_band <- as.matrix(band(grs, 1, 1))
  upper_band <- upper_band[upper_band != 0]
  upper_band <- c(grs[S, 1], upper_band)
  
  lower_band <- as.matrix(band(t(grs), 1, 1))
  lower_band <- lower_band[lower_band != 0]
  lower_band <- c(lower_band, grs[1, S])
  
  subset_grs <- cbind(diag(grs), upper_band, lower_band)
  sel_sub_grs <- apply(subset_grs, 1, FUN = min, na.rm = TRUE)
  inds <- which(grs == t(matrix(sel_sub_grs, S, S)), arr.ind = TRUE)
  sub_feas <- prod(inds[,1] == 1:S)
  
  return(data.frame(Feas = feas, Stab = stab, LLFeas = LLfeas,
                    SubFeas = sub_feas, LLSuccess = prod(feas, stab, LLfeas),
                    SubSuccess = prod(feas, stab, sub_feas)))
}

# Iterating over parameters

RunSimulation <- function(params, maxCd_mult, num_repl) {
  in_params <- do.call("rbind", replicate(num_repl, params, simplify = FALSE))
  out_data <- data.frame(EigCd = c(), Bound = c(), Feasibility = c(), Cd = c())
  
  for(cur_row in 1:nrow(in_params)) {
    cur_params <- in_params[cur_row,]
    pars <- with(cur_params, BuildPars(S = S, Ctype = Ctype, Cd = 0,
                                       C0 = C0, bC = bC, p = p, Ptype = Ptype,
                                       eps_val = eps_val, n = n, r = r))
    cur_data <- cbind(cur_params, FindCds(pars, maxCd_mult))
    out_data <- rbind(out_data, cur_data)
    if(cur_row %% num_repl == 0) {
      print(paste(cur_row /  num_repl, "/", nrow(in_params) / num_repl))
    }
  }
  return(out_data)
}

RunCorrSim <- function(params, maxCd_mult, num_repl) {
  in_params <- do.call("rbind", replicate(num_repl, params, simplify = FALSE))
  out_data <- data.frame(EigCd = c(), Bound = c(), Feasibility = c(), Cd = c())
  
  for(cur_row in 1:nrow(in_params)) {
    cur_params <- in_params[cur_row,]
    pars <- with(cur_params, BuildPars(S = S, Ctype = Ctype, Cd = 0,
                                       C0 = C0, bC = bC, p = p, Ptype = Ptype,
                                       eps_val = eps_val, n = n, r = r))
    #Ceigs <- eigen(pars$C, only.values = TRUE)$values
    #eig_C <- min(Re(Ceigs))
    #pars$C <- pars$C / eig_C
    Cds <- FindCds(pars, maxCd_mult)
    cur_Cd <- Cds$EigCd
    diag(pars$C) <- cur_Cd
    Ctilde <- pars$C
    diag(Ctilde) <- diag(pars$C) * (1 - diag(pars$eps))
    B <- - pars$C + t(pars$P) %*% Ctilde
    #print(B)
    Beigs <- eigen(B, only.values = TRUE)$values
    #print(eigen(pars$C, only.values = T)$values)
    #print(Beigs)
    Beig <- min(Re(Beigs))
    ImEig <- mean(abs(Im(Beigs)))
    VarC <- var(C)
    #print(ImEig)
    cur_data <- cbind(cur_params, FindCds(pars, maxCd_mult), data.frame(ImEig = ImEig))
    out_data <- rbind(out_data, cur_data)
    if(cur_row %% num_repl == 0) {
      print(paste(cur_row /  num_repl, "/", nrow(in_params) / num_repl))
    }
  }
  return(out_data)
}


RunHeatMap <- function(in_params, num_repl) {
  
  in_params <- do.call("rbind", replicate(num_repl, in_params, simplify = FALSE))
  out_data <- data.frame()
  
  for(cur_row in 1:nrow(in_params)) {
    cur_params <- in_params[cur_row,]
    
    eta <- rep(1, times = cur_params$S)
    
    pars <- with(cur_params, BuildPars(S = S, Ctype = Ctype, Cd = Cd,
                                       C0 = C0, bC = bC, Ptype = Ptype,
                                       eps_val = eps_val, n = 0, r = 0))
    pars$R <- GetR(eta, pars)
    rho_const <- rep(1 / cur_params$S, times = cur_params$S)
    const_data <- GetFeasAndStab(rho_const * cur_params$rho_mult, pars)
    
    pars <- with(cur_params, BuildPars(S = S, Ctype = Ctype, Cd = Cd,
                                       C0 = C0, bC = bC, Ptype = Ptype,
                                       eps_val = eps_val, n = 0, r = 0))
    pars$R <- GetR(eta, pars)
    rho_rand <- runif(cur_params$S)
    rho_rand <- rho_rand / sum(rho_rand)
    rand_data  <- GetFeasAndStab(rho_rand * cur_params$rho_mult, pars)
    
    pars <- with(cur_params, BuildPars(S = S, Ctype = Ctype, Cd = Cd,
                                       C0 = C0, bC = bC, Ptype = Ptype,
                                       eps_val = eps_val, n = 0, r = 0))
    pars$R <- GetR(eta, pars)
    rho_one <- c(1, rep(0, times = cur_params$S - 1))
    one_data   <- GetFeasAndStab(rho_one * cur_params$rho_mult, pars)
    
    cur_data <- rbind(const_data, rand_data, one_data)
    cur_data$Rho <- c("Constant", "Random", "One Resource")
    
    cur_data <- cbind(cur_params, cur_data)
    out_data <- rbind(out_data, cur_data)
    
    if(cur_row %% num_repl == 0) {
      print(paste(cur_row /  num_repl, "/", nrow(in_params) / num_repl))
    }
  }
  return(out_data)
}

# Dynamics

dynamics <- function(time, state, params) {
  
  S <- pars$S
  C <- params$C
  P <- params$P
  rho <- params$rho
  eta <- params$eta
  eps <- params$eps
  dyn <- params$dyn
  
  resources <- state[1:S]
  consumers <- state[(S+1):(2*S)]
  
  R <- diag(resources, S, S)
  grs <- R %*% (eps * C)
  
  if(dyn == "Constant") {
    sel_grs <- diag(grs)
  } else if(dyn == "Minimum") {
    grs[grs == 0] <- NA
    sel_grs <- apply(grs, 2, FUN = min, na.rm = TRUE)
  } else if(dyn == "Maximum") {
    grs[grs == 0] <- NA
    sel_grs <- apply(grs, 2, FUN = max, na.rm = TRUE)
  } else if(dyn == "SubMin") {
    
    upper_band <- grs[cbind(1:(S-1), 2:S)]
    upper_band <- c(grs[S, 1], upper_band)
    upper_band[upper_band == 0] <- NA
    
    lower_band <- grs[cbind(2:S, 1:(S-1))]
    lower_band <- c(lower_band, grs[1, S])
    lower_band[lower_band == 0] <- NA
    
    subset_grs <- cbind(diag(grs), upper_band, lower_band)
    sel_grs <- apply(subset_grs, 1, FUN = min, na.rm = TRUE)
  }

  inds <- which(grs == t(matrix(sel_grs, S, S)), arr.ind = TRUE)
  Ctilde <- C
  Ctilde[inds] <- C[inds] * (1 - eps[inds])
  
  dresdt  <- rho - (R %*% C - t(P) %*% R %*% Ctilde) %*% consumers
  dconsdt <- consumers * (sel_grs - eta)
  
  return (list(c(dresdt, dconsdt)))
}

integrate_dyn <- function(inistate, pars, endtime, fn, timestep = 1){
  times <- seq(0, endtime, by = timestep)
  timeseries <- as.data.frame(ode(inistate, times, fn, pars, maxsteps = 1e4))  
  return(timeseries)
}

# Plotting

GridM <- function(M, Title, in_color) {
  hm.palette <- colorRampPalette(brewer.pal(9, in_color), space='Lab')
  M.melted <- melt(t(M))
  gridM <- ggplot(M.melted, aes(x = Var1, y = dim(M)[1] - Var2, fill = value)) +
    geom_tile() +
    ggtitle(Title) +
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.grid = element_blank(),
          legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          plot.title = element_text(hjust = 0.5, size = 17.5)) +
    coord_equal() +
    scale_fill_gradientn(
      colours=hm.palette(100), limits = c(0, max(M))
    )
  return(gridM)
}

plotEigs <- function(J, Title = "Spectrum") {
  
  eigs <- eigen(J, only.values = TRUE)$values
  deJ <- data.frame(Real = Re(eigs), Imaginary = Im(eigs))
  
  plJ <- ggplot(data = deJ, aes(x = Real, y = Imaginary)) +
    geom_point(color = "black", alpha = 0.5, size = 2) +
    theme_bw() +
    ggtitle(Title) +
    geom_vline(xintercept = 0, color = "red")
  return(plJ)
}

plSeries <- function(series, S, title) {
  
  rseries <- series[,1:(S+1)] %>% mutate(Type = "Resources")
  cseries <- series[,-(2:(S+1))] %>% mutate(Type = "Consumers")
  names(cseries) <- names(rseries)
  series <- rbind(rseries, cseries)
  
  meltseries <- melt(series, id.vars = c("time", "Type"))
  pl <- ggplot() +
    geom_line(data = meltseries, aes(x = time, y = value, color = variable),
              size = 1.5, alpha = 0.75) +
    ggtitle(title) +
    xlab("Time") + ylab("Abundance") +
    theme_bw() +
    facet_wrap(~Type, scales = "free") +
    theme(legend.position = "none") +
    theme(text = element_text(size=20),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15),
          legend.text=element_text(size=15))
  return(pl)
}
