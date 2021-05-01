#########################################################################################
## BI.r
## Blinding Assessment Indices for Randomized, Controlled, Clinical Trials

## Copyright 2010 - 2020: Original 2010 R Code by Nate Mercaldo (nmercaldo@mgh.harvard.edu)
## Copyright 2021: Updates by Marc Schwartz (marc_schwartz@me.com)

## This software is distributed under the terms of the GNU General Public License
## Version 3, June 2007.

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see https://www.gnu.org/licenses/gpl-3.0.html
#########################################################################################




BI <- function(x, weights = NULL, conf.level = 0.95,
               alternative.J = c("two.sided", "less", "greater"),
               alternative.B = c("two.sided", "less", "greater"),
               group.names = c("Treatment", "Placebo")) {

  
  ## The format of the 3 x 2 'x' count matrix,
  ## or the 3 x 2 'weights' matrix, if specified, should be: 
  ##              Treatment  Placebo
  ## Treatment    xxx        xxx 
  ## Placebo      xxx        xxx 
  ## Don't Know   xxx        xxx
  
  ## where the rows are the assignment guesses by the surveyed party,
  ## and the columns are the actual assignments
  
  if (!identical(dim(x), c(3L, 2L))) {
    stop("'x' must be a 3 row by 2 column integer matrix of cross-tabulated counts. At present only 3 blinding query response levels and 2 treatment arms are supported.")
  }

  
  ## Use default 1996 James weights for 3 x 2 unless alternative
  ## weights are specified. Check dim(weights) otherwise.
  ## Correct guesses are assigned a weight of 0
  ## Incorrect guesses are assigned a weight of 0.5
  ## Don't Know guesses are assigned a weight of 1
  ## Thus, row and column ordering is critical here
  ## and must match the above structure 1:1 for the count matrix 'x'
  if (is.null(weights)) {
    weights <- matrix(c(0, 0.5, 0.5, 0, 1, 1), nrow = 3, ncol = 2, byrow = TRUE)
  } else if (!identical(dim(x), c(3L, 2L))) {
    stop("'weights' must be a 3 row by 2 column numeric matrix specifying alternative James weights for each cell in 'x'.")
  }

  
  if (length(conf.level) != 1 | (conf.level <= 0) | (conf.level >= 1)) { 
    stop("'conf.level' must be a scalar numeric value >0 and <1")
  }  

  alternative.J <- match.arg(alternative.J)
  alternative.B <- match.arg(alternative.B)

  if (alternative.J == "two.sided") {
    alpha.J <- (1 - conf.level) / 2
    Sided.J <- "2-Sided"
  } else if (alternative.J == "less") {
    alpha.J <- 1 - conf.level
    Sided.J <- "1-Sided"
  } else if (alternative.J == "greater") {
    alpha.J <- 1 - conf.level
    Sided.J <- "1-Sided"
  }

  if (alternative.B == "two.sided") {
    alpha.B <- (1 - conf.level) / 2
    Sided.B <- "2-Sided"
  } else if (alternative.B == "less") {
    alpha.B <- 1 - conf.level
    Sided.B <- "1-Sided"
  } else if (alternative.B == "greater") {
    alpha.B <- 1 - conf.level
    Sided.B <- "1-Sided"
  }

  CI.J <- qnorm(alpha.J, lower.tail = FALSE)
  CI.B <- qnorm(alpha.B, lower.tail = FALSE)

  
  #########################################################################################
  ## Compute James' Blinding Index
  ## Use 'x' here
  x1 <- addmargins(x)
  
  P <- x1 / max(x1)
  
  Pdk <- P[nrow(P) - 1, ncol(P)]
  
  Pdo <- Pde <- v <- term1.denom <- 0
  
  for (i in 1:(nrow(P) - 2)) {
    
    for (j in 1:(ncol(P) - 1)) {
      
      Pdo <- Pdo + ((weights[i, j] * P[i, j]) / (1 - Pdk))
      
      Pde <- Pde + ((weights[i, j] * P[i, ncol(P)] * (P[nrow(P), j] - P[nrow(P) - 1, j])) / (1 - Pdk) ^ 2)
      
      term1.denom <- term1.denom + weights[i,j] * P[i, ncol(P)] * (P[nrow(P), j] - P[nrow(P) - 1, j])
    }
  }
  
  Kd <- (Pdo - Pde) / Pde
  
  term1.denom <- 4 * term1.denom ^ 2
  
  term1.num <- 0
  
  for (i in 1:(nrow(P) - 2)) {
    
    for (j in 1:(ncol(P) - 1)) {
      
      extra <- 0
      
      for (r in 1:(ncol(P) - 1)) {
        
        extra <- extra + (weights[r, j] * P[r, ncol(P)] + weights[i, r] * (P[nrow(P), r] - P[nrow(P) - 1, r]))
      }
      
      term1.num <- term1.num + ((P[i, j] * (1 - Pdk) ^ 2 * ((1 - Pdk) * weights[i, j] - (1 + Kd) * extra) ^ 2))
    }
  }
  
  v.1 <- term1.num / term1.denom
  
  v.2 <- Pdk * (1 - Pdk) - (1 - Pdk) * (1 + Kd) * (Pdk + ((1 - Pdk) * (1 + Kd)) / 4)
  
  v <- (v.1 + v.2) / x1[nrow(x1), ncol(x1)]
  

  BI.est <- (1 + Pdk + (1 - Pdk) * Kd) / 2
  
  BI.se <- sqrt(v)

  if (alternative.J == "two.sided") {
    BI.CI.Lower <- BI.est - (CI.J * BI.se)
    BI.CI.Upper <- BI.est + (CI.J * BI.se)
  } else if (alternative.J == "less") {
    BI.CI.Lower <- 0
    BI.CI.Upper <- BI.est + (CI.J * BI.se)
  } else if (alternative.J == "greater") {
    BI.CI.Lower <- BI.est - (CI.J * BI.se)
    BI.CI.Upper <- 1
  }

  BI.James <- matrix(c(BI.est, BI.se,
                       BI.CI.Lower,
                       BI.CI.Upper),
                     nrow = 1)
  
 
  #########################################################################################
  ## Compute Bang Blinding Index
  ## Use t(x) here
  x2 <- addmargins(t(x))

  ## pre-allocate vectors for loop results
  BI.est <- numeric(2)
  BI.se <- numeric(2)
  
  for (i in 1:(nrow(x2) - 1)) {
    
    BI.est[i] <- (2 * (x2[i, i] / (x2[i, 1] + x2[i, 2])) - 1) *
                 ((x2[i, 1] + x2[i, 2]) / (x2[i, 1] + x2[i, 2] + x2[i, 3]))
    
    BI.se[i] <- sqrt(((x2[i, 1] / x2[i, ncol(x2)]) * (1 - (x2[i, 1] / x2[i, ncol(x2)])) +
                     (x2[i, 2] / x2[i, ncol(x2)]) * (1 - (x2[i, 2] / x2[i, ncol(x2)])) +
                     2 * (x2[i, 1] / x2[i, ncol(x2)]) * (x2[i, 2] / x2[i, ncol(x2)])) / x2[i, ncol(x2)])
  }

  if (alternative.B == "two.sided") {
    BI.CI.Lower <- BI.est - (CI.B * BI.se)
    BI.CI.Upper <- BI.est + (CI.B * BI.se)
  } else if (alternative.B == "less") {
    BI.CI.Lower <- c(-1, -1)
    BI.CI.Upper <- BI.est + (CI.B * BI.se)
  } else if (alternative.B == "greater") {
    BI.CI.Lower <- BI.est - (CI.B * BI.se)
    BI.CI.Upper <- c(1, 1)
  }
  
  
  BI.Bang <- matrix(c(BI.est,
                      BI.se,
                      BI.CI.Lower,
                      BI.CI.Upper),
                     nrow = 2)


  
  #########################################################################################
  ## Format for output
  
  rownames(BI.James) <- "Overall"
  
  colnames(BI.James) <- c("Estimate", "Std. Error",
                          paste(conf.level * 100, "% ", c("LCL", "UCL"), " (", Sided.J, ")", sep = ""))

  rownames(BI.Bang) <- group.names
  
  colnames(BI.Bang) <- c("Estimate", "Std. Error",
                          paste(conf.level * 100, "% ", c("LCL", "UCL"), " (", Sided.B, ")", sep = ""))

  
  ## list results 
  BI <- list(BI.James, BI.Bang)
  names(BI) <- c("JamesBI", "BangBI")

  ## Return list
  BI
}
