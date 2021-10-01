
ToxCrit <- function (npat, ptoxcrit, pstop, nsim, ptox, a=1, b=1) {


  #---------------------------------------------------------------------------
  # Fix problem concerning pseudo-random numbers:
  set.seed(as.integer(Sys.time()))

  # Preparatory commands:
  ptoxcrit <- ptoxcrit/100
  pstop <- pstop/100
  ptox <- ptox/100
  crit <- rep((npat + 1), npat)
  # Initially, no stopping is planned; i.e., the critical boundary is
  # > npat

  # Modification according to the posterior probability (Beta distribution):

  for (k in (1:npat)) {
    for (l in (0:k)) {
      z <- qbeta((1 - pstop), a + l, b + k - l)
      if (z > ptoxcrit) {
        crit[k] <- l
        break
      }
    }
  }

  # Simulation:
  nstop <- rep(0, npat)
  # nstop is a vector of length npat, whose m-th component contains the
  # number of simulation runs in which the trial is stopped after the
  # inclusion of the m-th patient.

  # Loop over simulation runs:
  for (i in (1:nsim)) {

    t <- rep(0, npat)
    # t is the vector containing the cumulative number of observed
    # toxicities.

    t[1] <- rbinom(size = 1, prob = ptox, n = 1)
    # observation for the first patient

    counttox <- t[1]
    # counttox counts the cumulative number of observed toxicities

    # Loop over the number of accrued patients:
    for (j in (2:npat)) {

      # Binomial random draw for the tox of the next  patient:
      tox <- rbinom(size = 1, prob = ptox, n = 1)
      t[j] <- t[j - 1] + tox
      counttox <- counttox + tox

      if ((t[j] > (crit[j] - 0.5)) && (counttox > 1)) {
        nstop[j] <- nstop[j] + 1
        break  # The trial is stopped because of tox, and the next
        # simulation run is started
      }


    } # End of loop over the number of accrued patients

  }  # End of loop of the simulation runs

  # Probability of stopping at the m-th patient:
  pstop.m <- round(100*nstop/nsim, digits = 2)

  # First patient:
  crit[1] <- NA

  # Vector of cumulative probabilities of stopping:
  pstop.cum <- cumsum(pstop.m)

  # Matrix of results:
  mres<-matrix(c((1:npat), crit, pstop.m, pstop.cum) ,npat, 4)
  colnames(mres) <- c('Pat.no m', '  Crit. boundary', '  Prob. stop (%)',   '  Cumulative prob. (%)')
  print(mres)

  #' ToxCrit
  #'
  #' Calculates critical values for continuous monitoring of toxicity using a Bayesian approach. It also
  #' simulates the probability for stopping the trial for each recruited patient given a true toxicity
  #' rate passed in this function. \cr
  #' To calculate the posterior distribution of the toxicity rate r, the binomial-beta model is used.
  #' Thus, if y cases of toxicities have been observed among the first k patients, the posterior distribution
  #' of r is beta(1+y, 1+k-y). The formal stopping criterion is reached if the posterior probability (Pr) that the
  #' true toxicity rate r exceeds the unacceptable toxicity rate. \cr
  #' Please note: The trial must be stopped or modified if the observed number of toxicities is >= the critical value for the respective number of recruited patients.
  #' This method is mostly used for early phase clinical trials and validated
  #' according to GCP. It was firstly presented at the 54th annual meeting of GMDS. \cr
  #' Assumptions:\cr
  #' 1. The prior distribution of the toxicity rate is assumed to be beta(a,b).\cr
  #' 2. No stopping at the first patient.\cr
  #' IMPORTANT INFORMATION:\cr
  #' The program is closed automatically after each execution on a server.
  #' Therefore, always the same pseudo random numbers are used and simulation based on the
  #' the same parameters will always yield identical output!
  #' In order to fix this, an additional random number is generated based on a time-dependent
  #' seed point.
  #'
  #' @param npat number of patients to be included in the trial.
  #' @param ptoxcrit critical (minimal unacceptable) toxicity rate (\%).
  #' @param pstop critical posterior probability of unacceptable toxicity; higher values lead to
  #'  the termination of the trial/accrual.
  #' @param nsim number of simulation runs (simulated studies).
  #' @param ptox assumed true toxicity rate used for simulations (\%).
  #' @param a first parameter of the Beta distribution.
  #' @param b second parameter of the Beta distribution.
  #' @return A table containing information on critical boundaries per patient, the probability
  #'  to stop the trial given the assumed true toxicity rate and the cumulative probability to stop
  #'  the trial given the assumed true toxicity rate.
  #' @examples
  #'  #Pr(r>30%)>=95%) --> The posterior probability that the true toxicity rate exceeds the
  #'  #unacceptable toxicity rate of 30% is at least 95% (with 9 patients, an assumed true
  #'  #toxicity rate of 20% and 100,000 simulation runs):
  #'  ToxCrit(9, 30, 95, 100000, 20)
  #' @export
  #' @importFrom stats rbinom qbeta


}

