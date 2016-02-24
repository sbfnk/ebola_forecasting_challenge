model ebola_fit {

  const e_rho_erlang = 2

  dim rho_erlang(e_rho_erlang)

  // durations
  param p_d_onset2notification // onset to notification
  param p_d_onset2outcome // onset to outcome
  param p_d_incubation // incubation

  // rates
  param p_r_gamma // loss of infectiousness after notifiaction
  param p_r_nu // notification
  param p_r_rho // rate of incubation == 1 / (incubation period)

  // other parameters
  param p_rep // reporting rate
  param p_phi // overdispersion in reporting
  param p_vol_R0 // volatility of R0

  param p_N

  param p_init_E
  param p_init_R0

  state S // susceptible
  state E[rho_erlang] // 2 exposed compartments
  state C // infectious in the community
  state Q // infectious but notified
  state R // recovered and immune
  state Z   // accumulator
  state R0
  state startR0
  state WalkR0
  state Import_track
  state next_obs

  input import

  obs Inc

  noise n_R0_walk
  noise n_notification

  sub transition {

    Z <- (t_next_obs > next_obs ? 0 : Z)
    next_obs <- (t_next_obs > next_obs ? t_next_obs : next_obs)

    Import_track <- import

    n_notification ~ gamma(shape = 100, scale = 0.01)
    
    ode (alg = "RK4", h = 1) {

      dS/dt =
      - R0 / (p_d_onset2outcome) * (C + Q) / (S + E[0] + E[1] + C + Q + R) * S // infection

      dE[rho_erlang]/dt =
      + (rho_erlang == 0 ? R0 / (p_d_onset2outcome) * (C + Q) / (S + E[0] + E[1] + C + Q + R) * S : 0)
      + (rho_erlang > 0 ? e_rho_erlang * p_r_rho * E[rho_erlang - 1] : 0) // moving through incubation compartments
      - e_rho_erlang * p_r_rho * E[rho_erlang] // moving through incubation compartments

      dC/dt =
      + e_rho_erlang * p_r_rho * E[e_rho_erlang - 1] // becoming infectious
      - p_r_nu * C // notification
      
      dQ/dt =
      + p_r_nu * C // notification
      - p_r_gamma * Q // loss of infectiousness

      dR/dt =
      + p_r_gamma * Q  // loss of infectiousness

      dZ/dt =
      + p_r_nu * C * n_notification // accumulator for incidence
    }

    E[rho_erlang] <- E[rho_erlang] + import * p_init_E / 2 // imports
    
    n_R0_walk ~ wiener()

    // walk
    WalkR0 <- startR0 * (WalkR0 * exp(p_vol_R0 * n_R0_walk))
    // init
    WalkR0 <- (startR0 == 0 && import > 0 ? p_init_R0 : WalkR0)
    // bound
    R0 <- max(WalkR0, 0)
    // start walking
    startR0 <- min(startR0 + import + WalkR0, 1)
  }

  sub observation {
    Inc ~ truncated_normal(p_rep * Z, sqrt(max(p_rep * (1 - p_rep) * Z + (p_rep * Z * p_phi)**2, 1)), lower = 0)
  }

  sub parameter {

    p_r_gamma <- 1 / (p_d_onset2outcome - p_d_onset2notification)
    p_r_nu <- 1 / p_d_onset2notification
    p_r_rho <- 1 / p_d_incubation

    p_phi ~ uniform(lower = 0, upper = 0.5)
    p_vol_R0 ~ uniform(lower = 0, upper = 0.5)
    p_init_E ~ uniform(lower = 0, upper = 5)
    p_init_R0 ~ uniform(lower = 0, upper = 5)
  }

  sub proposal_parameter {
    p_phi ~ truncated_normal(p_phi, 1, lower = 0, upper = 0.5)
    p_vol_R0 ~ truncated_normal(p_vol_R0, 0.1, lower = 0, upper = 0.5)
    p_init_E ~ truncated_normal(p_init_E, 1, lower = 0, upper = 5)
    p_init_R0 ~ truncated_normal(p_init_R0, 0.25, lower = 0, upper = 5)
  }
  
  sub initial {
    next_obs <- 0
    Z <- 0
    R <- 0
    C <- 0
    R0 <- 0
    Q <- 0
    E[rho_erlang] <- 0
    S <- p_N
    startR0 <- 0
    Import_track <- 0
    WalkR0 <- 0
  }
}
