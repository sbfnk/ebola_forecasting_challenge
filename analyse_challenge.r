library('RBi')
library('RBi.helpers')
library('data.table')
library('truncnorm')
library('scales')
library('binom')
library('cowplot')

time.points <- c(1, 2, 4, 5)

data.time.points <- data.table(point = 1:5, time = c(13, 20, NA, 35, 42))

code_dir <- path.expand("~/code/ebola_cmmid/")
data_dir <- paste(code_dir, "data/Challenge/", sep = "/")
fit_dir <- path.expand("~/Data/Ebola/Challenge")

libbi_dir <- paste0(code_dir, "libbi/")

## county_cases <- cases[county == "Bomi"]

## res_fit <- bi_read_file(paste(fit_dir, "bomi", "libbi", "fit_challenge_1.nc", sep = "/"))
## rdt_fit <- lapply(res_fit, data.table)
## rdt_fit$Z <- rdt_fit$Z[nr %% 7 == 0]

## p <- plot_libbi_results(rdt_fit, data = county_cases)

predictions <- list()
full_predictions <- list()
r0_trajectories <- list()
param_data <- list()
cases <- list()

for (time.point in time.points)
{
  tp <- paste("time", time.point, sep = ".")
  predictions[[tp]] <- list()
  full_predictions[[tp]] <- list()
  r0_trajectories[[tp]] <- list()
  param_data[[tp]] <- list()
  cases[[tp]] <- list()

  linelist_params <- readRDS(paste0(fit_dir, "/params_", time.point, ".rds"))
  for (scenario in 1:4)
  {
    cat("Scenario ", scenario, "\n")
    if ("died" %in% names(linelist_params[[scenario]]))
    {
      cfr <- sum(linelist_params[[scenario]][["died"]]) /
        length(linelist_params[[scenario]][["died"]])
    } else
    {
      cfr <- 0.88
    }

    predictions[[tp]][[scenario]] <- list()
    full_predictions[[tp]][[scenario]] <- list()
    r0_trajectories[[tp]][[scenario]] <- list()
    cases[[tp]][[scenario]] <- list()
    param_data[[tp]][[scenario]] <- list()
    if (scenario == 1)
    {
      data_file <- paste0("sc1_weekly_new_confirmed_EVD_cases_at_county_level_", time.point, ".csv")
    } else if (scenario %in% seq(2, 4))
    {
      data_file <- paste0("sc", scenario,
                          "_weekly_new_confirmed_EVD_cases_at_country_level_",
                          time.point, ".csv")
    }

    cases[[tp]][[scenario]] <- data.table(read.csv(paste(data_dir, data_file , sep = "/")))
    setnames(cases[[tp]][[scenario]], 2, "value")
    if (!("county" %in% colnames(cases[[tp]][[scenario]])))
    {
      region.nb <- 1
      cases[[tp]][[scenario]][, county := "Liberia"]
    } else if (!("Liberia" %in% unique(cases[[tp]][[scenario]][, county]))) {
      cases[[tp]][[scenario]] <-
        rbind(cases[[tp]][[scenario]],
              data.table(county = "Liberia",
                         cases[[tp]][[scenario]][, list(value = sum(value)), by = "week"]))
    }
    cases[[tp]][[scenario]][, time := week * 7 + 7] ## incidence is measured at the end of the week
    cases[[tp]][[scenario]][, state := "Inc"]
    cases[[tp]][[scenario]][, county := gsub(" ", "", county)]
    cases[[tp]][[scenario]] <- cases[[tp]][[scenario]][week <= data.time.points[point == time.point, time]]

    demographics <- data.table(read.csv(paste(data_dir, "demographics.csv", sep = "/")))
    demographics <-
      rbind(demographics,
            data.table(county.code = max(demographics[, county.code]) + 1,
                       county.name = "Liberia",
                       county.population = sum(demographics[, county.population]),
                       county.capital.population =
                         demographics[county.name == "Montserrado",
                                      county.capital.population]))

    geos <- intersect(unique(demographics[, county.name]),
                      unique(cases[[tp]][[scenario]][value > 0, county]))

    for (geo in geos)
    {
      cat("  ", geo, "\n")
      model <- bi_model(paste(fit_dir, geo,
                              paste0("fit_challenge_", scenario, "_", time.point, ".bi"),
                              sep = "/"))
      pred <- bi_read(paste(fit_dir, geo,
                            paste0("pred_fit_challenge_", scenario, "_", time.point, ".nc"),
                            sep = "/"),
                      vars = c("Z", "R0", "Q", "R"), thin = 10)
      pred <- lapply(pred, data.table)
      fit <- bi_read(paste(fit_dir, geo,
                           paste0("fit_challenge_", scenario, "_", time.point, ".nc"),
                           sep = "/"),
                     vars = model$get_vars("state"), thin = 10)
      fit <- lapply(fit, data.table)
      params <- bi_read(paste(fit_dir, geo,
                              paste0("fit_challenge_", scenario, "_", time.point, ".nc"),
                              sep = "/"),
                        vars = model$get_vars("param"), thin = 10)
      params <- lapply(params, data.table)

      p <- plot_libbi(params, model)

      param_data[[tp]][[scenario]][[geo]] <- p$data$params[distribution == "posterior" & parameter %in% c("p_phi", "p_vol_R0", "p_init_E", "p_init_R0")]
      param_data[[tp]][[scenario]][[geo]][parameter == "p_phi", parameter := "phi"]
      param_data[[tp]][[scenario]][[geo]][parameter == "p_vol_R0", parameter := "sigma"]
      param_data[[tp]][[scenario]][[geo]][parameter == "p_init_E", parameter := "paste(E, \"*\")"]
      param_data[[tp]][[scenario]][[geo]][parameter == "p_init_R0", parameter := "paste(R[0], \"*\")"]

      ## p <- ggplot(param_data, aes(x = value)) +
      ##   geom_density(adjust = 2, fill = "black") +
      ##   facet_wrap(~ parameter, scales = "free", labeller = label_parsed) 

      ## ggsave(paste0("challenge_parameters_", geo, "_", time.point, "_", scenario, ".pdf"), p, height = 7, width = 7)

      predictions[[tp]][[scenario]][[geo]] <- list()
      pred$Inc <- pred$Z[time %% 7 == 0] ## to be set to zero
      fit$Inc <- fit$Z[time %% 7 == 0 & time < max(time)] ## to be set to zero
      Inc <- rbind(fit$Inc, pred$Inc)
      Inc <- merge(Inc, params$p_rep[, list(np, rep = value)], by = "np")
      Inc <- merge(Inc, params$p_phi[, list(np, phi = value)], by = "np")
      setnames(Inc, "value", "Z")
      ## to be deleted
      ## Inc[, value := value * 0.6 / 0.7]
      Inc[, mean := rep * Z]
      Inc[, sd := sqrt(rep * (1 - rep) * Z + (rep * Z * phi)**2)]
      Inc[sd < 1, sd := 1]
      Inc[, value := rnorm(n = .N, mean = mean, sd = sd)]
      Inc[!is.finite(value), value := 0]
      full_predictions[[tp]][[scenario]][[geo]] <- Inc
      r0_trajectories[[tp]][[scenario]][[geo]] <- fit$R0
      CumInc <- copy(Inc)
      CumInc <- CumInc[, value := cumsum(value), by = "np"]

      p <- plot_libbi(list(Inc = Inc, Z = CumInc, R0 = fit$R0), model = model)
      ## p_fit_pred <- plot_libbi(list(Inc = Inc), model = model, quantiles = 0.5, d, data = cases %>% filter(county == geo), limit.to.data = FALSE)ata = cases %>% filter(county == geo), limit.to.data = FALSE)
      ## p_fit_r0 <- plot_libbi(list(R0 = fit$R0), model = model, quantiles = 0.5)

      ## peak week
      Inc.peak <- Inc[, list(peak.value = max(value)), by = np]
      Inc.peak <- merge(Inc, Inc.peak, by = "np")
      Inc.peak <- Inc.peak[value == peak.value]
      ## Inc.peak <- merge(Inc.peak[, list(np, nr, peak.value)], rdt$time, by = "time")
      setkey(Inc.peak, np)
      Inc.peak[, week := time %/% 7]
      predictions[[tp]][[scenario]][[geo]][["peak.week"]] <-
        list(median = median(Inc.peak[, week]),
             low = quantile(Inc.peak[, week], 0.25),
             high = quantile(Inc.peak[, week], 0.75))
      predictions[[tp]][[scenario]][[geo]][["peak"]] <-
        list(median = median(Inc.peak[, peak.value]),
             low = quantile(Inc.peak[, peak.value], 0.25),
             high = quantile(Inc.peak[, peak.value], 0.75))

      unresolved <- copy(pred$Q)
      setnames(unresolved, "value", "Q")
      unresolved <- merge(unresolved, pred$R, by = c("np", "time"))
      setnames(unresolved, "value", "R")
      unresolved[, proportion.unresolved := Q / (Q + R)]

      Inc <- merge(Inc, unresolved[, list(time, np, proportion.unresolved)],
                   by = c("np", "time"))
      Inc[, deaths := value * (1 - proportion.unresolved) * cfr]

      Death.peak <- Inc[, list(peak.deaths = max(deaths)), by = np]
      Death.peak <- merge(Inc, Death.peak, by = "np")
      Death.peak <- Death.peak[deaths == peak.deaths]
      setkey(Death.peak, np)
      Death.peak[, week := time %/% 7]
      predictions[[tp]][[scenario]][[geo]][["peak.deaths"]] <-
        list(median = median(Death.peak[, peak.deaths]),
             low = quantile(Death.peak[, peak.deaths], 0.25),
             high = quantile(Death.peak[, peak.deaths], 0.75))

      inc.min.time <- min(pred$Inc[, time])
      for (weeks in 1:4)
      {
        predictions[[tp]][[scenario]][[geo]][[paste("week", weeks, sep = ".")]] <-
          list(median = p$data$states[state == "Inc"][time == inc.min.time + 7 * weeks, value],
               low = p$data$states[state == "Inc"][time == inc.min.time + 7 * weeks, min.1],
               high = p$data$states[state == "Inc"][time == inc.min.time + 7 * weeks, max.1])
      }
      max_time <- p$data$states[state == "Z"][, max(time)]
      predictions[[tp]][[scenario]][[geo]][["final.size"]] <-
        list(median = p$data$states[state == "Z"][time == max_time, value],
             low = p$data$states[state == "Z"][time == max_time, min.1],
             high = p$data$states[state == "Z"][time == max_time, max.1])
      max_time <- p$data$states[state == "R0"][, max(time)]
      predictions[[tp]][[scenario]][[geo]][["R0"]] <-
        list(median = p$data$states[state == "R0"][time == max_time, value],
             low = p$data$states[state == "R0"][time == max_time, min.1],
             high = p$data$states[state == "R0"][time == max_time, max.1])
    }
  }
}

saveRDS(predictions, "predictions.rds")

## get peak timing
peak_timing <- matrix(unlist(lapply(predictions, function(x)
{
    x[["Liberia"]][["peak.week"]]
})), nrow = 4, ncol = 3, byrow = TRUE)
write.table(peak_timing, file = "peak_timing.csv", row.names = FALSE,
            col.names = FALSE, sep = ",")

## get prediction time points
counties <- setdiff(demographics[, county.name], "Liberia")
counties <- counties[order(counties)]
counties <- c("Liberia", counties)

incidence_count <- rbindlist(lapply(predictions, function(x)
{
    data.frame(t(sapply(counties, function(y)
    {
        vec <- c()
        for (var in c("peak", "peak.deaths", "week.1", "week.2", "week.3", "week.4",
                      "final.size"))
        {
            if (!(y %in% names(x)))
            {
                vec <- c(vec, NA, NA, NA)
            } else
            {
                value <- round(unname(unlist((x[[y]][[var]]))))
                vec <- c(vec, value)
            }
        }
        return(vec)
    })))
}))

write.table(incidence_count, file = "incidence_count.csv",
            row.names = FALSE, col.names = FALSE,
            sep = ",", na = "")

R0 <- matrix(unlist(lapply(predictions, function(x)
{
    x[["Liberia"]][["R0"]]
})), nrow = 4, ncol = 3, byrow = TRUE)

gen_time <- matrix(unlist(lapply(linelist_params, function(x)
{
    quantile(x[["generation"]], c(0.5, 0.25, 0.75))
})), nrow = 4, ncol = 3, byrow = TRUE)

cfr <- matrix(unlist(lapply(linelist_params, function(x)
{
    unname(unlist(binom.confint(sum(x[["died"]]), length(x[["died"]]), method = "lrt")[c("mean", "lower", "upper")]))
})), nrow = 4, ncol = 3, byrow = TRUE)
cfr[2, ] <- c(NA, NA, NA)

reporting <- matrix(nrow = 4, ncol = 3)

disease_params <- cbind(R0, gen_time, cfr, reporting)

write.table(disease_params, file = "disease_param.csv", row.names = FALSE,
            col.names = FALSE, sep = ",")

model_data <- list()
for (time.point in time.points)
{
    tp <- paste("time", time.point, sep = ".")

    for (name in names(full_predictions[[tp]][[1]]))
    {
      full_predictions[[tp]][[1]][[name]][, county := name]
    }
    predictions_point <- rbindlist(full_predictions[[tp]][[1]])

    flawed <- unique(predictions_point[value > 1e+10, list(np, county)])
    flawed[, flawed := TRUE]

    predictions_point <- merge(predictions_point, flawed, all.x = TRUE, by = c("np", "county"))
    predictions_point <- predictions_point[is.na(flawed)]
    predictions_point[, flawed := NULL]

    max.cases <- cases[["time.5"]][[1]][, list(max.value = max(value)), by = county]
    keep_counties <- max.cases[max.value > 5]
    plot.cases <- cases[["time.5"]][[1]][county %in% keep_counties[, county]]

    agg <- predictions_point[county %in% keep_counties[, county], list(median = median(value), min.1 = quantile(value, 0.25), max.1 = quantile(value, 0.75), min.2 = quantile(value, 0.025), max.2 = quantile(value, 0.975)), by = c("county", "time")]
    agg[, week := time / 7]

    agg <- merge(agg, max.cases, by = "county", all.x = TRUE)
    agg[max.2 > max.value * 2, max.2 := max.value * 2]
    agg[max.1 > max.value * 2, max.1 := max.value * 2]
    agg[median > max.value * 2, median := NA]
    agg[min.1 > max.value * 2, min.1 := max.value * 2]
    agg[min.2 > max.value * 2, min.2 := max.value * 2]
    agg[min.1 < 0, min.1 := 0]
    agg[min.2 < 0, min.2 := 0]
    agg[median < 0, median := 0]

    max.data <- data.time.points[point == time.point, time]

    plot.cases[, type := "fit"]

    county_levels <- c(setdiff(unique(plot.cases[, county]), "Liberia"), "Liberia")

    agg[, county := factor(county, levels = county_levels)]
    plot.cases[, county := factor(county, levels = county_levels)]

    plot.cases[week <= max.data, forecast := FALSE]
    plot.cases[is.na(forecast), forecast := TRUE]

    fit <- agg[week <= max.data]
    fit[, type := "fit"]
    pred <- agg[week > max.data]
    pred[, type := "forecast"]

    for (horizon in c(7, 52))
    {
      p <- ggplot(fit[week > max.data - horizon & week < max.data + horizon + 2], aes(x = week, y = median)) +
        geom_line() +
        geom_ribbon(aes(ymin = min.1, ymax = max.1), alpha = 0.5) +
        geom_ribbon(aes(ymin = min.2, ymax = max.2), alpha = 0.25) +
        facet_wrap(~ county, scales = "free") +
        geom_point(data = plot.cases[week > max.data - horizon & week < max.data + horizon + 2 & forecast == FALSE], aes(y = value), color = "red") +
        geom_point(data = plot.cases[week > max.data - horizon & week < max.data + horizon + 2 & forecast == TRUE], aes(y = value), color = "black") +
        scale_y_continuous("Incidence") +
        geom_line(data = pred[week > max.data - horizon & week < max.data + horizon + 2], color = "blue") +
        geom_ribbon(data = pred[week > max.data - horizon & week < max.data + horizon + 2], aes(ymin = min.1, ymax = max.1), alpha = 0.5, fill = "blue") +
        geom_ribbon(data = pred[week > max.data - horizon & week < max.data + horizon + 2], aes(ymin = min.2, ymax = max.2), alpha = 0.25, fill = "blue")

      ggsave(paste0("challenge_incidence_", time.point, "_", horizon, ".pdf"), p, height = 10, width = 10)

      for (name in names(r0_trajectories[[tp]][[1]]))
      {
        r0_trajectories[[tp]][[1]][[name]] <- r0_trajectories[[tp]][[1]][[name]][value > 0]
        r0_trajectories[[tp]][[1]][[name]][, county := name]
      }
      trajectories <- rbindlist(r0_trajectories[[tp]][[1]])

      agg.r0 <- trajectories[county %in% keep_counties[, county], list(median = median(value), min.1 = quantile(value, 0.25), max.1 = quantile(value, 0.75), min.2 = quantile(value, 0.025), max.2 = quantile(value, 0.975)), by = c("county", "time")]
      agg.r0[, week := time / 7]
      agg.r0[, county := factor(county, levels = county_levels)]
      agg.r0 <- agg.r0[week <= max.data]
      agg.r0 <- agg.r0[week > max.data - horizon]

      agg.pred <- rbind(agg.r0[week == max.data][, week := week + 1],
                        agg.r0[week == max.data][, week := week + horizon + 1])

      p <- ggplot(agg.r0, aes(x = week, y = median)) +
        geom_line() +
        geom_ribbon(aes(ymin = min.1, ymax = max.1), alpha = 0.5) +
        geom_ribbon(aes(ymin = min.2, ymax = max.2), alpha = 0.25) +
        facet_wrap(~ county, scales = "free") +
        geom_line(data = agg.pred, color = "blue") +
        geom_ribbon(data = agg.pred, aes(ymin = min.1, ymax = max.1), alpha = 0.5, fill = "blue") +
        geom_ribbon(data = agg.pred, aes(ymin = min.2, ymax = max.2), alpha = 0.25, fill = "blue") +
        facet_wrap(~ county, scales = "free") +
        scale_y_continuous(expression(R[0])) +
        geom_hline(yintercept = 1, linetype = "dashed")

      ggsave(paste0("challenge_r0_", time.point, "_", horizon,  ".pdf"), p, height = 10, width = 10)
    }

    model_data[[tp]] <- data.table(merge(pred, cases$time.5[[1]], by = c("week", "county")))
    model_data[[tp]][, add.weeks := week - min(week) + 1]
    model_data[[tp]][, point := time.point]

    for (name in names(param_data[[tp]][[1]]))
    {
      param_data[[tp]][[1]][[name]][, county := name]
    }
    params_point <- rbindlist(param_data[[tp]][[1]])
    params_point <- params_point[county %in% keep_counties[, county]]
    params_point[, county := factor(county, levels = county_levels)]
    params_point_summary <- params_point[, list(median = median(value),
                                                min = quantile(value, 0.25),
                                                max = quantile(value, 0.75)),
                                         by = list(parameter, distribution, varying, county)]

    p <- ggplot(params_point[parameter %in% c("phi", "sigma")],
                ## aes(x = county, y = mean, ymin = min, ymax = max, color = parameter)) +
                aes(x = county, y = value, color = parameter)) +
      geom_boxplot() +
      scale_color_brewer("", palette = "Set1", labels = c(expression(phi), expression(sigma))) +
      coord_flip() +
      scale_x_discrete("") +
      scale_y_continuous("") +
      theme(legend.position = "top")
    ## geom_point() +
    ## geom_errorbar()
    ggsave(paste0("challenge_error_", time.point, ".pdf"), p, heigh = 7, width = 7)
    save_plot(paste0("challenge_error_", time.point, ".pdf"), p)
}

compare <- rbindlist(model_data)

pred <- apply(unique(compare[, list(add.weeks, point)]), 1, function(x)
{
  data.frame(inside.50 = compare[add.weeks == x[["add.weeks"]] & point == x[["point"]], sum(min.1 < value & max.1 > value) / .N], 
            inside.95 = compare[add.weeks == x[["add.weeks"]] & point == x[["point"]], sum(min.2 < value & max.2 > value) / .N],
            greater.median = compare[add.weeks == x[["add.weeks"]] & point == x[["point"]], sum(median < value, na.rm = TRUE) / .N], 
            weeks = x[["add.weeks"]], time.point = x[["point"]])
})

pred <- rbindlist(pred)

mp <- melt(data.table(pred), id.vars = c("weeks", "time.point"))
mp$weeks <- unlist(mp$weeks)
mp$time.point <- unlist(mp$time.point)
mp$value <- unlist(mp$value)
mp[variable == "inside.50", variable := "inside 50% CI"]
mp[variable == "inside.95", variable := "inside 95% CI"]
mp[variable == "greater.median", variable := "greater than median"]

ideal <- data.frame(variable = unique(mp$variable),
                    value = c(0.5, 0.95, 0.5))

p <- ggplot(mp) +
  geom_point(aes(x = weeks, y = value, shape = factor(time.point))) +
  geom_hline(data = ideal, aes(yintercept = value)) +
  facet_wrap(~ variable, ncol = 3) +
  scale_y_continuous("proportion of observations", limits = c(0, 1))
ggsave("challenge_forecasting_performance.pdf", p, height = 3, width = 7)

