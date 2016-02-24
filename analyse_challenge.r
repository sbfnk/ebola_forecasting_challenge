library('RBi')
library('RBi.helpers')
library('data.table')
library('truncnorm')
library('scales')
library('binom')

time.point <- 5

code_dir <- path.expand("~/code/ebola/")
data_dir <- paste(code_dir, "data/Challenge/", sep = "/")
fit_dir <- path.expand("~/Data/Ebola/Challenge")

libbi_dir <- paste0(code_dir, "libbi/")

linelist_params <- readRDS(paste0(fit_dir, "/params_", time.point, ".rds"))
## county_cases <- cases[county == "Bomi"]

## res_fit <- bi_read_file(paste(fit_dir, "bomi", "libbi", "fit_challenge_1.nc", sep = "/"))
## rdt_fit <- lapply(res_fit, data.table)
## rdt_fit$Z <- rdt_fit$Z[nr %% 7 == 0]

## p <- plot_libbi_results(rdt_fit, data = county_cases)

predictions <- list()

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

    predictions[[scenario]] <- list()
    if (scenario == 1)
    {
        data_file <- paste0("sc1_weekly_new_confirmed_EVD_cases_at_county_level_", time.point, ".csv")
    } else if (scenario %in% seq(2, 4))
    {
        data_file <- paste0("sc", scenario,
                            "_weekly_new_confirmed_EVD_cases_at_country_level_",
                            time.point, ".csv")
    }

    cases <- data.table(read.csv(paste(data_dir, data_file , sep = "/")))
    setnames(cases, 2, "value")
    if (!("county" %in% colnames(cases)))
    {
        region.nb <- 1
        cases[, county := "Liberia"]
    } else if (!("Liberia" %in% unique(cases[, county]))) {
        cases <-
            rbind(cases,
                  data.table(county = "Liberia",
                             cases[, list(value = sum(value)), by = "week"]))
    }
    cases[, time := week * 7 + 7] ## incidence is measured at the end of the week
    cases[, state := "Inc"]
    cases[, county := gsub(" ", "", county)]

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
                      unique(cases[value > 0, county]))

    for (geo in geos)
    {
        cat("  ", geo, "\n")
        pred <- bi_read(paste(fit_dir, geo, "libbi",
                              paste0("pred_fit_challenge_", scenario, "_", time.point, ".nc"),
                              sep = "/"), vars = c("Z", "Q", "R", "R0"))
        pred <- lapply(pred, data.table)
        fit <- bi_read(paste(fit_dir, geo, "libbi",
                             paste0("fit_challenge_", scenario, "_", time.point, ".nc"),
                             sep = "/"), vars = c("Z", "Q", "R", "R0"))
        fit <- lapply(fit, data.table)
        predictions[[scenario]][[geo]] <- list()
        pred$Inc <- pred$Z[nr %% 7 == 0] ## to be set to zero
        fit$Inc <- fit$Z[nr %% 7 == 0 & nr < max(nr)] ## to be set to zero
        Inc <- rbind(fit$Inc, pred$Inc)
        ## to be deleted
        Inc[, value := value * 0.6 / 0.7]

        ## Inc[, value := rtruncnorm(nrow(Inc), 0, Inf, mean = rep * Z,
        ##                           sd = sqrt(rep * (1 - rep) * Z + (rep * Z * phi)**2))]
        ## Inc[!is.finite(value), value := 0]

        CumInc <- copy(Inc)
        CumInc <- CumInc[, value := cumsum(value), by = "np"]

        p <- plot_libbi(list(Inc = Inc, CumInc = CumInc, R0 = fit$R0))

        ## peak week
        Inc.peak <- Inc[, list(peak.value = max(value)), by = np]
        Inc.peak <- merge(Inc, Inc.peak, by = "np")
        Inc.peak <- Inc.peak[value == peak.value]
        ## Inc.peak <- merge(Inc.peak[, list(np, nr, peak.value)], rdt$time, by = "nr")
        setkey(Inc.peak, np)
        Inc.peak[, week := nr %/% 7]
        predictions[[scenario]][[geo]][["peak.week"]] <-
            list(median = median(Inc.peak[, week]),
                 low = quantile(Inc.peak[, week], 0.25),
                 high = quantile(Inc.peak[, week], 0.75))
        predictions[[scenario]][[geo]][["peak"]] <-
            list(median = median(Inc.peak[, peak.value]),
                 low = quantile(Inc.peak[, peak.value], 0.25),
                 high = quantile(Inc.peak[, peak.value], 0.75))

        unresolved <- copy(pred$Q)
        setnames(unresolved, "value", "Q")
        unresolved <- merge(unresolved, pred$R, by = c("np", "nr"))
        setnames(unresolved, "value", "R")
        unresolved[, proportion.unresolved := Q / (Q + R)]

        Inc <- merge(Inc, unresolved[, list(nr, np, proportion.unresolved)],
                     by = c("np", "nr"))
        Inc[, deaths := value * (1 - proportion.unresolved) * cfr]

        Death.peak <- Inc[, list(peak.deaths = max(deaths)), by = np]
        Death.peak <- merge(Inc, Death.peak, by = "np")
        Death.peak <- Death.peak[deaths == peak.deaths]
        setkey(Death.peak, np)
        Death.peak[, week := nr %/% 7]
        predictions[[scenario]][[geo]][["peak.deaths"]] <-
            list(median = median(Death.peak[, peak.deaths]),
                 low = quantile(Death.peak[, peak.deaths], 0.25),
                 high = quantile(Death.peak[, peak.deaths], 0.75))

        inc.min.time <- min(pred$Inc[, nr])
        for (weeks in 1:4)
        {
            predictions[[scenario]][[geo]][[paste("week", weeks, sep = ".")]] <-
                list(median = p$data$states[state == "Inc"][time == inc.min.time + 7 * weeks, value],
                     low = p$data$states[state == "Inc"][time == inc.min.time + 7 * weeks, min.1],
                     high = p$data$states[state == "Inc"][time == inc.min.time + 7 * weeks, max.1])
        }
        max_time <- p$data$states[state == "CumInc"][, max(time)]
        predictions[[scenario]][[geo]][["final.size"]] <-
            list(median = p$data$states[state == "CumInc"][time == max_time, value],
                 low = p$data$states[state == "CumInc"][time == max_time, min.1],
                 high = p$data$states[state == "CumInc"][time == max_time, max.1])
        max_time <- p$data$states[state == "R0"][, max(time)]
        predictions[[scenario]][[geo]][["R0"]] <-
            list(median = p$data$states[state == "R0"][time == max_time, value],
                 low = p$data$states[state == "R0"][time == max_time, min.1],
                 high = p$data$states[state == "R0"][time == max_time, max.1])
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

