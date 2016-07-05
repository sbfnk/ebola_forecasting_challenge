weeks_import <- 3 # number of Ebola-free weeks after which a new cases is treated as an import

onset2notification <- 7 ## assume a week delay in reporting

prediction_window <- 2 * 365 ## number of days to predict

code_dir <- path.expand("~/code/ebola_forecasting_challenge/")
data_dir <- paste(code_dir, "data/", sep = "/")
fit_dir <- path.expand("~/Data/Ebola/Challenge")

############################################################################
## Read command line options                                              ##
############################################################################

library('docopt')

"Script for fitting the Ebola model to data in libbi.

Usage: fit_regions.r [options]

Options:
  -r --region=<region.nb>        region number
  -i --timepoint=<time.point>    time point
  -n --samples=<num.sims>        number of samples
  -p --presamples=<pre.samples> number of preparatory samples to obtain
  -t --threads=<num.threads>     number of threads
  -d --dont-run                  don't run the model, just write the file(s)
  -s --scenario=<scenario.nb>    scenario to consider
  -o --output=<output>           output file name
  -h --help                      show this message" -> doc

opts <- docopt(doc)

if (opts[["help"]])
{
    print(opts)
    exit()
}

region.nb <- as.integer(opts[["region"]])
num.sims <- as.integer(opts[["samples"]])
pre.sims <- as.integer(opts[["presamples"]])
num.threads <- as.integer(opts[["threads"]])
time.point <- as.integer(opts[["timepoint"]])
scenario <- as.integer(opts[["scenario"]])
runModel <- !(opts[["dont-run"]])
output <- opts[["output"]]

if (length(time.point) == 0)
{
    stop("Must give time point using --timepoint")
}

if (length(num.sims) == 0)
{
    warning("Number of samples not given, setting to 10")
    num.sims <- 10
}

if (length(pre.sims) == 0)
{
    warning("Number of presamples not given, setting to 10")
    pre.sims <- 10
}

if (length(num.threads) == 0)
{
    num.threads <- 1
}

if (length(scenario) == 0)
{
    stop("Must provide scenario")
}

if (is.null(output))
{
    output <- paste0("fit_challenge_", scenario, "_", time.point)
}

if (!any(grep("\\.nc", output)))
{
    output <- paste(output, "nc", sep = ".")
}

library('RBi')
library('RBi.helpers')
library('data.table')

date()

if (scenario == 1) {
    data_file <- paste0("sc1_weekly_new_confirmed_EVD_cases_at_county_level_", time.point, ".csv")
} else if (scenario %in% seq(2, 4)) {
    data_file <- paste0("sc", scenario,
                        "_weekly_new_confirmed_EVD_cases_at_country_level_", time.point, ".csv")
}

cases <- data.table(read.csv(paste(data_dir, data_file , sep = "/")))
setnames(cases, 2, "cases")
if ("county" %in% colnames(cases))
{
    cases[, county := gsub(" ", "", county)]
}

demographics <- data.table(read.csv(paste(data_dir, "demographics.csv", sep = "/")))
demographics <-
    rbind(demographics,
          data.table(county.code = max(demographics[, county.code]) + 1,
                     county.name = "Liberia",
                     county.population = sum(demographics[, county.population]),
                     county.capital.population =
                         demographics[county.name == "Montserrado",
                                      county.capital.population]))

if (!("county" %in% colnames(cases)))
{
    if (length(region.nb) > 0)
    {
        warning("Region number given but national data. Will do a national fit")
    }
    region.nb <- 1
    cases[, county := "Liberia"]
} else if (!("Liberia" %in% unique(cases[, county]))) {
    cases <-
        rbind(cases,
              data.table(county = "Liberia",
                         cases[, list(cases = sum(cases)), by = "week"]))
}
cases[, county := gsub(" ", "", county)]

geos <- intersect(unique(demographics[, county.name]), unique(cases[cases > 0, county]))

if (length(region.nb) == 0) {
    print("No region given, considering all regions")
    region.nb <- seq_along(geos)
}

## mark new cases after weeks_import weeks as imports
setkey(cases, county, week)

min.weeks <- cases[cases > 0, list(min.week = min(week)), by = county]
cases <- merge(cases, min.weeks, by = "county")
max.weeks <- cases[cases > 0, list(max.week = max(week)), by = county]
cases <- merge(cases, max.weeks, by = "county")

for (i in seq_len(weeks_import))
{
    cases[, paste("prev", i, sep = ".") := 0]
    previous.cases <- cases[, list(prev = head(cases, -i)), by = county][, prev]
    cases[week > i, paste("prev", i, sep = ".") := as.numeric(previous.cases)]
}
prev.cases <- rowSums(cases[, paste("prev", seq_len(weeks_import), sep = "."), with = F], na.rm = TRUE)
cases[, prev.sum := prev.cases]
cases[, import := 0]
cases[week == min.week, import := 1]
cases[, min.week := NULL]
cases[, max.week := NULL]
cases[, paste("prev", seq_len(weeks_import), sep = ".") := NULL]

cases[(is.na(prev.sum) | prev.sum == 0) & cases > 0, import := 1]
cases[, prev.sum := NULL]

## preparation for libbi
cases[, obs.week := week + 1] # incidence is measured at the end of a week
cases[, import.week := week - 1] # imports happen a week before

## daily imports
import <- copy(cases[, list(county, week, import, import.week)])
import[, import.time := import.week * 7]
all_times <-
    data.table(expand.grid(import.time = seq_len(max(import[, import.time]) + 1),
                           county = geos))
import <- merge(import, all_times, by = c("import.time", "county"), all = TRUE)
import[is.na(import), import := 0]

## disease parameters
params  <- list()
linelist_params <- readRDS(paste0(fit_dir, "/params_", time.point, ".rds"))

if ("onset2outcome" %in% names(linelist_params[[scenario]]))
{
    params[["p_d_onset2outcome"]] <- mean(linelist_params[[scenario]][["onset2outcome"]])
} else
{
    cfr <- 0.88
    params[["p_d_onset2outcome"]] <- 7.49 * cfr + 10 * (1 - cfr) # from Epidemics paper
}

if ("incubation" %in% names(linelist_params[[scenario]]))
{
    params[["p_d_incubation"]] <- mean(linelist_params[[scenario]][["incubation"]])
} else
{
    params[["p_d_incubation"]] <- 5.99
}

params[["p_d_onset2notification"]] <- 7 # guessed
params[["p_rep"]] <- 0.7 # from camacho et al
params[["p_first_obs"]] <- 7 # first observation

for (nb in region.nb) {
    region <- geos[nb]
    print(region)

    region.lower <- tolower(gsub(" ", "_", region))
    dir.region.lower <- paste0(fit_dir, "/Regions/", region.lower, "/")
    dir.region.lower.libbi <- paste0(dir.region.lower, "libbi/")
    dir.region.lower.libbi.work <- paste0(dir.region.lower, "libbi/work_",
                                          scenario, "_", time.point)

    ## create directory for input files and libbi control files
    suppressWarnings(dir.create(dir.region.lower))
    suppressWarnings(dir.create(dir.region.lower.libbi))
    unlink(dir.region.lower.libbi.work, recursive = TRUE)
    suppressWarnings(dir.create(dir.region.lower.libbi.work))
    ## set region-specific parameters
    params[["p_N"]] <- demographics[county.name == region, county.population]

    ############################################################################
    ## Generate input and observation files                                   ##
    ############################################################################

    ## generate init file
    init_file_name <- paste0(dir.region.lower.libbi, "/init_", output)
    output_file_name <- paste0(dir.region.lower.libbi, "/", output)
    input_file_name <- paste0(dir.region.lower.libbi, "/input_", output)

    bi_write(filename = init_file_name, variables = params)
    ## check
    ## bi_file_ncdump(init_file_name)

    times <- import[county == region, import.time]
    ## generate input file
    time_var <- data.table(nr = times, value = times)
    import_var <-
        data.table(nr = times,
                   value = import[county == region, import])
    bi_write(variables = list(time_import = time_var, import = import_var),
             filename = input_file_name)
    ## check
    ## bi_file_ncdump(input_file_name)

    ## generate observation file
    obs_file_name <- paste0(dir.region.lower.libbi, "obs_", output)
    incidence <- cases[county == region, list(day = obs.week * 7, value = cases)]
    bi_write(filename = obs_file_name,
             variables = list(Inc = incidence),
             time_dim = "day")
    ## check
    ## bi_file_ncdump(obs_file_name)

    end_time <- max(cases[county == region, obs.week * 7])

    global_options <-
        list("obs-file" = obs_file_name,
             "init-file" = init_file_name,
             "input-file" = input_file_name,
             assert = FALSE,
             "end-time" = end_time,
             "start-time" = 0,
             noutputs = end_time,
             nsamples = pre.sims,
             nthreads = num.threads)

    if (num.threads == 1)
    {
        global_options <- c(global_options, list(openmp = FALSE))
    }

    ############################################################################
    ## fit the model                                                          ##
    ############################################################################

    ebola_model <- bi_model(paste(code_dir, "ebola_fit_challenge.bi",
                                  sep = "/"))
    ebola_deter <- ebola_model$clone()
    ebola_deter$fix(n_R0_walk = 0, n_notification = 1)
    ebola_deter_prior <- ebola_deter$propose_prior()

    ## run deterministic model,  sample from prior
    run_det_prior <- libbi(client = "sample", model = ebola_deter_prior,
                           global_options = global_options, run = TRUE,
                           working_folder = dir.region.lower.libbi.work)

    ## check
    ## res <- bi_read(run_det_prior$output_file_name)

    ## adapt deterministic model
    run_det_adapted <- adapt_mcmc(run_det_prior, min = 0.1, max = 0.5, max_iter = 10,
                                  scale = 2)

    ## adapt stochastic model
    ebola_model$add_block("proposal_parameter",
                          run_det_adapted$model$get_block("proposal_parameter"))
    ebola_model$add_block("proposal_initial",
                          run_det_adapted$model$get_block("proposal_initial"))

    run <- libbi(client = "sample", model = ebola_model,
                 global_options = global_options, run = TRUE,
                 working_folder = dir.region.lower.libbi.work)

    run_particle_adapted <- adapt_particles(run, min = nrow(incidence))

    run_adapted <- adapt_mcmc(run_particle_adapted, min = 0.1, max = 0.5,
                              max_iter = 10, scale = 2)

    run_adapted$run(add_options = list("init-file" = run_adapted$output_file_name,
                                       "init-np" = pre.sims - 1, 
                                       nsamples = num.sims),
                    output_file_name = output_file_name)

    run_adapted$model$write_model_file(paste0(sub("\\.nc$", ".bi", output_file_name)))

    ############################################################################
    ## predict
    ############################################################################

    prediction_options <- global_options
    prediction_options[["target"]] <- "prediction"
    prediction_options[["start-time"]] <- end_time
    prediction_options[["end-time"]] <- end_time + prediction_window
    prediction_options[["noutputs"]] <- prediction_window
    prediction_options[["init-file"]] <- output_file_name

    prediction_file_name <- paste0(dir.region.lower.libbi, "pred_", output)

    prediction_model <- ebola_model$clone()
    prediction_model$fix(n_R0_walk = 0)
    prediction_model$remove_block("initial")

    lines <- prediction_model$get_lines()
    Z_line <- min(grep("^[[:space:]]*Z[[:space:]]*<-", lines))
    new_Z_line <- c("Z <- (mod(t_now, 7) == 0 ? 0 : Z)")
    prediction_model$update_lines(Z_line, new_Z_line)
    prediction_model$remove_lines(Z_line + 1)

    pred <- libbi(client = "sample", model = prediction_model,
                  global_options = prediction_options, run = TRUE,
                  working_folder = dir.region.lower.libbi.work,
                  output_file_name = prediction_file_name)

    sample_observations(pred,
                        output_file_name = sub("\\.nc$", "_joint.nc",
                                               prediction_file_name))
}

date()
