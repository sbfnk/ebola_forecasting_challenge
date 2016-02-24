library('data.table')

fit_dir <- path.expand("~/Data/Ebola/Challenge")
code_dir <- path.expand("~/code/ebola/")
data_dir <- paste(code_dir, "data/Challenge/", sep = "/")

time.point <- 5

params <- list()

for (scenario in 1:4)
{
    params[[scenario]] <- list()
    if (scenario %in% c(1, 3, 4))
    {
        ## disease parameters
        linelist <- NULL
        if (scenario %in% c(1, 3))
        {
            linelist <- data.table(read.csv(paste0(data_dir, "/sc", scenario, "_patient_records_", time.point, ".csv")))
        } else if (scenario %in% c(3, 4))
        {
            ttree_files <- list.files(path = data_dir,
                                      pattern = paste0("scenario0", scenario, "_transmission_tree_.*.csv"),
                                      full.names = TRUE)
            ttrees <- list()
            for (file in ttree_files)
            {
                ttrees[[file]] <- data.table(read.csv(file))
            }
            if (is.null(linelist))
            {
                linelist <- rbindlist(ttrees)
            } else {
                linelist <- rbind(linelist, rbindlist(ttrees), fill = TRUE)
            }
        }
        for (time_column in grep("\\.Time$", colnames(linelist), value = TRUE)) {
            linelist[, paste(time_column) := as.numeric(as.character(get(time_column)))]
        }
        params[[scenario]][["died"]] <- as.integer(linelist[, !is.na(Death.Time)])
        params[[scenario]][["onset2outcome"]] <- linelist[!is.na(Death.Time) & !is.na(Symptom.Time), Death.Time - Symptom.Time]
        params[[scenario]][["incubation"]] <- linelist[!is.na(Infection.Time) & !is.na(Symptom.Time), Symptom.Time - Infection.Time]
        params[[scenario]][["generation"]] <- linelist[!is.na(Death.Time) & !is.na(Infection.Time), Death.Time - Infection.Time]
    }
}

saveRDS(params, paste0(fit_dir, "/params_", time.point, ".rds"))
