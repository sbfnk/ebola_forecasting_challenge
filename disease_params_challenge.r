library('data.table')

fit_dir <- path.expand("~/Data/Ebola/Challenge")
code_dir <- path.expand("~/code/ebola_forecasting_challenge/")
data_dir <- paste(code_dir, "data/", sep = "/")

params <- list()

for (time.point in 1:5)
{
  for (scenario in 1:4)
  {
    params[[scenario]] <- list()
    if (scenario %in% c(1, 3, 4))
    {
      ## disease parameters
      linelist <- NULL
      linelist.file <- paste0(data_dir, "/sc", scenario, "_patient_records_", time.point, ".csv")
      if (file.exists(linelist.file))
      {
        linelist <- data.table(read.csv(linelist.file))
      } else 
      {
        ttree_files <- list.files(path = data_dir,
                                  pattern = paste0("scenario0", scenario, "_transmission_tree_.*.csv"),
                                  full.names = TRUE)
        ttrees <- list()
        for (file in ttree_files)
        {
          ttrees[[file]] <- data.table(read.csv(file))
        }
        if (length(ttrees) > 0)
        {
          if (is.null(linelist))
          {
            linelist <- rbindlist(ttrees)
          } else {
            linelist <- rbind(linelist, rbindlist(ttrees), fill = TRUE)
          }
        }
      }

      if (!is.null(linelist))
      {
        for (time_column in grep("\\.Time$", colnames(linelist), value = TRUE)) {
          linelist[, paste(time_column) := as.numeric(as.character(get(time_column)))]
        }
        params[[scenario]][["died"]] <- as.integer(linelist[, !is.na(Death.Time)])
        params[[scenario]][["onset2death"]] <- linelist[!is.na(Death.Time) & !is.na(Symptom.Time), Death.Time - Symptom.Time]
        params[[scenario]][["incubation"]] <- linelist[!is.na(Infection.Time) & !is.na(Symptom.Time), Symptom.Time - Infection.Time]
        params[[scenario]][["generation"]] <- linelist[!is.na(Death.Time) & !is.na(Infection.Time), Death.Time - Infection.Time]
      }
    }
  }
  saveRDS(params, paste0(fit_dir, "/params_", time.point, ".rds"))
}

