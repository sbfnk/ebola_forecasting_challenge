## script to convert last batch to common format

library('readr')
library('dplyr')

last_batch <- read_csv(path.expand("~/code/ebola_forecasting_challenge/data/country_level_6.csv"))

for (scenario in unique(last_batch$Scenario))
{
  country_scenario_data <- last_batch %>%
    filter(Scenario == scenario) %>%
    select(-New.Deaths, -Scenario) %>%
    rename(`# of new_confirmed_EVD_cases` = New.Confirmed.Cases)
  write_csv(country_scenario_data, path.expand(paste0("~/code/ebola_forecasting_challenge/data/sc", scenario, "_weekly_new_confirmed_EVD_cases_at_country_level_6.csv")))
}

