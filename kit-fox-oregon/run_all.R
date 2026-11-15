# Run the full pipeline from the project root: Rscript run_all.R
stages <- file.path("pipeline", c("01_occurrences.R", "02_env_layers.R", "03_model_data.R", "04_fit_models.R",
                                  "05_project.R", "06_habitat_viability.R", "07_export_web.R"))
for (s in stages) {
  message("\n######## ", s)
  t0 <- Sys.time()
  source(s, local = new.env())
  message(sprintf("%s done in %.1f min", s, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
if (interactive()) shiny::runApp("shiny")
