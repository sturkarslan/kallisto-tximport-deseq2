library(tidyverse)
##### create_conditions:  for comparison #####
create_conditions <- function(org){

  cat("Creating conditions for DEG comparison \n\n")
  ## Get list of abundance files from kallisto results
  cat("\t - Getting list of abundance files from kallisto \n\n")
  my.dir <- paste("../rnaseq_results_kallisto/",org,sep="")
  sample_id <- dir(file.path(my.dir))

  ## list of abundance files.
  files <- file.path(my.dir, sample_id, "abundance.h5")
  names(files) <- sample_id
  cat("\t\t - Collected abundance files for ", length(files), "files in ", org, "\n\n")


  ## separate meta data info
  meta1 <- as_tibble(sample_id) %>%
    separate(value, into = c("No","pH1","pH2","Rep","Time"), remove = F, sep = "_") %>%
    mutate(path = paste("../rnaseq_results_kallisto/",org,"/",value,"/abundance.h5",sep = "")) %>%
    mutate(condition = paste(pH1,pH2,Time, sep="_")) %>%
    rename(sample_name = value) %>%
    mutate(treatment = if_else(Time == "Tminus4", "control", "treatment"))

  ## Create conditions for DEG analysis
  samples <- unique(meta1$condition)
  conditions.1 <- expand.grid(samples, samples) %>%
    filter(Var1 != Var2) %>%
    filter(grepl("pH_7", Var2)) %>%
    separate(Var1, into = c("pH11","pH12","Time1"), remove = F, sep = "_") %>%
    separate(Var2, into = c("pH21","pH22","Time2"), remove = F, sep = "_") %>%
    filter(Time1 == Time2 & pH11 != pH22) %>%
    mutate(conditions = paste("condition_",Var1,"_vs_",Var2, sep = ""))

  conditions.2 <- expand.grid(samples, samples) %>%
    separate(Var1, into = c("pH11","pH12","Time1"), remove = F, sep = "_") %>%
    separate(Var2, into = c("pH21","pH22","Time2"), remove = F, sep = "_") %>%
    filter(pH12 == 6  & pH22 == 6 & Time2 == "Tminus4") %>%
    filter(Var1 != Var2) %>%
    mutate(conditions = paste("condition_",Var1,"_vs_",Var2, sep = ""))

  conditions.3 <- expand.grid(samples, samples) %>%
    separate(Var1, into = c("pH11","pH12","Time1"), remove = F, sep = "_") %>%
    separate(Var2, into = c("pH21","pH22","Time2"), remove = F, sep = "_") %>%
    filter(pH12 == 7  & pH22 == 7 & Time2 == "Tminus4") %>%
    filter(Var1 != Var2) %>%
    mutate(conditions = paste("condition_",Var1,"_vs_",Var2, sep = ""))

  conditions.final <- bind_rows(conditions.1, conditions.2, conditions.3)

  ### Uncomment to run specific conditions
  conditions <- conditions.final %>%
    #filter(timepoint1 == timepoint2) %>%
    #filter(Var1 != Var2) %>%
    pull(conditions) %>%
    unique()

  cat("\t Returned ", length(conditions), " conditions \n")
  print(conditions)
  return(conditions)
}
