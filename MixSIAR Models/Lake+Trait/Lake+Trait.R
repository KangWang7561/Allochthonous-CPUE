
rm(list = ls())
set.seed(756199513)

library(tidyr)
library(plyr) 
library(dplyr)
library(nlstools)
library(reshape2) 
library(ggplot2)
library(ggnewscale)
library(ggrepel)
library(cowplot)
library(grid)
library(R2WinBUGS)
library(splancs)
library(gridExtra)
library(rjags)
library(MixSIAR)
library(readxl)
library(writexl)

# Trait + Lake model ####
mix = load_mix_data(filename = "Consumer_data.csv",
                    iso_names = "d2H",
                    factors = c("Lake", "Trait"),
                    fac_random = c(TRUE, TRUE),
                    fac_nested = c(FALSE, FALSE),
                    cont_effects = "PC1")

source = load_source_data(filename = "Source_data.csv", 
                          source_factors = "Lake",
                          conc_dep = FALSE,
                          data_type = "means", mix)

discr = load_discr_data(filename = "TEF.csv", mix)

output_options = list(summary_save = TRUE,                 
                      summary_name = "summary_Lake+Trait", 
                      sup_post = TRUE,        
                      plot_post_save_pdf = FALSE,           
                      plot_post_name = "posterior_density",
                      sup_pairs = TRUE,   
                      plot_pairs_save_pdf = TRUE,    
                      plot_pairs_name = "pairs_plot",
                      sup_xy = TRUE,       
                      plot_xy_save_pdf = FALSE,
                      plot_xy_name = "xy_plot",
                      gelman = TRUE,
                      heidel = FALSE,  
                      geweke = TRUE,   
                      diag_save = TRUE,
                      diag_name = "diagnostics_Lake+Trait",
                      indiv_effect = FALSE,       
                      plot_post_save_png = FALSE, 
                      plot_pairs_save_png = FALSE,
                      plot_xy_save_png = FALSE,
                      diag_save_ggmcmc = FALSE,
                      return_obj = TRUE)

model_filename = "MixSIAR_model.txt"
resid_err = TRUE
process_err = TRUE
write_JAGS_model(model_filename, process_err, resid_err, mix, source)
jags.mod = run_model(run = "very long", mix, source, discr, model_filename,
                     alpha.prior=1, process_err = TRUE, resid_err = TRUE)

options(max.print = 99999)
saveRDS(jags.mod, file = paste0("MixSIAR_Lake+Trait", ".rds"))
saveRDS(mix, file = paste0("mix_Lake+Trait",".rds"))
saveRDS(source, file = paste0("source_Lake+Trait",".rds"))
saveRDS(discr, file = paste0("discr_Lake+Trait",".rds"))

# Posterior ####
mix = readRDS(file = paste0("mix_Lake+Trait",".rds"))
source = readRDS(file = paste0("source_Lake+Trait",".rds"))
jags = readRDS(file = paste0("MixSIAR_Lake+Trait",".rds"))

Posterior = jags[["BUGSoutput"]][["sims.list"]][["p.ind"]]
Posterior = reshape2::melt(Posterior)
colnames(Posterior) = c("Draw", "Individual", "Source", "Proportion")
Posterior$Source2 = factor(Posterior$Source, labels = source[["source_names"]])
Posterior$Source2 = factor(Posterior$Source2, levels = source[["source_names"]])

number.sim = jags[["BUGSoutput"]][["n.sims"]]
number.sources = length(source[["source_names"]])
gg = cbind(rep(rep(mix[["data"]][["Lake"]], each = number.sim), times = number.sources),  
           rep(rep(mix[["data"]][["Trait"]], each = number.sim), times = number.sources))
colnames(gg) = c("Lake", "Trait")
Posterior = cbind(Posterior, gg)

Posterior$Lake = factor(Posterior$Lake, levels = unique(mix[["data"]][["Lake"]]))
Posterior$Species = factor(Posterior$Species, levels = unique(mix[["data"]][["Species"]])[order(unique(mix[["data"]][["Species"]]))])
Posterior$Source2 = factor(Posterior$Source2, levels = source[["source_names"]])

Posterior = subset(Posterior, select = -c(Individual, Source))
colnames(Posterior) = c("Draw", "Proportion", "Source", "Lake", "Trait")

Posterior.Lake = Posterior %>%
  group_by(Lake, Source) %>% 
  summarise(y0.025 = quantile(Proportion, 0.025, na.rm = TRUE),
            y0.1= quantile(Proportion, 0.1, na.rm = TRUE),
            y0.25= quantile(Proportion, 0.25, na.rm = TRUE),
            y0.375= quantile(Proportion, 0.375, na.rm = TRUE),
            y0.4= quantile(Proportion, 0.4, na.rm = TRUE),
            y0.5= quantile(Proportion, 0.5, na.rm = TRUE),
            y0.625= quantile(Proportion, 0.625, na.rm = TRUE),
            y0.75= quantile(Proportion, 0.75, na.rm = TRUE),
            y0.8= quantile(Proportion, 0.8, na.rm = TRUE),
            y0.9= quantile(Proportion, 0.9, na.rm = TRUE),
            y0.975= quantile(Proportion,0.975, na.rm = TRUE))

Posterior.Trait = Posterior %>%
  group_by(Lake, Trait, Source) %>%  
  summarise(y0.025 = quantile(Proportion, 0.025, na.rm = TRUE),
            y0.1= quantile(Proportion, 0.1, na.rm = TRUE),
            y0.25= quantile(Proportion, 0.25, na.rm = TRUE),
            y0.375= quantile(Proportion, 0.375, na.rm = TRUE),
            y0.4= quantile(Proportion, 0.4, na.rm = TRUE),
            y0.5= quantile(Proportion, 0.5, na.rm = TRUE),
            y0.625= quantile(Proportion, 0.625, na.rm = TRUE),
            y0.75= quantile(Proportion, 0.75, na.rm = TRUE),
            y0.8= quantile(Proportion, 0.8, na.rm = TRUE),
            y0.9= quantile(Proportion, 0.9, na.rm = TRUE),
            y0.975= quantile(Proportion,0.975, na.rm = TRUE))

# write_xlsx(Posterior.Lake, path = "Posterior.Lake.xlsx")
# write_xlsx(Posterior.Trait, path = "Posterior.Trait.xlsx")

