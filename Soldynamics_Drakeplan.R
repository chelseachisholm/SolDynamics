#############################
## SOLDYNAMICS DRAKE PLAN ###
#############################

rm(list=ls())
# Load library
library("drake")
library("tidyverse")
library("readxl")
library("lubridate")
library("e1071")
library("rjags")

# drake configurations
pkgconfig::set_config("drake::strings_in_dots" = "literals")

# trick
pn <- . %>% print(n = Inf)

# source scripts
source("./data/readdata_jags_allyears.R")
source("./code/JAGS/occupancymodel_jags_07052019.R")

# Import Data
ImportDrakePlan <- drake_plan(
  jags_data = Importandclean_soldat()
)

ModelDrakePlan <- drake_plan(
  jags_model = Run_JAGSmodel(jags_data)
)

#AnalyzeDrakePlan <- drake_plan(
#  JAGS_output = Run_JAGSmodel()
#)

MyPlan <- bind_rows(ImportDrakePlan, ModelDrakePlan)

conf <- drake_config(MyPlan)
make(MyPlan)
loadd()
vis_drake_graph(conf, targets_only = TRUE)
