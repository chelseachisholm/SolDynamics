#############################
## SOLDYNAMICS DRAKE PLAN ###
#############################

rm(list=ls())
# Load library
library("drake")
library("tidyverse")
library("rjags")
library(geosphere)
library(rgdal)

# drake configurations
pkgconfig::set_config("drake::strings_in_dots" = "literals")

# trick
pn <- . %>% print(n = Inf)

# source scripts
source("./data/readdata_jags_dem.R")
#source("./code/JAGS/occupancymodel_jags_07052019.R")

# Import Data
ImportDrakePlan <- drake_plan(
  jags_data = Importsoldat()
)

# ModelDrakePlan <- drake_plan(
#   jags_model = Run_JAGSmodel(jags_data)
# )

#AnalyzeDrakePlan <- drake_plan(
#  JAGS_output = Run_JAGSmodel()
#)

MyPlan <- ImportDrakePlan

conf <- drake_config(MyPlan)
make(MyPlan)
loadd()
vis_drake_graph(conf, targets_only = TRUE)
