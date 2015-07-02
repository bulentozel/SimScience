rm(list = ls())
source("LibrarySimScience.R")
ALFA = 0.5
WORK_DIR = getwd()
PLOTS_DIR = "../plots/"
DATA_DIR = "../data/"


AxA.data <- read.table("../data/AxA.csv", header = T, sep = ";")
Community <- initialize.community(AxA.data)
KxK.data <- read.table("../data/KxK.csv", header = T, sep = ";")
Community <- add.knowledge(KxK.data, Community)


run.simulation(Community, T.it = 10)



KxK.Edgelist <- form.network(Community, type = "K")
plot.as.network(KxK.Edgelist, type = "K", f.name = "kxk.0")

AxA.Edgelist <- form.network(Community, type = "A")
plot.as.network(AxA.Edgelist, type = "A", f.name = "axa.0")



