# Mon Mar 07 14:32:23 2022 ------------------------------
# Loading packages
library(RRphylo)
library(ape)
library(letsR)
library(visdat)
library(tidyverse)

# Extract data
dir()
anura_traits <- read.csv2("traits_amphibia.csv")
head(anura_traits)

vis_miss(anura_traits)

# Imputação dados faltantes - SVL/Body mass - Número de prole
# Extrair váriaveis climáticas - LetsR
# Taxa de evolução dos atributos
# modelos PGLS