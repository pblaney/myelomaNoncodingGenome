#!/usr/bin/env Rscript

# This script runs the tiled genome FishHook analysis

#########################
#####   Libraries   #####

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressPackageStartupMessages(library(gUtils))
suppressPackageStartupMessages(library(fishHook))

Sys.setenv(DEFAULT_BSGENOME = 'BSgenome.Hsapiens.UCSC.hg38::Hsapiens')
options(scipen = 999)

#########################
#####   Execution   #####

# Accept command line arguments as input
input_args <- commandArgs(trailingOnly = T)

muts <- input_args[1]

trophy_fish <- input_args[2]

covariates <- "covars_bait.rds"
hypothesis_tiles <- "genome_tiles.rds"
eligible_background <- "eligible_region.rds"

# From https://gitlab.com/sanjanalab/circle/-/blob/main/CIRCLE_notebook.ipynb#analysis7
# Creating an empty FishHook object
message("Loading the FishHook bait (covariates): ", covariates, " ...")
fish_cov <- readRDS(covariates)
message("Loading the hypothesis: ", hypothesis_tiles, " ...")
fish_hypo <- readRDS(hypothesis_tiles)
message("Loading the eligible region: ", eligible_background, " ...")
fish_eligible <- readRDS(eligible_background)
message("Loading the mutations: ", muts, " ...")
fish_muts <- readRDS(muts)

# Run basic model without covariates on each tile with the SNVs
message("Here fishy fishy ...")
fishing_pole <- fishHook::Fish(hypotheses = fish_hypo,
                               events = fish_muts,
                               eligible = fish_eligible,
                               idcol = "Matched_Norm_Sample_Barcode",
                               use_local_mut_density = FALSE,
                               nb = TRUE)

# Merging in the covariates
fishing_pole$merge(fish_cov)
message("Reeling the fish in ...")

saveRDS(object = fishing_pole, file = trophy_fish)
message("We caught a big one ...")
