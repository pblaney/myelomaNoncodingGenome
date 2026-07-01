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

fish <- input_args[1]

outfile_base <- input_args[2]

message("Let's make some sushi ...")
# Pull the catch from the freezer
tuna <- readRDS(file = fish)
message("Tuna ready to be made into sushi: ", fish," ...")
# Peak at full model before removing low eligibility tiles
message("Base fishHook model (Tuna) metrics ...")
tuna
# Remove nonsense and low eligibility tiles
sashimi <- tuna[!is.na(tuna$data$frac.eligible),]
sashimi <- sashimi[sashimi$data$frac.eligible > 0.75,]
message("High coverage tile fishHook model (Sashimi) metrics ...")
sashimi
# Compute p-values
message("Making sushi rolls ...")
sashimi$score()
# Save tiles in data.table
sushi_results <- gr2dt(sashimi$res)[order(p), ][, p:= as.character(p)]
saveRDS(object = sushi_results, file = paste0(outfile_base, "_results.sushi.rds"))
saveRDS(object = sashimi, file = paste0(outfile_base, ".sushi.rds"))
message("Sushi was eaten and enjoyed.....")

message("Let's make a QQ plot ...")
pdf(file = paste0(outfile_base, "_qq.sushi.pdf"))
sashimi$qqp(plotly = F)
dev.off()

message("Finally, a peak at the model ...")
saveRDS(object = summary(sashimi$model), file = paste0(outfile_base, "_model.sushi.rds"))
summary(sashimi$model)
