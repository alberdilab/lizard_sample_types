#!/usr/bin/env Rscript

####
# Render final report
# By Antton Alberdi (antton.alberdi@sund.ku.dk)
# 03/03/2024
# Dependencies: pandoc
####

#Load libraries
suppressPackageStartupMessages(library(bookdown))
suppressPackageStartupMessages(library(htmlwidgets))
suppressPackageStartupMessages(library(webshot))

#Install install_phantomjs
#webshot::install_phantomjs()

#Render it as github pages ()
render_book(input = ".", output_format = "bookdown::gitbook", output_dir = "docs")

