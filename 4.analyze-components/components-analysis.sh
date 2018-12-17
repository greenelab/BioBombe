#!/bin/bash

# The following scripts will reproduce the figures and results for the
# compression component analysis module

Rscript --vanilla scripts/nbconverted/1.visualize-reconstruction.r
Rscript --vanilla scripts/nbconverted/2.visualize-sample-correlation.r
Rscript --vanilla scripts/nbconverted/3.sample-correlation-main-figures.r

