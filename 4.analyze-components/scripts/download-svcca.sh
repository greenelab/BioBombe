#!/bin/bash
#
# Download the a script to perform Singular Vector Canonical Correlation Analysis
# Script retrieved from https://github.com/google/svcca
#
# Citation:
# Maithra Raghu, Justin Gilmer, Jason Yosinski, Jascha Sohl-Dickstein (2017).
# "SVCCA: Singular Vector Canonical Correlation Analysis for Deep Learning Dynamics and Interpretability".
# Neural Information Processing Systems (NIPS) 2017.

# The current version is slightly modified for python 3 compatibility
wget --output-document='scripts/cca_core.py' https://github.com/gwaygenomics/svcca/raw/82d26dffd6b65784206ace17bfcf9e481fd3c65c/cca_core.py
