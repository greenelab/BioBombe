#!/bin/bash

# Gregory Way 2019
# Build top feature matrices for downstream analysis comparing prediction features

# Top 200 features
python scripts/compile-top-matrices.py

# Top 1 feature
python scripts/compile-top-matrices.py --num_top_features 1

# Top 200 features, but randomize their selection
python scripts/compile-top-matrices.py --random_features
