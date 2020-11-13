#!/bin/bash

set -ev

#python3 ../plotting/parse-binning.py sad*.yaml && python3 ../plotting/parse-replicas.py r-pieces.yaml && python3 ../plotting/analyze-boundaries.py ./r-pieces ./sad-pieces-0.001 ./sad-pieces-0.01 ./sad-pieces-0.1

#Parse the Data Files
#python3 ../plotting/parse-binning.py sad*.yaml
python3 ../plotting/parse-replicas.py r-pieces.yaml
