#!/bin/bash

set -ev

#Parse the Data Files
python3 ../plotting/parse-replicas.py r-*.yaml r-*/*.cbor

python3 ../plotting/parse-binning.py sad-*1.yaml sad-*1/*.cbor *wl*1.yaml *wl*1/*.cbor