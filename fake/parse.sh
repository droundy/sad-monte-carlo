#!/bin/bash

set -ev

python3 ../plotting/parse-replicas.py r-*.yaml r-*/*.cbor

python3 ../plotting/parse-binning.py sad-*.yaml sad-*/*.cbor *wl*.{yaml,cbor} *wl*/*.cbor