#!/bin/sh

set -ev

../plotting/parse-replicas.py --reparse z-*.cbor
../plotting/parse-replicas.py z-*/*.cbor

# ../plotting/parse-binning.py sad-*.cbor
# ../plotting/parse-binning.py sad-*/*.cbor

../plotting/parse-binning.py itwl-*.cbor
../plotting/parse-binning.py itwl-*/*.cbor
