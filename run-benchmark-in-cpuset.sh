#!/bin/bash

set -ev

DIR=`pwd`
echo Working in $DIR

if test $EUID -ne 0; then
   echo "This script must be run as root"
   exit 1
fi

USERSET='1'
for i in `seq 2 15`; do
    if test -e "/sys/devices/system/cpu/cpu$i/topology/thread_siblings_list" && ! grep 0 "/sys/devices/system/cpu/cpu$i/topology/thread_siblings_list"; then
        USERSET="$USERSET,$i"
    fi
done

cset shield --cpu=$USERSET
cset shield -e -- su - droundy -c "cd $DIR && cargo bench" || true
cset shield --reset
