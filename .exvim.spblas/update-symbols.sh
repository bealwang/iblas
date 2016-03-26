#!/bin/bash
export DEST="./.exvim.spblas"
export TOOLS="/home/zzh/.vim/tools/"
export TMP="${DEST}/_symbols"
export TARGET="${DEST}/symbols"
sh ${TOOLS}/shell/bash/update-symbols.sh
