#!/bin/bash
export DEST="./.exvim.spblas"
export TOOLS="/home/zzh/.vim/tools/"
export EXCLUDE_FOLDERS=""
export TMP="${DEST}/_ID"
export TARGET="${DEST}/ID"
sh ${TOOLS}/shell/bash/update-idutils.sh
