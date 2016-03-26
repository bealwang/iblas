#!/bin/bash
export DEST="./.exvim.spblas"
export TOOLS="/home/zzh/.vim/tools/"
export FOLDERS="./include,./src"
export FILE_SUFFIXS=".*"
export TMP="${DEST}/_files"
export TARGET="${DEST}/files"
sh ${TOOLS}/shell/bash/update-filelist.sh
