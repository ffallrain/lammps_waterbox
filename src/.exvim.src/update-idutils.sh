#!/bin/bash
export DEST="./.exvim.src"
export TOOLS="/home/zengping/.vim/tools/"
export TMP="${DEST}/_ID"
export TARGET="${DEST}/ID"
sh ${TOOLS}/shell/bash/update-idutils.sh
