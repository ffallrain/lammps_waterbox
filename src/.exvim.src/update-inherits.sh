#!/bin/bash
export DEST="./.exvim.src"
export TOOLS="/home/zengping/.vim/tools/"
export TMP="${DEST}/_inherits"
export TARGET="${DEST}/inherits"
sh ${TOOLS}/shell/bash/update-inherits.sh
