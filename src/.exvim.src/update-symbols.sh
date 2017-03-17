#!/bin/bash
export DEST="./.exvim.src"
export TOOLS="/home/zengping/.vim/tools/"
export TMP="${DEST}/_symbols"
export TARGET="${DEST}/symbols"
sh ${TOOLS}/shell/bash/update-symbols.sh
