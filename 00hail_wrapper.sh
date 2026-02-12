#!/bin/bash

#command
usage() {
    echo "Usage: $0 --phecode <PHECODE> --pop <POP>"
    exit 1
}

#command-line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --phecode)
            PHECODE=$2
            shift 2
            ;;
        --pop)
            POP=$2
            shift 2
            ;;
        *)
            echo "unknown flag: $1"
            usage
            ;;
    esac
done

#check for required arguments
if [[ -z "$PHECODE" || -z "$POP" ]]; then
    usage
fi

#github repo path
REPO=$HOME/mesa_coloc

#download hail table
#python "$REPO/01pull_data.py" --phecode "$PHECODE" --pop "$POP"

#format hail tables
Rscript "$REPO/02table_format.R" --phecode "$PHECODE" --pop "$POP"
