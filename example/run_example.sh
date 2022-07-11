#!/bin/bash

env MOL_MASK=":1" REFTOP="$PWD/top_traj_files/drug_25.top" \
    TRAJDIR="$PWD/top_traj_files" ALIAS_MDCRD="drug_25"  \
    DO_CC_MLA=1  CUTOFF=" -1 "  \
    SCRATCH="/scratch" ../calc_sconform.sh

