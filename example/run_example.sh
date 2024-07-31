#!/bin/bash

env NPROCS=6 MOL_MASK=":1" REFTOP="$PWD/top_traj_files/drug_25.top" \
    TRAJDIR="$PWD/top_traj_files" PREFIX_MDCRD="drug_25" SUFFIX_MDCRD=".mdcrd" \
    DO_CC_MLA=1  CUTOFF=" -1 "  \
    SCRATCH="/scratch" ../calc_sconform.sh

