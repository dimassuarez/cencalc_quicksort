#!/bin/bash

# make sure that calc_sconform.sh adn get_tor.py have execution permissions

env NPROCS=6 MOL_MASK=":1" REFTOP="$PWD/top_traj_files/drug_25.top" \
    TRAJDIR="$PWD/top_traj_files" MDCRD_PREFIX="drug_25" MDCRD_SUFFIX=".mdcrd" \
    DO_CC_MLA=1  CUTOFF=" -1 "  \
    SCRATCH="/scratch" ../calc_sconform.sh

