#!/bin/bash

label=${1}

dependencias='mod_rndgen_multiple.f90 geraRede.f90 mod_tools_redes.f90 mod_types.f90' # mod_pqmf.f90'
principal='main_SIRS_PQMF_simplificado.f90'

executavel='pqmfg23_10k_1'

#flags='-check all -traceback'

#flags='-heap-arrays -O3 -g -fp-stack-check'

#flags='-traceback -check all -heap-arrays -O3 -g -fp-stack-check -check uninit -check bounds -ftrapuv'

flags=''

ifort ${dependencias} ${principal} ${flags} -o ${executavel}

rm -r *.mod

time ./${executavel} ${label}
