#!/bin/bash

label1=${1}
label2=${2}

dependencias='mod_rndgen_multiple.f90 geraRede.f90 mod_tools_redes.f90 mod_types.f90' # mod_pqmf.f90'
principal='sirs_pqmf_rede_real_est.f90'

executavel='ex_pqmfest_rede_real'

#flags='-check all -traceback'

#flags='-heap-arrays -O3 -g -fp-stack-check'

#flags='-traceback -check all -heap-arrays -O3 -g -fp-stack-check'

flags=''

ifort ${dependencias} ${principal} ${flags} -o ${executavel}

rm -r *.mod

time ./${executavel} ${label1} ${label2}
