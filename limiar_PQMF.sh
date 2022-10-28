#!/bin/bash

##########################################################

amostra=${1}

tam=100000; grau_min=3; gama=3.5;
##########################################################
# dlambda = 0.0125/divisor

##########################################################
alp=0.5;
##########################################################

#flags='-check all -traceback'

#flags='-heap-arrays -O3 -g -fp-stack-check'

#flags='-traceback -check all -heap-arrays -O3 -g -fp-stack-check'

flags=''

##########################################################
# Nao mexer daqui pra baixo!
nzeros=$( echo "scale=0; l(${tam})/l(10) " | bc -l)

IFS='.'

arr=( ${gama} )

gama_sp=${arr[0]}${arr[1]}

alp_arr=( ${alp} )

IFS=' '
##########################################################
dependencias='mod_rndgen_multiple.f90 geraRede.f90 mod_tools_redes.f90'
principal='limiar_SIRS_PQMF.f90'
executavel='n'1${nzeros}g${gama_sp}'a'${alp_arr[0]}p${alp_arr[1]}'_'${amostra}
#############################################################################
rm ${executavel} &> erroRmExe.log

ifort ${dependencias} ${principal} ${flags} -o ${executavel}

rm -r *.mod

time ./${executavel} ${amostra} ${tam} ${grau_min} ${gama} ${alp}
