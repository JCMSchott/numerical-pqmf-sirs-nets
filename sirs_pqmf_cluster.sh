#!/bin/bash

ind_amostra=${1}

n_base=1

tam_rede=${n_base}0000

grau_RRN=6; grau_Star=998

lamb0=0.0

divisor=5.0

lambdaf=1.0

alp=0.1
ind_lamb=${2}   
ind_soCalculaIPR=0
#=======================================================================
#flags='-check all -traceback'

#flags='-heap-arrays -O3 -g -fp-stack-check'

flags='-traceback -check all -heap-arrays -O3 -g -fp-stack-check'

#flags=''
#=======================================================================
nzeros=$( echo "scale=0; l(${tam_rede}/${n_base})/l(10) " | bc -l)

IFS='.'

arr=( ${gama_exp} )

gama_sp=${arr[0]}${arr[1]}

alp_arr=( ${alp} )

IFS=' '
##########################################################



dependencias='mod_rndgen_multiple.f90 geraRede.f90 mod_tools_redes.f90 mod_types.f90'
principal='main_SIRS_PQMF_RRN_Plus_Star_simplificado.f90'

executavel='pqmf'${n_base}${nzeros}g${gama_sp}'a'${alp_arr[0]}'p'${alp_arr[1]}'ams'${ind_amostra}'indLambd'${ind_lamb}

rm ${executavel} > erroRmExe.log

ifort ${dependencias} ${principal} ${flags} -o ${executavel}

#ifx ${dependencias} ${principal} ${flags} -o ${executavel}

rm -r *.mod

time ./${executavel} ${ind_amostra} ${tam_rede} ${grau_RRN} ${grau_Star} ${lamb0} ${divisor} ${lambdaf} ${alp} ${ind_lamb} ${ind_soCalculaIPR}

#qsub -N ${executavel} -cwd executa_pqmf.sh ${executavel} ${ind_amostra} ${tam_rede} ${grau_RRN} ${grau_Star} ${lamb0} ${divisor} ${lambdaf} ${alp} ${ind_lamb} ${ind_soCalculaIPR}
