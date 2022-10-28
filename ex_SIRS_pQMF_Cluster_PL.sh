#!/bin/bash

ind_amostra=${1}
ind_lamb=${2}
# ind_noh=${3}
ind_noh=$((${ind_lamb} % 10))
num_base=1

tam_rede=${num_base}0000000

grau_min=3; gama_exp=3.5

lamb0=0.0

divisor=5.0

lambdaf=1.0

alp=0.1

ind_soCalculaIPR=0
#=======================================================================
#maquina='cluster'
#=======================================================================
maquina='cluster_ind'
maq=4
noh=(0-0 0-1 0-12 0-14 0-2 0-20 0-3 0-5 0-6 0-8)
no=${noh[${ind_noh}]}
#no=4-3
echo ${no}
#-----------------------------------------------------------------------
#compute-0-0             linux-x64      24  4.00   63.0G    3.1G   72.3G     0.0
#compute-0-1             linux-x64      24  4.00   63.0G    1.1G   72.3G     0.0
#compute-0-10            linux-x64      12  2.00   47.2G    1.0G   72.3G   11.6M
#compute-0-11            linux-x64      24  3.00   31.4G    1.0G   72.3G   17.0M
#compute-0-12            linux-x64      24  4.01   63.0G    8.3G   72.3G   15.7M
#compute-0-13            linux-x64      12     -   31.3G       -   72.3G       -
#compute-0-14            linux-x64      24  4.01   63.0G    1.1G   72.3G   14.1M
#compute-0-15            linux-x64      24     -   62.8G       -   72.3G       -
#compute-0-16            linux-x64      24  3.00   63.0G    5.2G   72.3G   14.2M
#compute-0-17            linux-x64      12  2.00   31.2G  988.5M   72.3G   10.5M
#compute-0-2             linux-x64      24  4.01   63.0G    1.1G   72.3G   12.0M
#compute-0-20            linux-x64      24  5.00   62.8G    1.2G   72.3G   13.7M
#compute-0-3             linux-x64      24  4.00   63.0G    2.8G   72.3G   13.4M
#compute-0-4             linux-x64      24     -   39.3G       -   72.3G       -
#compute-0-5             linux-x64      24  3.56   63.0G    4.9G   72.3G     0.0
#compute-0-6             linux-x64      24  2.98   63.0G    1.1G   72.3G   14.5M
#compute-0-7             linux-x64      24  3.00   47.2G  835.6M   72.3G   13.0M
#compute-0-8             linux-x64      24  4.00   63.0G    1.1G   72.3G   14.2M
#compute-0-9             linux-x64      24  3.07   47.2G    5.2G   72.3G     0.0
#compute-1-18            linux-x64      24     -   23.5G       -   72.3G       -
#compute-1-19            linux-x64      24     -   31.5G       -   72.3G       -
#compute-1-21            linux-x64      24     -   23.5G       -   72.3G       -
#=======================================================================
#maquina='casa'
#=======================================================================
#flags='-check all -traceback'

#flags='-heap-arrays -O3 -g -fp-stack-check'

#flags='-traceback -check all -heap-arrays -O3 -ipo -g -fp-stack-check'

flags='-Ofast -ipo'
#=======================================================================
nzeros=$( echo "scale=0; l(${tam_rede}/${num_base})/l(10) " | bc -l)

IFS='.'

arr=( ${gama_exp} )

gama_sp=${arr[0]}${arr[1]}

alp_arr=( ${alp} )

IFS=' '
##########################################################
dependencias='mod_rndgen_multiple.f90 geraRede.f90 mod_tools_redes.f90 mod_types.f90'
principal='main_SIRS_PQMF_simplificado.f90'

executavel='pQMF_n'${num_base}${nzeros}g${gama_sp}'a'${alp_arr[0]}'p'${alp_arr[1]}'_ams_'${ind_amostra}'_indLambd_'${ind_lamb}


rm ${executavel} > erroRmExe.log

ifort ${dependencias} ${principal} ${flags} -o ${executavel}

rm -r *.mod

#=======================================================================
if [[ ${maquina} == 'casa' ]];then
   #--------------------------------------------------------------------
   time ./${executavel} ${ind_amostra} ${tam_rede} ${grau_min} ${gama_exp} ${lamb0} ${divisor} ${lambdaf} ${alp} ${ind_lamb} ${ind_soCalculaIPR}
   #--------------------------------------------------------------------
elif [[ ${maquina} == 'cluster' ]];then
   #--------------------------------------------------------------------
   qsub -N ${executavel} -cwd executa_sirs_pqmf_cluster_PL.sh ${executavel} ${ind_amostra} ${tam_rede} ${grau_min} ${gama_exp} ${lamb0} ${divisor} ${lambdaf} ${alp} ${ind_lamb} ${ind_soCalculaIPR}
   #--------------------------------------------------------------------
elif [[ ${maquina} == 'cluster_ind' ]];then
   qsub -N ${executavel} -q all.q@compute-${no} -cwd executa_sirs_pqmf_cluster_PL.sh ${executavel} ${ind_amostra} ${tam_rede} ${grau_min} ${gama_exp} ${lamb0} ${divisor} ${lambdaf} ${alp} ${ind_lamb} ${ind_soCalculaIPR}
   echo "Escolhemos o noh "${no}
else
   #--------------------------------------------------------------------
   echo "Maquina nao-reconhecida. Abortando script..."
   #--------------------------------------------------------------------
   exit
fi
