#$ -S /bin/sh

programa=${1}
ind_amostra=${2}
tam_rede=${3}
grau_RRN=${4}
grau_Star=${5}
lamb0=${6}
divisor=${7}
lambdaf=${8}
alp=${9}
ind_lamb=${10}

ind_soCalculaIPR=${11}

echo $@ > argumentos.log

time ./${programa} ${ind_amostra} ${tam_rede} ${grau_RRN} ${grau_Star} ${lamb0} ${divisor} ${lambdaf} ${alp} ${ind_lamb} ${ind_soCalculaIPR}
