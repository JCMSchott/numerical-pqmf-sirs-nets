#$ -S /bin/sh

programa=${1}

label2=${2}

echo $@ > argumentos.log

time ./${programa} ${label2}
