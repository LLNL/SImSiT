#!/bin/sh

path=$1
env_path=$2
config_file=$3
sim_type=$4
sim_name=$5

cd ${path}
source ${env_path}/bin/activate
python3 ${path}/scripts/simulate.py ${config_file}
wait

mv ${path}/${sim_type}_${sim_name}.out ${path}/branches/${sim_type}/${sim_name}
mv ${path}/config/${sim_type}_${sim_name}.yaml ${path}/branches/${sim_type}/${sim_name}