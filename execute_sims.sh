#!/bin/sh

# --------------------------------------
# Set all parameters for pipeline here:
# ${sim_name} MUST match final 'outdir' directory listed in yaml file for specific run
sim_name=test
sim_type=target
debug=True
delete=True
ndx=1800 # 225 # 300 # 1800 # 75
ndy=1536 # 192 # 256 # 1536 # 64
n_obs=2
n_sat=1
seed=$RANDOM
time="24:00:00" # "00:05:00"
account=sdadp
path=/p/lustre2/pruett6/xfiles/xfiles_sim
env_path=/p/lustre2/pruett6/xfiles/venv
# --------------------------------------

config_file=${path}/config/${sim_type}_${sim_name}.yaml

if [[ -e ${path}/branches/${sim_type}/${sim_name} ]] 
then 
    if [[ $delete == "True" ]]
    then
        echo "Deleting ${path}/branches/${sim_type}/${sim_name}"
        rm -rf ${path}/branches/${sim_type}/${sim_name}
    elif [[ $delete != "True" ]]
    then
        echo "Directory ${path}/branches/${sim_type}/${sim_name} already exists. Please remove the directory or change 'sim_name'."
        exit
    fi
fi

cd ${path}
source ${env_path}/bin/activate
mkdir ${path}/branches/${sim_type}/${sim_name}
python3 setup.py install
wait

python3 scripts/edit_yaml.py $path $sim_name $sim_type $n_obs $seed $n_sat $ndx $ndy

if [[ $debug == "True" ]]
then
    sbatch --time="00:30:00" --partition=pdebug --output=${sim_type}_${sim_name}.out run_sim.sh $path $env_path $config_file $sim_type $sim_name
elif [[ $debug == "False" ]]
then 
    sbatch --time=${time} --account=${account} --output=${sim_type}_${sim_name}.out run_sim.sh $path $env_path $config_file $sim_type $sim_name
fi
