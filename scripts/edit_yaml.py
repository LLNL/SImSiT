import yaml
from configparser import ConfigParser
import sys

path = sys.argv[1]
sim_name = sys.argv[2]
sim_type = sys.argv[3]
n_obs = sys.argv[4]
seed = sys.argv[5]
n_sat = sys.argv[6]
ndx = sys.argv[7]
ndy = sys.argv[8]

print("Creating config files")
config = yaml.safe_load(open(path+'/config/'+sim_type+'.yaml'))
config['outdir'] = 'branches/'+str(sim_type)+'/'+str(sim_name)
config['n_obs'] = int(n_obs)
config['n_sat'] = int(n_sat)
config['seed'] = int(seed)
config['instrument']['image_shape'] = [int(ndx), int(ndy)]

with open('config/'+sim_type+'_'+sim_name+'.yaml', "w+") as yaml_f:
    yaml.dump(config, yaml_f, indent=4, default_flow_style=False)





