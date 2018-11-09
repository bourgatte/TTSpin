#!/usr/bin/env python

import os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output-file",help="output file prefix. [Default: %(default)s] ", action="store", default = 'TEST')
parser.add_argument("-c", "--config-parameters",help="Configuration file:  pythia_H.conf or pythia_Z.conf   [Default: %(default)s] ",  type=str, action="store", default = 'pythia_H.conf')
parser.add_argument("-t", "--collision-type",help="Collision type; 0 - ee, 1 - pp   [Default: %(default)s] ",   action="store", default = '1')
parser.add_argument("-n", "--events-number",help="Number of events   [Default: %(default)s] ",  action="store", default = '1000')
parser.add_argument("-ma","--mixing-angle",help="CP Higgs mixing angle (valid for Higgs config)   [Default: %(default)s] ",  action="store", default = '0.7853')
args = parser.parse_args()
cmd = './spin_correlation_pythia_tauola.exe ' + args.config_parameters     + ' ' + args.collision_type + ' ' +  args.events_number     + ' ' +  args.mixing_angle + ' ' + args.output_file
print cmd
os.system(cmd) 
