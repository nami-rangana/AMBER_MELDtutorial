#!/usr/bin/env python3
"""
Generate REMD input files for each replica from a template.
Usage: python setup_remd.py -i <template_file> -t <temperatures_file> -p <topology> -r <restart_path>
"""

import argparse
import random
import os

def setup_remd(template_file, temp_file, topology, restart_path):
    with open(temp_file, 'r') as f:
        temperatures = [float(line.strip()) for line in f if line.strip()]
    
    nrep = len(temperatures)
    print(f"Setting up {nrep} replicas for REMD...")
    
    with open(template_file, 'r') as f:
        template = f.read()
    
    with open('remd.groupfile', 'w') as groupfile:
        for count, temp in enumerate(temperatures, start=1):
            rep = f"{count:03d}"
            print(f"TEMPERATURE: {temp} K ==> FILE: remd.mdin.{rep}")
            
            content = template.replace('XXXXX', str(temp))
            content = content.replace('RANDOM_NUMBER', str(random.randint(1, 999999)))
            
            with open(f'remd.mdin.{rep}', 'w') as f:
                f.write(content)
            
            restart_file = restart_path.replace('{rep}', rep)
            
            if not os.path.exists(restart_file):
                print(f"WARNING: Restart file not found: {restart_file}")
            
            groupfile.write(
                f"-O -rem 1 -remlog rem.log -i remd.mdin.{rep} "
                f"-o remd.mdout.{rep} -c {restart_file} "
                f"-r remd.rst.{rep} -x remd.mdcrd.{rep} "
                f"-inf remd.mdinfo.{rep} -p {topology}\n"
            )
        
        groupfile.write("#\n")
    
    print(f"N REPLICAS = {nrep}")
    print("Done.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate REMD input files for replica exchange')
    parser.add_argument('-i', '--input', required=True, help='Template input file (mdin)')
    parser.add_argument('-t', '--temperatures', required=True, help='Temperatures file')
    parser.add_argument('-p', '--topology', required=True, help='Topology file (prmtop)')
    parser.add_argument('-r', '--restart', required=True, 
                        help='Restart file path pattern (use {rep} for replica number, e.g., ../2_equilibrate/equilibrate.rst.{rep})')
    
    args = parser.parse_args()
    setup_remd(args.input, args.temperatures, args.topology, args.restart)