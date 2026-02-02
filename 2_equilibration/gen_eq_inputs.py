#!/usr/bin/env python3
"""
Generate equilibration input files for each replica from a template.
Usage: python setup_equilibrate.py -i <template_file> -t <temperatures_file> -p <topology> -r <restart>
"""

import argparse
import random

def setup_equilibrate(template_file, temp_file, topology, restart):
    # Read temperatures
    with open(temp_file, 'r') as f:
        temperatures = [float(line.strip()) for line in f if line.strip()]
    
    nrep = len(temperatures)
    print(f"Setting up {nrep} replicas...")
    
    # Read template
    with open(template_file, 'r') as f:
        template = f.read()
    
    # Create groupfile
    with open('equilibrate.groupfile', 'w') as groupfile:
        for count, temp in enumerate(temperatures, start=1):
            rep = f"{count:03d}"
            print(f"TEMPERATURE: {temp} K ==> FILE: equilibrate.mdin.{rep}")
            
            # Replace placeholders
            content = template.replace('XXXXX', str(temp))
            content = content.replace('RANDOM_NUMBER', str(random.randint(1, 999999)))
            
            # Write input file
            with open(f'equilibrate.mdin.{rep}', 'w') as f:
                f.write(content)
            
            # Add to groupfile
            groupfile.write(
                f"-O -rem 0 -i equilibrate.mdin.{rep} "
                f"-o equilibrate.mdout.{rep} -c {restart} "
                f"-r equilibrate.rst.{rep} -x equilibrate.mdcrd.{rep} "
                f"-inf equilibrate.mdinfo.{rep} -p {topology}\n"
            )
        
        groupfile.write("#\n")
    
    print(f"N REPLICAS = {nrep}")
    print("Done.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate equilibration input files for REMD replicas')
    parser.add_argument('-i', '--input', required=True, help='Template input file (mdin)')
    parser.add_argument('-t', '--temperatures', required=True, help='Temperatures file')
    parser.add_argument('-p', '--topology', required=True, help='Topology file (prmtop)')
    parser.add_argument('-r', '--restart', required=True, help='Input restart/coordinate file')
    
    args = parser.parse_args()
    setup_equilibrate(args.input, args.temperatures, args.topology, args.restart)