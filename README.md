# AMBER_MELDtutorial
This repository has a tutorial on how to implement Replica Exchange MELD simulation in AMBER molecular dynamics software package.

## 1. System setup

In this section we will build the system via Leap and run minimization via Sander. Here a brief description of the system and the procedure used to generate the topology and coordinate files.

For this tutorial we will generate structure of the 3GB1 protein using the sequence FASTA file. Use following command to download the sequence.
```
curl -o 3gb1.fasta https://www.rcsb.org/fasta/entry/3GB1
```
Now lets bild the system to simulate using Leap. Here we will use ff19SB forcefield and mbondii2 radii that are appropriate for the igb=5 option in sander. We use [makeLeap.x](1_system_setup/makeLeap.x) to make the Leap input file ([leap.in](1_system_setup/leap.in)) using the [3gb1.fasta](1_system_setup/3gb1.fasta) file we downloaded before.

* [makeLeap.x](1_system_setup/makeLeap.x)

```
#!/bin/bash

sequence=$(grep -v "^>" 3gb1.fasta | tr -d '\n')
convert_aa() {
    case $1 in
        A) echo "ALA" ;;
        C) echo "CYS" ;;
        D) echo "ASP" ;;
        E) echo "GLU" ;;
        F) echo "PHE" ;;
        G) echo "GLY" ;;
        H) echo "HIS" ;;
        I) echo "ILE" ;;
        K) echo "LYS" ;;
        L) echo "LEU" ;;
        M) echo "MET" ;;
        N) echo "ASN" ;;
        P) echo "PRO" ;;
        Q) echo "GLN" ;;
        R) echo "ARG" ;;
        S) echo "SER" ;;
        T) echo "THR" ;;
        V) echo "VAL" ;;
        W) echo "TRP" ;;
        Y) echo "TYR" ;;
        *) echo "UNK" ;;
    esac
}

three_letter_seq=""
for (( i=0; i<${#sequence}; i++ )); do
    aa="${sequence:$i:1}"
    three_letter_seq="$three_letter_seq $(convert_aa $aa)"
done

# Create the leap.in file
cat > leap.in << EOF
source leaprc.protein.ff19SB
set default PBradii mbondi2
pro = sequence { ACE$three_letter_seq NHE }
saveamberparm pro 3gb1.prmtop 3gb1.inpcrd
quit
EOF

echo "leap.in file created successfully!"
```
Execute this with following commands:
```
chmod +x makeLeap.x
./makeLeap.x
```
This will create leap.in:
```
source leaprc.ff99SB
set default PBradii mbondi2
pro = sequence { ACE MET THR TYR LYS LEU ILE LEU ASN GLY LYS THR LEU LYS GLY GLU THR THR THR GLU ALA VAL ASP ALA ALA THR ALA GLU LYS VAL PHE LYS GLN TYR ALA ASN ASP ASN GLY VAL ASP GLY GLU TRP THR TYR ASP ASP ALA THR LYS THR PHE THR VAL THR GLU NHE }
saveamberparm pro 3gb1.prmtop 3gb1.inpcrd
quit
```
Now we use this inputfile to generate topology ([3gb1.prmtop](1_system_setup/3gb1.prmtop)) and coordinate file ([3gb1.inpcrd](1_system_setup/3gb1.inpcrd)) with: ```tleap -f leap.in```<br>
Check the leap output and log file for any errors.


## 2. Minimization
Now we can move on to minimization.<br>
Here, out main goal is to use MELD to guid our protein to fold into it's native structure by imposing some restrains we already know. Therefore we do a simple minimization.
* [min.in](1_system_setup/min.in)
```
energy minimization
 &cntrl
  imin=1, maxcyc=1000, ncyc=500,
  ntwr = 1000, ntpr = 100,
  cut = 999.0, rgbmax = 999.0,
  ntb = 0, igb = 5, saltcon = 0.0,
 /
 ```
 ```
 $AMBERHOME/bin/sander -O -i min.in -o min.out -p 3gb1.prmtop -c 3gb1.inpcrd -r min.rst 2> min.err
 ```
## 3. Other requirements
### Number of replicas & Temperature scaling
Before we dive into replica exchange simulations, we need to figure out two things: how many replicas we'll need and what temperature each one should be at.

The key is making sure neighboring replicas have overlapping potential energy distributions – this gives us good exchange probabilities. Using temperatures that are closer together increases this overlap, but it also means you'll need more replicas to cover your temperature range. Finding the sweet spot for temperature distribution can be tricky, and there are lots of different approaches out there.

To get started, you'll typically need to know how many atoms are in your system. This is easy to find – just look at your topology or coordinate files. In fact, the very first number in any coordinate file tells you exactly how many atoms you're working with! -- in our case 861 atoms"
```
head -n2 3gb1.inpcrd | tail -n1
```

Now that we know our system has 861 atoms, we can figure out how many replicas to use and what temperatures to assign them. As a general rule of thumb, the number of replicas is usually related to the square root of the number of atoms, and temperatures are typically set up in a geometric progression (where each temperature is a constant multiple of the previous one).

For this tutorial we choose to have 30 replicas between 270 K - 600 K scaled using geometric tempreature scaler in [temperature_scale.py](Helpers/temperature_scale.py).
```
#!/usr/bin/env python3
import sys
import math

def generate_constant(num_replicas, temperature):
    temps = []
    for i in range(num_replicas):
        alpha = i / (num_replicas - 1) if num_replicas > 1 else 0.0
        temps.append(temperature)
    return temps

def generate_linear(num_replicas, alpha_min, alpha_max, temp_min, temp_max):
    temps = []
    delta_alpha = alpha_max - alpha_min
    delta_temp = temp_max - temp_min
    
    for i in range(num_replicas):
        alpha = i / (num_replicas - 1) if num_replicas > 1 else 0.0
        
        if alpha <= alpha_min:
            temps.append(temp_min)
        elif alpha <= alpha_max:
            frac = (alpha - alpha_min) / delta_alpha
            temp = temp_min + frac * delta_temp
            temps.append(temp)
        else:
            temps.append(temp_max)
    
    return temps

def generate_geometric(num_replicas, alpha_min, alpha_max, temp_min, temp_max):
    temps = []
    delta_alpha = alpha_max - alpha_min
    
    for i in range(num_replicas):
        alpha = i / (num_replicas - 1) if num_replicas > 1 else 0.0
        
        if alpha <= alpha_min:
            temps.append(temp_min)
        elif alpha <= alpha_max:
            frac = (alpha - alpha_min) / delta_alpha
            delta = math.log(temp_max) - math.log(temp_min)
            temp = math.exp(delta * frac + math.log(temp_min))
            temps.append(temp)
        else:
            temps.append(temp_max)
    
    return temps

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: ./generate_temperatures.py <method> <num_replicas> <args...>")
        sys.exit(1)
    
    method = sys.argv[1]
    num_replicas = int(sys.argv[2])
    
    if method == "constant":
        temperature = float(sys.argv[3])
        temps = generate_constant(num_replicas, temperature)
    elif method == "linear":
        alpha_min, alpha_max = float(sys.argv[3]), float(sys.argv[4])
        temp_min, temp_max = float(sys.argv[5]), float(sys.argv[6])
        temps = generate_linear(num_replicas, alpha_min, alpha_max, temp_min, temp_max)
    elif method == "geometric":
        alpha_min, alpha_max = float(sys.argv[3]), float(sys.argv[4])
        temp_min, temp_max = float(sys.argv[5]), float(sys.argv[6])
        temps = generate_geometric(num_replicas, alpha_min, alpha_max, temp_min, temp_max)
    else:
        print(f"Unknown method: {method}")
        sys.exit(1)
    
    with open("temperatures.dat", "w") as f:
        for temp in temps:
            f.write(f"{temp:.1f}\n")
    
    print(f"temperatures.dat generated with {num_replicas} replicas using {method} method")
```
```
chmod +x temperature_scale.py
./temperature_scale.py geometric 30 0.0 1.0 270.0 600.0
```
Now you have scaled temperatures for each replica. This generates a file containing the temperature in Kelvin for each replica: [temperatures.dat](Helpers/temperatures.dat).

> ### Note: <br>
> * Constant temperature : ```./generate_temperatures.sh constant <num_replicas> <temperature>```<br>
> * Linear scaling : ```./generate_temperatures.sh linear <num_replicas> <alpha_min> <alpha_max>  <temp_min> <temp_max>```<br>
> * Geometric scaling : ```./generate_temperatures.sh geometric <num_replicas> <alpha_min> <alpha_max>  <temp_min> <temp_max>```

### Restraints 

Now we need to create restraints for our system. Restraints help guide the protein folding simulation by defining which atoms should be close together based on what we know about protein structure.

You could create these restraints manually, but that would be tedious! Instead, we'll use a Python script that automatically generates them for us. The script [gen_restraints_prot.py](Helpers/gen_restraints_prot.py) needs three input files: your topology file (prmtop), your sequence file (fasta), and a secondary structure file [ss.dat](Helpers/ss.dat).

The script creates two output files: an AMBER-compatible restraint file (restraints.rst) and a restraint index file [restraint_index.txt](Helpers/restraint_index.txt) that tracks the restraints -- collection & group .

These restraints are organized into two collections: strand-pair restraints that keep beta-sheet structures together, and hydrophobic restraints that ensure hydrophobic residues interact appropriately. This combination helps the protein fold into its correct structure during the simulation.


