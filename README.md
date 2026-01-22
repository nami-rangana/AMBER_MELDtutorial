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

Now we can move on to minimization.<br>
Here, out main goal is to use MELD to guid our protein to fold into it's native structure by imposing some restrains we already know. Therefore we do a simple minimization.
* [min.in](1_system_setup/min.in)
{{< include file="1_system_setup/leap.in" language="bash" >}}
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
## 2. Equlibration
### Detrminig number of replicas
Before we dive into replica exchange simulations, we need to figure out two things: how many replicas we'll need and what temperature each one should be at.

The key is making sure neighboring replicas have overlapping potential energy distributions – this gives us good exchange probabilities. Using temperatures that are closer together increases this overlap, but it also means you'll need more replicas to cover your temperature range. Finding the sweet spot for temperature distribution can be tricky, and there are lots of different approaches out there.

To get started, you'll typically need to know how many atoms are in your system. This is easy to find – just look at your topology or coordinate files. In fact, the very first number in any coordinate file tells you exactly how many atoms you're working with!"
