#!/usr/bin/env python3
"""
Generate AMBER restraint files for MELD protocol.
Creates distance restraints for secondary structure elements and hydrophobic interactions.
"""

import argparse
import numpy as np

HP_ATOMS = {
    "A": ["CA", "CB"],
    "V": ["CA", "CB", "CG1", "CG2"],
    "L": ["CA", "CB", "CG", "CD1", "CD2"],
    "I": ["CA", "CB", "CG1", "CG2", "CD1"],
    "F": ["CA", "CB", "CG", "CD1", "CE1", "CZ", "CE2", "CD2"],
    "W": ["CA", "CB", "CG", "CD1", "NE1", "CE2", "CZ2", "CH2", "CZ3", "CE3", "CD2"],
    "M": ["CA", "CB", "CG", "SD", "CE"],
    "P": ["CD", "CG", "CB", "CA"]
}

THREE_LETTER_CODE = {
    "A": "ALA", "V": "VAL", "L": "LEU", "I": "ILE",
    "F": "PHE", "W": "TRP", "M": "MET", "P": "PRO"
}

def read_atoms_from_prmtop(prmtop_file):
    """Read atom names from AMBER prmtop file."""
    print(f"Reading atoms from prmtop file: {prmtop_file}")
    atom_names = []
    atom_flag = False
    
    with open(prmtop_file, 'r') as f:
        for line in f:
            if not atom_flag:
                if line.startswith("%FLAG ATOM_NAME"):
                    atom_flag = True
            else:
                if line.startswith("%FORMAT"):
                    continue
                if line.startswith("%FLAG"):
                    break
                
                # Parse atom names (4 characters each)
                for i in range(0, len(line), 4):
                    if i + 4 > len(line):
                        break
                    atom_name = line[i:i+4].strip()
                    if atom_name:
                        atom_names.append(atom_name)
    
    return atom_names

def get_atom_index(atom_names, residue_index, atom_name, atom_count):
    """Get 1-based atom index for a specific residue and atom."""
    res_count = 0
    for i in range(atom_count):
        if atom_names[i] == "N":
            res_count += 1
            if atom_name == "N" and res_count == residue_index:
                return i + 1  # 1-based index
            
            if res_count == residue_index:
                for j in range(i + 1, atom_count):
                    if atom_names[j] == "N":
                        break
                    if atom_names[j] == atom_name:
                        return j + 1  # 1-based index
                break
    
    raise ValueError(f"Could not find the {residue_index}-th occurrence of atom '{atom_name}'")

def make_ss_groups(seq_file, ss_file):
    """
    Parse secondary structure file to find beta-strand segments (marked as 'E').
    Returns list of (start, end) tuples (1-based indexing).
    """
    with open(seq_file, 'r') as f:
        seq = f.readline().strip()
    
    with open(ss_file, 'r') as f:
        ss = f.readline().strip()
    
    starts = []
    ends = []
    active = 0
    extended = False
    start = 0
    
    for i, c in enumerate(ss):
        if c not in ['H', 'E', '.']:
            continue
        
        # Close a run if we leave 'E'
        if c != 'E' and extended:
            starts.append(start + 1)  # 1-based
            ends.append(i)  # 1-based (i is already the position after the last E)
            extended = False
        
        # When we see 'E', count it and possibly start a new run
        if c == 'E':
            active += 1
            if not extended:
                start = i
                extended = True
    
    # Close run if it extends to end of string
    if extended:
        starts.append(start + 1)
        ends.append(len(ss))
    
    return starts, ends, active

def write_restraints(seq_file, rst_file, index_file, atom_names, starts, ends):
    """Write restraint files for secondary structure and hydrophobic contacts."""
    atom_count = len(atom_names)
    
    with open(rst_file, 'w') as fo, open(index_file, 'w') as fi:
        # Strand-pair collection
        ss_group = 0
        n_segments = len(starts)
        
        if n_segments > 1:
            for i in range(n_segments - 1):
                for j in range(i + 1, n_segments):
                    # Generate all pairs between segments i and j
                    for res_i in range(starts[i], ends[i] + 1):
                        for res_j in range(starts[j], ends[j] + 1):
                            ss_group += 1
                            
                            # N to O restraint
                            atom_i = get_atom_index(atom_names, res_i, "N", atom_count)
                            atom_j = get_atom_index(atom_names, res_j, "O", atom_count)
                            fo.write(f" &rst\tiat = {atom_i}, {atom_j},\n")
                            fo.write(f"\tr1 = {1.0:.3f}, r2 = {2.0:.3f}, r3 = {4.0:.3f}, r4 = {5.0:.3f},\n")
                            fo.write(f"\trk2 = {100.0:8.3f}, rk3 = {100.0:8.3f},\n\t/\n")
                            fi.write(f"{1:2d} {ss_group:4d} {'SS':8s}\n")
                            
                            # O to N restraint
                            atom_i = get_atom_index(atom_names, res_i, "O", atom_count)
                            atom_j = get_atom_index(atom_names, res_j, "N", atom_count)
                            fo.write(f" &rst\tiat = {atom_i}, {atom_j},\n")
                            fo.write(f"\tr1 = {1.0:.3f}, r2 = {2.0:.3f}, r3 = {4.0:.3f}, r4 = {5.0:.3f},\n")
                            fo.write(f"\trk2 = {100.0:8.3f}, rk3 = {100.0:8.3f},\n\t/\n\n")
                            fi.write(f"{1:2d} {ss_group:4d} {'SS':8s}\n")
        else:
            print("No strand-pairing found; skipping strand-pair restraints.")
        
        # Hydrophobic collection
        hp_group = 0
        
        with open(seq_file, 'r') as f:
            seq = f.readline().strip()
        
        seq_len = len(seq)
        
        if seq_len > 8:
            for i in range(seq_len):
                if seq[i] in HP_ATOMS:
                    for j in range(i + 7, seq_len):  # At least 7 residues apart
                        if seq[j] in HP_ATOMS:
                            hp_group += 1
                            
                            atoms_i = HP_ATOMS[seq[i]]
                            atoms_j = HP_ATOMS[seq[j]]
                            
                            for atom_i in atoms_i:
                                for atom_j in atoms_j:
                                    idx_i = get_atom_index(atom_names, i + 1, atom_i, atom_count)
                                    idx_j = get_atom_index(atom_names, j + 1, atom_j, atom_count)
                                    
                                    fo.write(f" &rst\tiat = {idx_i}, {idx_j},\n")
                                    fo.write(f"\tr1 = {1.0:.3f}, r2 = {2.0:.3f}, r3 = {4.0:.3f}, r4 = {5.0:.3f},\n")
                                    fo.write(f"\trk2 = {100.0:8.3f}, rk3 = {100.0:8.3f},\n\t/\n")
                                    fi.write(f"{2:2d} {hp_group:4d} {'HP':8s}\n")
                            
                            fo.write("\n")
        else:
            print("Sequence length is less than 8 residues; skipping hydrophobic restraints.")

def main():
    parser = argparse.ArgumentParser(
        description="Generate AMBER restraint files for MELD protocol"
    )
    parser.add_argument('-p', '--prmtop', required=True, help="AMBER prmtop file")
    parser.add_argument('-s', '--sequence', required=True, help="Sequence file (1-letter code)")
    parser.add_argument('-ss', '--ss', required=True, help="Secondary structure file")
    parser.add_argument('-o', '--output', default='restraints', help="Output file prefix")
    
    args = parser.parse_args()
    
    rst_file = f"{args.output}.rst"
    index_file = f"{args.output}_index.txt"
    
    # Parse secondary structure
    starts, ends, active = make_ss_groups(args.sequence, args.ss)
    print(f"Found {len(starts)} beta-strand segments with {active} total residues")
    
    # Read atoms from prmtop
    atom_names = read_atoms_from_prmtop(args.prmtop)
    print(f"Read {len(atom_names)} atoms from prmtop")
    
    # Write restraints
    write_restraints(args.sequence, rst_file, index_file, atom_names, starts, ends)
    
    print(f"Restraints written to {rst_file}")
    print(f"Index written to {index_file}")

if __name__ == "__main__":
    main()