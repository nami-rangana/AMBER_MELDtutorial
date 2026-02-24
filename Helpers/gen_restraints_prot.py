#!/usr/bin/env python3
"""
Generate AMBER restraint files for MELD protocol.
Creates restraints for:
1. Beta-sheets (SS) - Distance restraints (N-O, O-N)
2. Helices (HELIX) - Torsion (Phi/Psi) and Distance (CA-CA) restraints
3. Hydrophobic interactions (HP) - Distance restraints
"""
import argparse

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

def count_hydrophobic_residues(seq_file):
    """Count the number of hydrophobic residues in the sequence."""
    with open(seq_file, 'r') as f:
        lines = f.readlines()
        seq = lines[1].strip() if lines[0].startswith('>') else lines[0].strip()
    hp_count = sum(1 for aa in seq if aa in HP_ATOMS)
    return hp_count, len(seq)

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
    """Get 1-based atom index for a specific residue and atom.
    Handles ACE/NHE caps properly.
    """
    residue_starts = []
    # Check for ACE cap
    has_ace = True
    for i in range(min(6, atom_count)):
        if atom_names[i] == "N":
            has_ace = False
            break
    ace_offset = 0
    if has_ace:
        ace_offset = 1
        residue_starts.append(0) # Placeholder for ACE start

    for i in range(atom_count):
        if atom_names[i] == "N":
            residue_starts.append(i)
            
    # Adjust residue index (1-based input -> 0-based list index)
    actual_index = residue_index - 1 + ace_offset
    
    if actual_index >= len(residue_starts):
         raise ValueError(f"Residue index {residue_index} out of range")

    start_idx = residue_starts[actual_index]
    if actual_index + 1 < len(residue_starts):
        end_idx = residue_starts[actual_index + 1]
    else:
        end_idx = atom_count

    for i in range(start_idx, end_idx):
        if atom_names[i] == atom_name:
            return i + 1
    raise ValueError(f"Could not find atom '{atom_name}' in residue {residue_index}")

def make_ss_groups(seq_file, ss_file):
    """Parse secondary structure file to find beta-strand segments (E)."""
    with open(seq_file, 'r') as f:
        lines = f.readlines()
    with open(ss_file, 'r') as f:
        ss = f.readline().strip()
    starts = []
    ends = []
    active = 0
    extended = False
    start = 0
    for i, c in enumerate(ss):
        if c not in ['H', 'E', '.']: continue
        if c != 'E' and extended:
            starts.append(start + 1)
            ends.append(i)
            extended = False
        if c == 'E':
            active += 1
            if not extended:
                start = i
                extended = True
    if extended:
        starts.append(start + 1)
        ends.append(len(ss))
    return starts, ends, active

def make_helix_groups(seq_file, ss_file, run_length=5, min_match=4):
    """
    Parse secondary structure file to find Helix segments (marked as 'H').
    Uses the sliding window approach from parse.py.
    Returns list of (start, end) tuples (0-based indexing).
    """
    with open(ss_file, 'r') as f:
        ss_content = f.readline().strip()
        
    # Mark elements that are 'H'
    has_correct_type = [1 if ss == 'H' else 0 for ss in ss_content]
    length = len(ss_content)
    
    if length < run_length:
        return []

    # Sliding window sum
    totals = [0] * (length - run_length + 1)
    for index in range(length - run_length + 1):
        totals[index] = sum(has_correct_type[index : index + run_length])
        
    results = []
    for index in range(len(totals)):
        if totals[index] >= min_match:
            # Store as (start, end) indices
            results.append((index, index + run_length))
            
    return results

def write_restraints(seq_file, rst_file, index_file, atom_names, 
                    starts_ss, ends_ss, helix_windows,
                    cfrac_ss, cfrac_hp, cfrac_helix,
                    gnum_ss, gnum_hp, gnum_helix):
    
    atom_count = len(atom_names)
    ss_restraints = []
    hp_restraints = []
    helix_restraints = []

    # --- 1. Beta-Sheet (SS) Restraints ---
    n_segments = len(starts_ss)
    ss_group_count = 0
    if n_segments > 1:
        for i in range(n_segments - 1):
            for j in range(i + 1, n_segments):
                for res_i in range(starts_ss[i], ends_ss[i] + 1):
                    for res_j in range(starts_ss[j], ends_ss[j] + 1):
                        ss_group_count += 1
                        try:
                            # N to O
                            ai = get_atom_index(atom_names, res_i, "N", atom_count)
                            aj = get_atom_index(atom_names, res_j, "O", atom_count)
                            ss_restraints.append((ai, aj, ss_group_count))
                            # O to N
                            ai = get_atom_index(atom_names, res_i, "O", atom_count)
                            aj = get_atom_index(atom_names, res_j, "N", atom_count)
                            ss_restraints.append((ai, aj, ss_group_count))
                        except ValueError: continue

    # --- 2. Helix Restraints ---
    # Derived from parse.py [1]:
    # Torsions: 2.5e-2 kJ/mol/deg^2 force constant.
    # Distances: 2500 kJ/mol/nm^2 force constant.
    helix_group_count = 0
    
    for start, end in helix_windows:
        helix_group_count += 1
        
        # A. Torsions (Phi/Psi)
        # Loop indices: start+1 to end-1 (non-inclusive)
        for i in range(start + 1, end - 1):
            res_idx = i + 1  # Current residue (1-based)
            res_prev = i     # Previous residue (1-based)
            res_next = i + 2 # Next residue (1-based)
            
            try:
                # Phi Restraint: C(prev) - N - CA - C
                # Target: -62.5 deg (range -122.5 to -2.5)
                at1 = get_atom_index(atom_names, res_prev, "C", atom_count)
                at2 = get_atom_index(atom_names, res_idx, "N", atom_count)
                at3 = get_atom_index(atom_names, res_idx, "CA", atom_count)
                at4 = get_atom_index(atom_names, res_idx, "C", atom_count)
                
                # Format: (TYPE, a1, a2, a3, a4, r1, r2, r3, r4, k2, k3, group)
                # Using r2/r3 = target +/- width. 
                # parse.py uses broad flat wells: target - 60 to target + 60 usually for min/max
                # Here we map: r2=min, r3=max. 
                # Values: -62.5 +/- 20 seems tight? 
                # parse.py inputs: -62.5 (center), 17.5 (width?) -> [-80, -45]? 
                # Let's use standard broad helix definitions:
                helix_restraints.append(('TORSION', at1, at2, at3, at4, -122.5, -82.5, -42.5, -2.5, 20.0, 20.0, helix_group_count)) # Degrees, 20 kcal/mol/rad^2 ** how amber defines torsion force constants

                # Psi Restraint: N - CA - C - N(next)
                # Target: -42.5 deg
                at1 = get_atom_index(atom_names, res_idx, "N", atom_count)
                at2 = get_atom_index(atom_names, res_idx, "CA", atom_count)
                at3 = get_atom_index(atom_names, res_idx, "C", atom_count)
                at4 = get_atom_index(atom_names, res_next, "N", atom_count)
                helix_restraints.append(('TORSION', at1, at2, at3, at4, -102.5, -62.5, -22.5, 17.5, 20.0, 20.0, helix_group_count))
            except ValueError: pass

        # B. Distances (CA-CA)
        # parse.py uses nm. AMBER uses Angstroms. We multiply by 10.
        try:
            r_start = start + 1 # 1-based start index
            
            # d1: i to i+3
            # 0.485 nm -> 4.85 A. quadratic cut 0.2nm -> 2.0 A
            at1 = get_atom_index(atom_names, r_start, "CA", atom_count)
            at2 = get_atom_index(atom_names, r_start+3, "CA", atom_count)
            helix_restraints.append(('DIST', at1, at2, 0.0, 4.85, 5.61, 7.61, 0.6, 0.6, helix_group_count)) # Angstroms, 0.6 kcal/mol/A^2

            # d2: i+1 to i+4
            at1 = get_atom_index(atom_names, r_start+1, "CA", atom_count)
            at2 = get_atom_index(atom_names, r_start+4, "CA", atom_count)
            helix_restraints.append(('DIST', at1, at2, 0.0, 4.85, 5.61, 7.61, 0.6, 0.6, helix_group_count))

            # d3: i to i+4
            at1 = get_atom_index(atom_names, r_start, "CA", atom_count)
            at2 = get_atom_index(atom_names, r_start+4, "CA", atom_count)
            helix_restraints.append(('DIST', at1, at2, 0.0, 5.81, 6.84, 8.84, 0.6, 0.6, helix_group_count))
            
        except ValueError: pass

    # --- 3. Hydrophobic (HP) Restraints ---
    with open(seq_file, 'r') as f:
        lines = f.readlines()
        seq = lines[1].strip() if lines[0].startswith('>') else lines[0].strip()
    seq_len = len(seq)
    hp_group_count = 0
    if seq_len > 8:
        for i in range(seq_len):
            if seq[i] not in HP_ATOMS: continue
            for j in range(i + 7, seq_len):
                if seq[j] not in HP_ATOMS: continue
                hp_group_count += 1
                atoms_i = HP_ATOMS[seq[i]]
                atoms_j = HP_ATOMS[seq[j]]
                for atom_i in atoms_i:
                    for atom_j in atoms_j:
                        try:
                            ai = get_atom_index(atom_names, i + 1, atom_i, atom_count)
                            aj = get_atom_index(atom_names, j + 1, atom_j, atom_count)
                            hp_restraints.append((ai, aj, hp_group_count))
                        except ValueError: continue

    # --- Write Files ---
    collections = []
    if len(ss_restraints) > 0: collections.append('SS')
    if len(helix_restraints) > 0: collections.append('HELIX')
    if len(hp_restraints) > 0: collections.append('HP')
    
    num_collections = len(collections)
    
    with open(rst_file, 'w') as fo, open(index_file, 'w') as fi:
        def write_rst_entry(atoms, r1, r2, r3, r4, k2, k3):
            # Join atoms with commas
            atom_str = ",".join(map(str, atoms))
            fo.write(f" &rst\n\tiat = {atom_str},\n")
            fo.write(f"\tr1 = {r1:.3f}, r2 = {r2:.3f}, r3 = {r3:.3f}, r4 = {r4:.3f},\n")
            fo.write(f"\trk2 = {k2:8.3f}, rk3 = {k3:8.3f},\n\t/\n")

        # Header Line 1: Num collections
        fi.write(f"{num_collections}\n")
        
        # Header Line 2: Groups per collection
        counts = []
        if 'SS' in collections: counts.append(str(ss_group_count))
        if 'HELIX' in collections: counts.append(str(helix_group_count))
        if 'HP' in collections: counts.append(str(hp_group_count))
        fi.write(",".join(counts) + "\n")
        
        # Header Line 3: Active per group
        gnums = []
        if 'SS' in collections: gnums.append(str(gnum_ss))
        if 'HELIX' in collections: gnums.append(str(gnum_helix))
        if 'HP' in collections: gnums.append(str(gnum_hp))
        fi.write(",".join(gnums) + "\n")
        
        # Header Line 4: Fractions
        cfracs = []
        if 'SS' in collections: cfracs.append(f"{cfrac_ss:.2f}")
        if 'HELIX' in collections: cfracs.append(f"{cfrac_helix:.2f}")
        if 'HP' in collections:
            hp_count, _ = count_hydrophobic_residues(seq_file)
            target = 1.2 * hp_count
            calc_cfrac = target / (hp_group_count * gnum_hp) if hp_group_count > 0 else 0
            if cfrac_hp is not None: calc_cfrac = cfrac_hp
            if calc_cfrac > 1.0: calc_cfrac = 1.0
            cfracs.append(f"{calc_cfrac:.2f}")
        fi.write(",".join(cfracs) + "\n")

        # Write Data
        col_idx = 1
        
        # SS Data
        if 'SS' in collections:
            for at1, at2, grp in ss_restraints:
                write_rst_entry([at1, at2], 0.0, 0.0, 3.5, 5.5, 0.6, 0.6) # Angstroms, 0.6 kcal/mol/A^2
                fi.write(f"{col_idx:2d} {grp:4d} SS\n")
            col_idx += 1
            fo.write("#\n#\n") # collection breaker
            
        # HELIX Data
        if 'HELIX' in collections:
            for r in helix_restraints:
                rtype = r[0]
                # Dynamic unpacking
                if rtype == 'TORSION':
                    at1, at2, at3, at4 = r[1], r[2], r[3], r[4]
                    r1, r2, r3, r4 = r[5], r[6], r[7], r[8]
                    k2, k3 = r[9], r[10]
                    grp = r[11]
                    write_rst_entry([at1, at2, at3, at4], r1, r2, r3, r4, k2, k3)
                else: # DIST
                    at1, at2 = r[1], r[2]
                    r1, r2, r3, r4 = r[3], r[4], r[5], r[6]
                    k2, k3 = r[7], r[8]
                    grp = r[9]
                    write_rst_entry([at1, at2], r1, r2, r3, r4, k2, k3)
                
                fi.write(f"{col_idx:2d} {grp:4d} HX\n")
            col_idx += 1
            fo.write("#\n#\n") # collection breaker
            
        # HP Data
        if 'HP' in collections:
            for at1, at2, grp in hp_restraints:
                 write_rst_entry([at1, at2], 0.0, 0.0, 5.0, 7.0, 0.6, 0.6) # Angstroms, 0.6 kcal/mol/A^2
                 fi.write(f"{col_idx:2d} {grp:4d} HP\n")
            fo.write("#\n#\n") # collection breaker

    return 0

def main():
    parser = argparse.ArgumentParser(description="Generate AMBER restraints for MELD")
    parser.add_argument('-p', '--prmtop', required=True, help="AMBER prmtop file")
    parser.add_argument('-s', '--sequence', required=True, help="Sequence file (FASTA)")
    parser.add_argument('-ss', '--ss', required=True, help="Secondary structure file (H/E)")
    parser.add_argument('-o', '--output', default='restraints', help="Output prefix")
    
    parser.add_argument('--gnum-ss', type=int, default=1)
    parser.add_argument('--gnum-hp', type=int, default=1)
    parser.add_argument('--gnum-helix', type=int, default=7, help="Restraints per Helix group (default 7)")
    
    parser.add_argument('--cfrac-ss', type=float, default=0.45)
    parser.add_argument('--cfrac-hp', type=float, default=None)
    parser.add_argument('--cfrac-helix', type=float, default=0.8, help="Fraction for Helix (default 0.8)")
    
    args = parser.parse_args()
    
    # 1. Parse SS (Beta)
    starts, ends, active = make_ss_groups(args.sequence, args.ss)
    
    # 2. Parse Helix
    helix_windows = make_helix_groups(args.sequence, args.ss)
    
    # 3. Read Atoms
    atom_names = read_atoms_from_prmtop(args.prmtop)
    
    # 4. Write
    try:
        write_restraints(args.sequence, f"{args.output}.rst", f"{args.output}_index.txt", atom_names,
                        starts, ends, helix_windows,
                        args.cfrac_ss, args.cfrac_hp, args.cfrac_helix,
                        args.gnum_ss, args.gnum_hp, args.gnum_helix)
        print("Done.")
    except Exception as e:
        print(f"Error: {e}")
        return 1
    return 0

if __name__ == "__main__":
    main()
    