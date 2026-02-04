#!/usr/bin/env python3
"""
Generate AMBER restraint files for MELD protocol.
Creates distance restraints for secondary structure elements and hydrophobic interactions.
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
    
    Handles ACE/NHE caps properly. ACE cap (if present) is residue 0,
    first real amino acid is residue 1.
    """
    residue_starts = []
    
    # Check for ACE cap at the beginning (no N atom in first few atoms)
    has_ace = True
    for i in range(min(6, atom_count)):
        if atom_names[i] == "N":
            has_ace = False
            break
    
    ace_offset = 0
    if has_ace:
        for i in range(atom_count):
            if atom_names[i] == "N":
                residue_starts.append(0)
                ace_offset = 1
                break
    
    for i in range(atom_count):
        if atom_names[i] == "N":
            residue_starts.append(i)
    
    # Adjust residue index for ACE cap
    actual_index = residue_index - 1 + ace_offset
    
    if actual_index < 0 or actual_index >= len(residue_starts):
        raise ValueError(f"Residue index {residue_index} out of range (1-{len(residue_starts)-ace_offset})")
    
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
    """
    Parse secondary structure file to find beta-strand segments (marked as 'E').
    Returns list of (start, end) tuples (1-based indexing).
    """
    with open(seq_file, 'r') as f:
        lines = f.readlines()
        seq = lines[1].strip() if lines[0].startswith('>') else lines[0].strip()
    
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
        
        if c != 'E' and extended:
            starts.append(start + 1)  # 1-based
            ends.append(i)  # 1-based
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

def write_restraints(seq_file, rst_file, index_file, atom_names, starts, ends, 
                    cfrac_ss, cfrac_hp, gnum_ss, gnum_hp):
    """Write restraint files for secondary structure and hydrophobic contacts.
    
    Returns:
        tuple: (cfrac_hp_final, num_hp_groups, num_group_active_hp)
    """
    atom_count = len(atom_names)
    
    ss_restraints = []
    hp_restraints = []
    
    # Collect SS restraints
    n_segments = len(starts)
    ss_group = 0
    
    if n_segments > 1:
        for i in range(n_segments - 1):
            for j in range(i + 1, n_segments):
                for res_i in range(starts[i], ends[i] + 1):
                    for res_j in range(starts[j], ends[j] + 1):
                        ss_group += 1
                        
                        try:
                            # N to O restraint
                            atom_i = get_atom_index(atom_names, res_i, "N", atom_count)
                            atom_j = get_atom_index(atom_names, res_j, "O", atom_count)
                            ss_restraints.append((res_i, res_j, atom_i, atom_j, True, ss_group))
                            
                            # O to N restraint
                            atom_i = get_atom_index(atom_names, res_i, "O", atom_count)
                            atom_j = get_atom_index(atom_names, res_j, "N", atom_count)
                            ss_restraints.append((res_i, res_j, atom_i, atom_j, False, ss_group))
                        except ValueError as e:
                            print(f"Warning: Skipping SS restraint {res_i}-{res_j}: {e}")
                            continue
    
    # Collect HP restraints
    with open(seq_file, 'r') as f:
        lines = f.readlines()
        seq = lines[1].strip() if lines[0].startswith('>') else lines[0].strip()
    
    seq_len = len(seq)
    hp_group = 0
    
    if seq_len > 8:
        for i in range(seq_len):
            if seq[i] not in HP_ATOMS:
                continue
            
            for j in range(i + 7, seq_len):
                if seq[j] not in HP_ATOMS:
                    continue
                
                hp_group += 1
                atoms_i = HP_ATOMS[seq[i]]
                atoms_j = HP_ATOMS[seq[j]]
                
                for atom_i in atoms_i:
                    for atom_j in atoms_j:
                        try:
                            idx_i = get_atom_index(atom_names, i + 1, atom_i, atom_count)
                            idx_j = get_atom_index(atom_names, j + 1, atom_j, atom_count)
                            hp_restraints.append((i, j, atom_i, atom_j, idx_i, idx_j, hp_group))
                        except ValueError:
                            continue
    
    # Calculate statistics
    num_collections = 2 if (len(ss_restraints) > 0 and len(hp_restraints) > 0) else 1
    num_ss_groups = ss_group
    num_hp_groups = hp_group
    
    # Calculate number of group-active restraints after Step 1
    num_group_active_ss = num_ss_groups * gnum_ss
    num_group_active_hp = num_hp_groups * gnum_hp
    
    print(f"\nGroup-level activation:")
    print(f"  SS: {num_ss_groups} groups × {gnum_ss} = {num_group_active_ss} group-active restraints")
    print(f"  HP: {num_hp_groups} groups × {gnum_hp} = {num_group_active_hp} group-active restraints")
    
    # Validate and calculate cfrac_hp
    hp_count, _ = count_hydrophobic_residues(seq_file)
    target_hp_active = 1.2 * hp_count
    
    print(f"\nTarget HP restraints: {target_hp_active:.1f} (1.2 × {hp_count} HP residues)")
    
    # Validation
    if target_hp_active > num_group_active_hp:
        raise ValueError(
            f"\n{'='*60}\n"
            f"ERROR: Invalid HP activation parameters!\n"
            f"{'='*60}\n"
            f"Target HP restraints: {target_hp_active:.1f}\n"
            f"Group-active HP restraints: {num_group_active_hp}\n"
            f"\nThe target exceeds available group-active restraints!\n"
            f"\nSolutions:\n"
            f"  1. Increase --gnum-hp (currently {gnum_hp})\n"
            f"     Minimum required: {int(target_hp_active / num_hp_groups) + 1}\n"
            f"  2. Manually set --cfrac-hp to a lower value\n"
            f"     Maximum possible: {num_group_active_hp / hp_count:.2f}\n"
            f"{'='*60}"
        )
    
    if cfrac_hp is None:
        cfrac_hp = target_hp_active / num_group_active_hp
        
        # Clamp to [0.0, 1.0]
        if cfrac_hp > 1.0:
            print(f"Warning: Calculated cfrac_hp={cfrac_hp:.2f} > 1.00, clamping to 1.00")
            cfrac_hp = 1.00
        
        print(f"Auto-calculated cfrac_hp: {cfrac_hp:.2f}")
    
    with open(rst_file, 'w') as fo, open(index_file, 'w') as fi:
        # Line 1: Number of collections
        fi.write(f"{num_collections}\n")
        
        # Line 2: Number of groups per collection
        if num_collections == 2:
            fi.write(f"{num_ss_groups},{num_hp_groups}\n")
        elif len(ss_restraints) > 0:
            fi.write(f"{num_ss_groups}\n")
        else:
            fi.write(f"{num_hp_groups}\n")
        
        # Line 3: Number to activate per group
        if num_collections == 2:
            fi.write(f"{gnum_ss},{gnum_hp}\n")
        elif len(ss_restraints) > 0:
            fi.write(f"{gnum_ss}\n")
        else:
            fi.write(f"{gnum_hp}\n")
        
        # Line 4: Collection fractions
        if num_collections == 2:
            fi.write(f"{cfrac_ss:.2f},{cfrac_hp:.2f}\n")
        elif len(ss_restraints) > 0:
            fi.write(f"{cfrac_ss:.2f}\n")
        else:
            fi.write(f"{cfrac_hp:.2f}\n")
        
        # Write SS restraints
        for res_i, res_j, atom_i, atom_j, is_N_to_O, group_id in ss_restraints:
            fo.write(f" &rst\tiat = {atom_i}, {atom_j},\n")
            fo.write(f"\tr1 = {1.0:.3f}, r2 = {2.0:.3f}, r3 = {4.0:.3f}, r4 = {5.0:.3f},\n")
            fo.write(f"\trk2 = {0.6:8.3f}, rk3 = {0.6:8.3f},\n\t/\n")
            if is_N_to_O:
                fo.write("\n")
            fi.write(f"{1:2d} {group_id:4d} {'SS':8s}\n")
        
        # Write HP restraints
        current_group = None
        for i, j, atom_i, atom_j, idx_i, idx_j, group_id in hp_restraints:
            if current_group is not None and group_id != current_group:
                fo.write("\n")
            
            fo.write(f" &rst\tiat = {idx_i}, {idx_j},\n")
            fo.write(f"\tr1 = {1.0:.3f}, r2 = {2.0:.3f}, r3 = {4.0:.3f}, r4 = {5.0:.3f},\n")
            fo.write(f"\trk2 = {0.6:8.3f}, rk3 = {0.6:8.3f},\n\t/\n")
            fi.write(f"{2:2d} {group_id:4d} {'HP':8s}\n")
            
            current_group = group_id
    
    return cfrac_hp, num_hp_groups, num_group_active_hp

def main():
    parser = argparse.ArgumentParser(
        description="Generate AMBER restraint files for MELD protocol",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -p system.prmtop -s sequence.fasta -ss secondary.dat
  %(prog)s -p system.prmtop -s sequence.fasta -ss secondary.dat --gnum-hp 2
  
Default values:
  gnum-ss: 1 (activate 1 restraint per SS group)
  gnum-hp: 1 (activate 1 restraint per HP group)
  cfrac-ss: 0.45 (activate 45%% of group-active SS restraints)
  cfrac-hp: (1.2 × num_HP_residues) / (num_HP_groups × gnum_hp)
        """
    )
    parser.add_argument('-p', '--prmtop', required=True,
                        help="AMBER prmtop file")
    parser.add_argument('-s', '--sequence', required=True,
                        help="Sequence file (FASTA format, 1-letter code)")
    parser.add_argument('-ss', '--ss', required=True,
                        help="Secondary structure file (H/E/. notation)")
    parser.add_argument('-o', '--output', default='restraints',
                        help="Output file prefix (default: restraints)")
    parser.add_argument('--gnum-ss', type=int, default=1,
                        help="Number of restraints to activate per SS group (default: 1)")
    parser.add_argument('--gnum-hp', type=int, default=1,
                        help="Number of restraints to activate per HP group (default: 1)")
    parser.add_argument('--cfrac-ss', type=float, default=0.45,
                        help="Collection fraction for SS restraints (default: 0.45)")
    parser.add_argument('--cfrac-hp', type=float, default=None,
                        help="Collection fraction for HP restraints (default: auto-calculated)")
    
    args = parser.parse_args()
    
    rst_file = f"{args.output}.rst"
    index_file = f"{args.output}_index.txt"
    
    # Count hydrophobic residues
    hp_count, seq_len = count_hydrophobic_residues(args.sequence)
    print(f"Sequence length: {seq_len}")
    print(f"Hydrophobic residues: {hp_count} ({hp_count/seq_len*100:.1f}%)")
    
    # Parse secondary structure
    starts, ends, active = make_ss_groups(args.sequence, args.ss)
    print(f"Found {len(starts)} beta-strand segments with {active} total residues")
    
    # Read atoms from prmtop
    atom_names = read_atoms_from_prmtop(args.prmtop)
    print(f"Read {len(atom_names)} atoms from prmtop")
    
    # Check for ACE cap
    has_ace = True
    for i in range(min(6, len(atom_names))):
        if atom_names[i] == "N":
            has_ace = False
            break
    if has_ace:
        print("Detected ACE cap at N-terminus")
    
    # Write restraints (will calculate and validate cfrac_hp)
    try:
        cfrac_hp_final, num_hp_groups, num_group_active_hp = write_restraints(
            args.sequence, rst_file, index_file, atom_names, 
            starts, ends, args.cfrac_ss, args.cfrac_hp, args.gnum_ss, args.gnum_hp)
        
        print(f"\nRestraints written to {rst_file}")
        print(f"Index written to {index_file}")
        print(f"\nActivation parameters:")
        print(f"  Group level:")
        print(f"    SS: {args.gnum_ss} restraint(s) per group")
        print(f"    HP: {args.gnum_hp} restraint(s) per group")
        print(f"  Collection level:")
        print(f"    SS: cfrac={args.cfrac_ss:.2f}")
        print(f"    HP: cfrac={cfrac_hp_final:.2f}")
        
        return 0
        
    except ValueError as e:
        print(f"\n{e}")
        print("\nRestraint generation failed. No files written.")
        return 1

if __name__ == "__main__":
    exit(main())