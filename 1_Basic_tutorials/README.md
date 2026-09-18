# Creating **MELD** `INDXF` from scratch:

## Learning Outcomes:
* Be able to write an **AMBER** compatible MELD `INDXF` file and its **AMBER** restraint file `DISANG`
* Be able to explain what a MELD **restraint**, **group** and **collection** are, and decide what belongs in each
* Be able to build a hydrophobic-core collection for a peptide from its sequence alone, and check it before spending GPU time on it

## Introduction:
In this tutorial, we use chignolin (PDB ID: 1UAO), a small protein made up of just 10 amino acids. It folds into a compact β-hairpin, with two short antiparallel β-strands connected by a tight turn. Although chignolin is much smaller than most proteins, it captures the basic features of β-sheet folding and is therefore commonly used to test protein-folding simulations.

Chignolin’s small size and simple topology make it a convenient model for understanding how MELD restraints are organized. The peptide contains three residues with substantial nonpolar character—**Tyr2, Pro4, and Trp9**—whose side chains contribute to the compact cluster that stabilizes the folded β-hairpin. Because this hydrophobic core involves only a few residues, the complete restraint hierarchy can be constructed manually and inspected one restraint at a time. In this tutorial, we therefore focus on a single MELD collection: the hydrophobic-core collection.

> *For a complete walkthrough of building and using CPI restraint collections, see the [protein folding tutorial](../2_Protein_Folding/README.md).*


## Theory: MELD Restraints, Groups, and Collections

MELD combines an atomistic force field with structural information that may be incomplete, ambiguous, or partly incorrect. In Bayesian terms the AMBER force field supplies the prior and the restraints supply the likelihood of the data given a structure.<sup>[1]</sup>

Rather than enforcing every restraint simultaneously, MELD treats the information as a hierarchy and activates only the subset that is most compatible with the current structure. This allows several possible interpretations of the input data to create different low-energy regions in the conformational landscape.

A MELD hierarchy has three levels:

| Level | Description |
| :-- | :-- |
| **Restraint** | A single geometrical condition, such as a distance or torsion that is flexibly restricted rather than locked in place. <br> *ex:* the distance CZ(Tyr2)–CZ2(Trp9), free below 5 Å and penalized beyond it. |
| **Group** | A set of related restraints representing **one structural statement**, together with a `KEEP` count saying how much of that set has to be true for the statement to count as satisfied. <br> *ex:* all 88 Tyr2–Trp9 atom–atom distances. The statement is "Tyr2 touches Trp9", and *any one* of the 88 contacts makes it true — so `KEEP 1`. |
| **Collection** | A set of groups derived from the same type of information, ranked against each other, together with an `ACTIVE` count saying how many of them are enforced at any instant. <br> *ex:* the three residue-pair groups Tyr2–Pro4, Tyr2–Trp9 and Pro4–Trp9. They are all "hydrophobic contact" statements, so they are comparable, and we enforce the best-satisfied `ACTIVE` of them. |

For the collection built in this tutorial, the two levels stack like this:

```
hy.disang -- 164 &rst records, in a fixed order
      |
      |  GROUPS 32 88 44   sliced by position, in DISANG order
      v
   group 1  Tyr2-Pro4   restraints   1 -  32
   group 2  Tyr2-Trp9   restraints  33 - 120
   group 3  Pro4-Trp9   restraints 121 - 164
      |
      |  KEEP 1            each group keeps its lowest-energy restraint
      v
   3 candidate restraints, one per group
      |
      |  ACTIVE 2          the collection keeps its 2 lowest-energy groups
      v
   2 restraints exert force this timestep
```

### Where the groups come from


Nothing in MELD decides this for you. A group is a modelling decision, and the decision is always the same question: **what is the ambiguity I cannot resolve?** Whatever you cannot resolve goes *inside* a group, and `KEEP` says how much of it you expect to be true.
 
For a hydrophobic contact the ambiguity is atomic. The statement "Tyr2 and Trp9 are in the same hydrophobic core" is about two side chains, not two atoms. Tyr2 contributes 8 apolar atoms, Trp9 contributes 11, and any of the 8 × 11 = 88 pairings would satisfy the statement equally well — a ring stacked face-to-face and a ring packed edge-on are both packed. So all 88 distances go into one group with `KEEP 1`, and MELD switches on whichever single distance the current structure has come closest to satisfying. As the chain moves, the kept restraint changes; the *statement* does not.
 
Write those same 88 distances as 88 groups of one instead and you have said something entirely different. The collection would now rank *atom pairs* against each other rather than side chains, and its `ACTIVE` count would be a statement about how many individual atoms touch — at `ACTIVE all`, that all 88 pairs are within 5 Å at once, which no two side chains can manage. The group is what keeps the ambiguity from leaking up into the level above it.

### Where the collections come from

A collection is the second decision: **which statements am I willing to rank against each other, and how many of them do I actually believe?**

Ranking only means something if the groups are comparable — same kind of data, same shape, same number of kept restraints. The three chignolin pairs qualify: each is one hydrophobic contact, each keeps exactly one distance. `ACTIVE` then says how much of the collection survives the ranking. `ACTIVE all` is ordinary restrained MD — every statement enforced, no tolerance for a wrong one. Anything less buys tolerance: the worst-satisfied groups are switched off completely, and if one of your contacts is wrong, that is the one MELD will drop.

At every force evaluation MELD resolves the two levels in order:

### How the two levels resolved:
 
At every force evaluation, MELD resolves the hierarchy in order:
 
1. **Inside each group**, the `KEEP` lowest-energy restraints are kept; the group's energy is the **sum** of those kept.
2. **Inside each collection**, the `ACTIVE` lowest-energy groups are kept, and only their kept restraints exert any force. A group the collection drops contributes nothing at all.

The summation in step 1 explains why groups within the same collection should use a consistent `KEEP` value. Although our three groups contain 32, 88, and 44 candidate restraints, each group uses `KEEP 1`. MELD therefore evaluates each group using only the energy of its best-satisfied restraint, so every group contributes one restraint energy to the collection-level ranking. This selection is repeated as the structure changes during the simulation [1].

If we used `KEEP 2`, the energy assigned to each group would instead be the sum of its two best-satisfied restraints. The ranking could then become more sensitive to differences in how the groups were constructed—for example, to the number of atom-pair distances generated by side chains of different sizes. Using `KEEP 1` gives each hydrophobic-contact group the same basic interpretation: **find the single atom pair that best represents the side-chain contact**.

> **`ALWAYS`.** MELD also accepts restraints that bypass both levels and stay on for the whole run — a disulfide, a covalent link, a metal site, things chemistry guarantees. A hydrophobic contact is never one of them, and chignolin has none, so no `ALWAYS` record appears in this tutorial.

## Method
### 1. System setup

This is the earliest stage where we create the topology and the initial coordinate file for the intended system 1UAO. As in the folding tutorial, MELD begins from an extended chain and lets the restraints and the force field guide the fold, so the starting geometry carries no structural information.
 
For this tutorial we will generate the initial peptide chain using the sequence alone. Use the following command to download the FASTA entry for 1UAO from the PDB.

```bash
curl -o 1uao.fasta https://www.rcsb.org/fasta/entry/1UAO
```

> **Note:** Only the amino-acid sequence is taken from the PDB entry. The deposited coordinates are not used to build the starting structure. They are kept separate until the final validation stage, when the structures produced by MELD can be compared with the experimental structure.

#### Build the initial chain

The peptide is built with `tleap` using the `ff19SB` protein force field. We also select `mbondi3` atomic radii because the simulations in this tutorial use the `igb=8` generalized Born implicit-solvent model.

Instead of writing the complete `tleap` input manually, the helper script [makeLeap.x](1_system_setup/makeLeap.x) reads the one-letter FASTA sequence and converts it into the three-letter residue names expected by `tleap`.

Make the script executable and run it for 1UAO:

```bash
chmod +x makeLeap.x
./makeLeap.x -s 1uao.fasta -o 1uao
```

This creates [leap.in](1_system_setup/leap.in). Execute it, capturing the output verbose so that nothing is missed:

```bash
tleap -f leap.in > leap.log 2>&1
```

Before continuing, inspect `leap.log` carefully and confirm that the system was built without errors. Chignolin contains two acidic residues, `Asp3` and `Glu5`, so `tleap` may report a nonzero net charge. This is not necessarily a problem, particularly for an implicit-solvent simulation. However, the reported charge should match the protonation states and terminal groups defined in `leap.in`.

Warnings about unrecognized residues, missing atoms, or unavailable force-field parameters must be resolved before proceeding, as they may prevent the simulation from running correctly. A successful `tleap` run produces two files: [1uao.prmtop](1_system_setup/1uao.prmtop), which stores the topology, and [1uao.inpcrd](1_system_setup/1uao.inpcrd), which contains the coordinates of the initial extended structure.


#### Apply Hydrogen Mass Repartitioning

Next we apply Hydrogen Mass Repartitioning (HMR)<sup>[[3](https://doi.org/10.1021/ct5010406)]</sup>. HMR shifts mass from heavy atoms to the hydrogens bonded to them, slowing the fastest bond vibrations without changing the total mass or thermodynamics of the system. This allows a longer timestep of about 4 fs instead of the usual 2 fs. The repartitioning is handled by `ParmED`, which ships with AMBER, using the input file [HMR.in](1_system_setup/HMR.in):

HMR is applied with `ParmED`, which is included with AMBER. Create [HMR.in](1_system_setup/HMR.in) with the following contents:

```text
parm 1uao.prmtop
hmassrepartition
outparm 1uao_HMR.prmtop
quit
```

To execute simply use command:

```
parmed -i HMR.in
```

This creates [1uao_HMR.prmtop](1_system_setup/1uao_HMR.prmtop), which is used for the remainder of the tutorial. Note that the coordinate file has not changed: only the masses have been redistributed.

#### Minimize the starting structure

Finally, we minimize the initial structure. The purpose of this step is not to move chignolin toward its native fold or introduce any structural bias. Instead, minimization removes steric clashes and relaxes strained local geometries introduced when `tleap` builds the peptide from its sequence. This provides a stable starting point for molecular dynamics, while MELD subsequently uses the physical model and supplied restraints to guide conformational sampling. A short minimization in implicit solvent is sufficient for this purpose.


```text
energy minimization
 &cntrl
  imin = 1, maxcyc = 2000, ncyc = 500,
  ntwr = 1000, ntpr = 100,
  cut = 999.0, rgbmax = 999.0,
  ntb = 0, igb = 8,
  saltcon = 0.0,
 /
```

he above settings follow from implicit-solvent treatment. The flag `ntb = 0` removes periodic boundaries while `igb = 8` selects the generalized Born model, so there is no solvent box and no cutoff is needed (`cut = 999.0`, `rgbmax = 999.0`). Keep `igb = 8` here: it is the model matched to the `mbondi3` radii used to build the topology, and it is the model used in production simulation.

Run the minimization with:

```
srun $AMBERHOME/bin/pmemd -O \
    -i min.in \
    -p 1uao_HMR.prmtop \
    -c 1uao.inpcrd \
    -o min.out \
    -r min.rst \
    -inf min.info
```

Before proceeding, examine `min.out` and confirm that the minimization completed successfully. The energy should decrease steadily and begin to level off, with no `NaN` values or numerical errors. The resulting restart file, [min.rst](1_system_setup/min.rst), contains the minimized coordinates that will be used to initialize every replica in the MELD simulation.

It is good practice to inspect this structure visually before moving on to the more computationally expensive stages. At this point, chignolin should remain extended and unfolded, but its local geometry should be well relaxed, with no obvious steric clashes or overlapping atoms. It should not yet resemble the native β-hairpin, because the purpose of minimization is only to prepare a stable starting structure—not to fold the peptide.

### 2. Choosing the hydrophobic pairs

Nothing decides this step for you. A collection begins as a list of statements you are willing to make about the structure, and for chignolin we make three:

| | Contact | The statement |
| :-- | :-- | :-- |
| group 1 | Tyr2 – Pro4 | the tyrosine ring and the proline ring are packed against each other |
| group 2 | Tyr2 – Trp9 | the tyrosine ring and the indole are packed against each other |
| group 3 | Pro4 – Trp9 | the proline ring and the indole are packed against each other |

Three statements, three groups. That is the entire input to this tutorial, and where it comes from is a modelling judgement: chignolin's hairpin buries these three side chains against one another, and they are the only residues in the sequence with a side chain large and apolar enough to be worth restraining. Someone else might make two statements, or five. The machinery below does not care — it turns whatever list you write into restraints.

A contact is between *side chains*, and a restraint is between *atoms*, so each statement needs a set of atoms to stand for each residue. We take the apolar heavy atoms of each side chain:

| Residue | Atoms | *n* | prmtop atom numbers |
| :-- | :-- | :-- | :-- |
| Tyr2 | `CA CB CG CD1 CE1 CZ CE2 CD2` | 8 | 10, 12, 15, 16, 18, 20, 23, 25 |
| Pro4 | `CD CG CB CA` | 4 | 42, 45, 48, 51 |
| Trp9 | `CA CB CG CD1 NE1 CE2 CZ2 CH2 CZ3 CE3 CD2` | 11 | 107, 109, 112, 113, 115, 117, 118, 120, 122, 124, 126 |

#### Every atom pair becomes one restraint

Each statement is turned into restraints by pairing every atom of one residue with every atom of the other, because we are not claiming to know *which* atoms are nearby:

| Group | Contact | Atom pairs | Restraints | `DISANG` range|
| :-- | :-- | :-- | :-- | :-- |
| 1 | `Tyr2 – Pro4` | 8 × 4 | 32 | 1 – 32 |
| 2 | `Tyr2 – Trp9` | 8 × 11 | 88 | 33 – 120 |
| 3 | `Pro4 – Trp9` | 4 × 11 | 44 | 121 – 164 |
| **TOTAL :** | | | **164** | |

The last column (`DISANG` order) is fixed from this moment on. The `INDXF` file addresses restraints by their position in the `DISANG` and by nothing else, so deciding the groups has also decided the order the `DISANG` must be written in.

For this tutorial, we will not cover the full procedure for writing an AMBER-compatible NMR restraint file (`DISANG`). Instead, we use `CPPTRAJ` to generate [hy.disang](1_system_setup/hy.disang) from the commands in [rst.cpptraj.in](1_system_setup/rst.cpptraj.in).

The input begins by loading the **1UAO** topology. Each `rst` command then defines one atom-pair distance restraint and appends it to `hy.disang`. The restraints are written in three consecutive blocks corresponding to the Tyr2–Pro4, Tyr2–Trp9, and Pro4–Trp9 contact groups.

```text
parm 1uao.prmtop

# ---- group 1: Tyr2–Pro4; 8 × 4 = 32 restraints (DISANG 1–32) ----
rst :2@CA :4@CD r1 0.000 r2 0.000 r3 5.000 r4 7.000 rk2 0.29876 rk3 0.29876 out hy.disang   # 1
rst :2@CA :4@CG r1 0.000 r2 0.000 r3 5.000 r4 7.000 rk2 0.29876 rk3 0.29876 out hy.disang   # 2
rst :2@CA :4@CB r1 0.000 r2 0.000 r3 5.000 r4 7.000 rk2 0.29876 rk3 0.29876 out hy.disang   # 3
...
...
rst :4@CA :9@CE3 r1 0.000 r2 0.000 r3 5.000 r4 7.000 rk2 0.29876 rk3 0.29876 out hy.disang   # 163
rst :4@CA :9@CD2 r1 0.000 r2 0.000 r3 5.000 r4 7.000 rk2 0.29876 rk3 0.29876 out hy.disang   # 164

# 164 restraints in total: GROUPS 32 88 44
quit
```

Here, `r2`–`r3` defines the zero-penalty distance range. Because both are set to 0 and 5 Å, respectively, an atom pair is unpenalized up to 5 Å. Between `r3 = 5 Å` and `r4 = 7 Å`, the penalty increases quadratically; beyond 7 Å, it increases linearly. The values of `rk2` and `rk3` set the restraint force constants.

The complete file contains 164 candidate distance restraints divided into groups of 32, 88, and 44. These group sizes must match the `GROUPS 32 88 44` declaration in the MELD index file.

Run `CPPTRAJ` with:

```bash
cpptraj -i rst.cpptraj.in
```

A detailed description of the `rst` command and its parameters is available from the [AMBER-Hub restraint documentation](https://amberhub.chpc.utah.edu/rst/).

#### Writing the `INDXF`

The index file is where MELD actually lives. It holds the group and collection hierarchy, the force constants, and the replica ladder — and it identifies restraints purely by their position in the `DISANG`.

> **`DISANG` is a numbered list of geometry. The `INDXF` is a hierarchy over those numbers. The only thing that links them is position.**

```
# hy.indxf -- the hierarchy over hy.disang's 164 restraints.
# Pair this with hy.disang.  Nothing but POSITION links the two: "restraint 33" below means the 33rd &rst record in that file.
#

UNITS  meld           # K below is kJ/mol/nm^2; pmemd converts per restraint type
INFO   on             # every replica writes meld.info.<rank>

# ---- the ladder ----------------------------------------------------------
TSCALE  0.0 1.0  300.0 450.0  "geometric"

SCALER  prot   0.4 1.0 4.0        "nonlinear"
RAMP    warmup 1 200 1e-3 1 4.0   "nonlinear_ramp"

# ---- the hierarchy -------------------------------------------------------
# hydrophobic core: 3 residue pairs, every apolar atom pair, keep one
# 3 groups, 164 restraints
COLL HY  ACTIVE 2  KEEP 1  SCALER prot  RAMP warmup  K distance 250  GROUPS 32 88 44
```

Reading it from the top:

**`UNITS meld`** declares that the `K` values below are in MELD's own units — kJ/mol/nm² for a distance, kJ/mol/deg² for a torsion — and pmemd converts them on the way in. It **must appear before the collections**, because each restraint's force constant is converted at the moment its group claims it; a `UNITS` record at the bottom of the file converts nothing. Without it the default is `amber`, and `K distance 250` would be read as 250 kcal/mol/Å², roughly 800 times too strong.

**`INFO on`** makes every replica write a `meld.info.<rank>` file recording what it selected — the only direct evidence of what the collection is doing, and worth having on the first run.

**`TSCALE 0.0 1.0 300.0 450.0 "geometric"`** is the temperature ladder: replica *n* of *N* sits at alpha = (n−1)/(N−1), and its bath temperature is read off this curve, geometrically spaced from 300 K at alpha 0 to 450 K at alpha 1. **This overrides `temp0` in the mdin.** Several `TSCALE` records make a piecewise ladder; they must tile their range with no gap and no overlap, and outside the range they cover the ladder clamps.

**`SCALER prot 0.4 1.0 4.0 "nonlinear"`** is the force-constant ladder: full strength from alpha 0 up to 0.4, then decaying nonlinearly to a thousandth of it by alpha 1.0. The hot replicas are effectively free of the hydrophobic restraints, which is what lets them explore; the cold ones feel them fully. The name `prot` is arbitrary — the collection refers to it by name.

**`RAMP warmup 1 200 1e-3 1 4.0 "nonlinear_ramp"`** fades the restraints in over the first 200 **exchange** steps, starting at a thousandth of their force constant, so an extended starting chain is not yanked into a ball on step one.

**`COLL HY …`** is the hierarchy, and everything before `GROUPS` is a setting that applies to the groups after it:

| Field | Meaning |
| :-- | :-- |
| `HY` | the collection's name; required on a one-line `COLL` record, and it is what `meld.info` reports |
| `ACTIVE 2` | enforce the 2 lowest-energy groups of the 3 |
| `KEEP 1` | each group keeps its 1 lowest-energy restraint |
| `SCALER prot` | scale this collection's force constants by `prot`; must be defined **above** this line |
| `RAMP warmup` | fade them in with `warmup`; likewise defined above |
| `K distance 250` | the base force constant, per restraint type. **There is no default — an unset `K` is fatal.** |
| `GROUPS 32 88 44` | three groups, of 32, 88 and 44 restraints |


So `32 88 44` builds {1…32}, {33…120}, {121…164}, in the order the generator wrote them. Nothing states residue names or atom numbers anywhere in this file; it is 164 restraints sliced into three pieces, and the slicing is only correct because the `DISANG` was written in the planned order. Insert one `&rst` record in the middle of that file and every group afterwards silently points at different physics.

Two shorthands are worth knowing. Repeated sizes can be written `size x count` — `GROUPS 5 x 3` is three groups of five restraints each — and every restraint the `DISANG` must be claimed exactly once: claimed twice is fatal, unclaimed is fatal, and there is no default collection.

The same collection can also be written as a block, which is the spelling to prefer once a record gets long:

```
COLLECTION HY
  ACTIVE 2
  KEEP 1
  SCALER prot
  RAMP warmup
  K distance 250
  GROUPS 32 88 44
END
```

> **Note:** The two are read identically, with one difference that matters: a one-line `COLL` is read settings-first, so the order of fields on it is free, while a block is read **top to bottom**, so every setting must appear *above* the `GROUPS` line it is meant to reach. Move `K distance 250` below `GROUPS` in the block above and the file is no longer valid.

If a group is not a contiguous run of restraints, `GROUPS` cannot express it; use `MEMBERS`, which names indices outright and where **one record is exactly one group**.

The quickest way to see what it does is to write a collection you already know. Here is `HY` — the same three groups, the same hierarchy, the same physics — with `GROUPS 32 88 44` expanded into the indices it stands for:

```
COLLECTION HY
  ACTIVE 2
  KEEP 1
  SCALER prot
  RAMP warmup
  K distance 250
  MEMBERS 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32    # group 1: Tyr2 - Pro4
  MEMBERS 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120    # group 2: Tyr2 - Trp9
  MEMBERS 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164    # group 3: Pro4 - Trp9
END
```

Three `MEMBERS` records, three groups, of 32, 88 and 44 restraints. Read against the `GROUPS` spelling above, the correspondence is exact: `GROUPS` counts off the next unclaimed restraints, `MEMBERS` names them, and for a hierarchy whose groups are contiguous those are two ways of writing one thing. pmemd reads either file into the same collection.

Our groups are contiguous by construction, which is the whole reason the generator writes the `DISANG` in hierarchy order.

### 3. Launch MELD simulations

Everything MELD reads now exists: `hy.disang` holds the 164 restraints, `hy.indxf` holds the hierarchy over them, and minimized coordinates in `min.rst`. What is left is ordinary Amber replica exchange — one input file per replica, a groupfile listing them, and one `pmemd` command.

#### The sample input

Write one `mdin` and let it stand for the whole ladder. MELD is switched on with `meld=1` alongside `nmropt=1`; the index file is named by `indxf=` inside `&cntrl`, and the `DISANG` goes where Amber always puts it, after the `&wt TYPE='END'` block:

```
chignolin -- MELD with a hydrophobic-core collection, 50 ns per replica
 &cntrl
   irest=0, ntx=1,
   nstlim=5000, dt=0.004, numexchg=2500,
   ntt=3, gamma_ln=1.0,
   temp0=300.0, ig=-1,
   ntc=2, ntf=2,
   ntb=0, igb=8, cut=999.0, rgbmax=999.0,
   ntpr=5000, ntwx=5000, ntwr=5000,
   nmropt=1, meld=1, indxf='hy.indxf',
 /
 &wt TYPE='END'
 /
DISANG=hy.disang
```

The three lines that set the length of the run work together: `nstlim=5000` steps of `dt=0.004` ps is 20 ps of dynamics between exchange attempts, and `numexchg=2500` of those attempts is **50 ns per replica**. The 4 fs timestep is what the HMR topology from §1 bought, and `ntc=2, ntf=2` is the `SHAKE` setting it requires. `igb=8` matches the `mbondi3` radii the topology was built with — the same model as the minimization, which is what you want.

Two entries are easy to misread. `temp0` is written because Amber requires it and then **ignored**: the `TSCALE` record in `hy.indxf` sets every replica's bath temperature, and `mdout` reports the value actually used. `ig=-1` asks for a random seed, and the next subsection explains why each replica needs its own.

> **`DISANG=` must be the line immediately after the `&wt TYPE='END'` namelist.** Amber stops reading redirections at the first line that is not one, so a blank line or a comment in that gap silently drops every restraint in the run — the job starts, finishes, and has simply done unrestrained MD.

 **A replica's position in the groupfile is its position on the ladder**: line *n* of *N* gets alpha = (*n*−1)/(*N*−1), and its temperature and force-constant scaling are read off the `TSCALE` and `SCALER` curves at that alpha. Line 1 is the cold, fully-restrained replica; line *N* is the hot, nearly-free one.

Six files can be copied by hand. Past that it is worth generating them, which is what [gen_md_inputs.py](../2_Protein_Folding/2_meld_remd/gen_md_inputs.py) from the folding tutorial does — it copies your sample once per replica, changes those two fields and nothing else, and writes the groupfile to match:

```
../2_Protein_Folding/2_meld_remd/gen_md_inputs.py \
    -i meld.mdin --indxf hy.indxf -n 6 \
    -c ../1_system_setup/min.rst \
    -p ../1_system_setup/1uao_HMR.prmtop
```

This writes `meld.mdin.001` … `meld.mdin.006` and [meld.groupfile](2_meld_remd/meld.groupfile) given the template [meld.mdin](2_meld_remd/meld.mdin).

> **The replica count must be even.** pmemd refuses an odd count.

#### Run it

`-ng` is the number of groups and must equal the number of lines in the groupfile:

```
srun --mpi=pmix_v5 $AMBERHOME/bin/pmemd.cuda.MPI -ng 6 -groupfile meld.groupfile
```

## Where to go next

This tutorial built one collection from one kind of information. The [protein-folding tutorial](../2_Protein_Folding/README.md) builds all three CPI collections for a 56-residue domain and runs the replica ladder.
