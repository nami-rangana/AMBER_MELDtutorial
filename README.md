# MELD in AMBER — Tutorials

**Modeling Employing Limited Data (MELD)** is a physics-based, Bayesian framework that combines atomistic molecular dynamics with structural information that may be sparse, ambiguous, or only partly reliable <sup>[[1](https://www.pnas.org/doi/full/10.1073/pnas.1506788112)]</sup>. To explore alternative interpretations of these data while avoiding kinetic traps, MELD uses combined Hamiltonian and temperature replica-exchange molecular dynamics. This repository provides a practical series of AMBER tutorials that follows the complete workflow: starting from a FASTA sequence, constructing the restraints, running the GPU-accelerated replica ladder, and analyzing the resulting ensemble to determine whether the protein adopted its expected fold.

| Tutorial | System | What you build |
| :-- | :-- | :-- |
| [1_Basic_tutorials](1_Basic_tutorials/README.md) | Chignolin, 10 residues (1UAO) | One collection, written by hand, restraint by restraint |
| [2_Protein_Folding](2_Protein_Folding/README.md) | Protein G B1 domain, 56 residues (3GB1) | The three CPI collections; a 30-replica, 200 ns/replica fold |

*New to MELD? Read the theory below, then work through the tutorials in order.*

---
## 1. The problem MELD addresses

MELD addresses two closely related challenges in protein structure prediction: efficiently exploring conformational space and making useful structural inferences from imperfect data.

**Conformational sampling.** In principle, atomistic molecular dynamics can describe protein folding using a physical energy model. Unlike methods that produce only a single structure, molecular dynamics can also provide information about conformational populations, motions, stability, and folding mechanisms. In practice, however, the number of possible conformations is enormous, and folding often occurs on timescales that are difficult to reach with conventional simulations. Sampling and forcefield limitations also become more important as the size and complexity of the system increase.

**Semireliable structural information.** Experimental measurements and bioinformatic predictions can greatly reduce the search space, but the available information is not always complete or fully reliable. MELD distinguishes three common forms of limited data:

| Class | Description | Example |
| :-- | :-- | :-- |
| **Sparse** | Only a small number of specific, reliable restraints are available | A few solid-state NMR methyl–methyl contacts |
| **Ambiguous** | The general structural relationship is reliable, but its exact atomic interpretation is not known | EPR/DEER measurements or nonspecific hydrophobic contacts |
| **Uncertain** | Some entries are correct and others are incorrect, but their identities are not known in advance | Evolutionary-coupling contact predictions with limited accuracy |

Conventional restrained molecular dynamics typically assumes that every supplied restraint should be satisfied. If a substantial fraction of the input is incorrect, enforcing all restraints simultaneously can distort the structure or trap the simulation in a nonnative state.

MELD takes a different approach. It incorporates sparse, ambiguous, and uncertain information as semireliable evidence within a physics-based Bayesian framework <sup>[[1](https://www.pnas.org/doi/full/10.1073/pnas.1506788112)]</sup>. Rather than requiring every restraint to be satisfied at once, MELD activates only a defined subset of the available information. The restraints narrow the conformational search and promote a funnel toward data-compatible structures, while the atomistic force field helps distinguish physically realistic conformations from alternative structures that may satisfy the same limited data. In this way, the same information that makes sampling more efficient can be used without assuming that every input restraint is correct.

---

## 2. The Bayesian framework

MELD applies Bayesian inference to molecular simulation. Let \(\mathbf{x}\) represent the atomic coordinates of a protein and \(D\) the available structural data. Bayes’ rule gives

```math
p(\mathbf{x}\mid D)
=
\frac{p(D\mid \mathbf{x})\,p(\mathbf{x})}{p(D)}
\propto
p(D\mid \mathbf{x})\,p(\mathbf{x})
```

where $\(p(\mathbf{x}\mid D)\)$ is the **posterior distribution**, $\(p(\mathbf{x})\)$ is the **prior**, and $\(p(D\mid\mathbf{x})\)$ is the **likelihood**. MELD samples from this posterior distribution rather than searching only for a single structure with the best score.

### The prior: the physical model

The prior describes which conformations are physically plausible before the external structural information is considered. In MELD, this distribution is defined by the Amber potential energy:

$$
p(\mathbf{x})
\propto
\exp\!\left[-\beta E_{\mathrm{Amber}}(\mathbf{x})\right],
\qquad
\beta=\frac{1}{k_{\mathrm B}T}.
$$

The force field provides the physical description of the protein, including bonded interactions, steric packing, electrostatics, and solvation. It therefore distinguishes physically realistic conformations from structures that may satisfy the data geometrically but contain unfavorable molecular interactions.

### The likelihood: agreement with the data

The likelihood describes how compatible a conformation is with the available structural information. Individual observations are represented as nonnegative restraint energies involving geometrical quantities such as distances, angles, or torsions [1]. For an observation \(D_i\), the corresponding likelihood can be written as

$$
p(D_i\mid\mathbf{x})
\propto
\exp\!\left[-\beta E_i^{\mathrm{rest}}(\mathbf{x})\right].
$$

The restraint energy \(E_i^{\mathrm{rest}}\) is small when the conformation agrees with the observation and increases as the disagreement becomes larger.

Combining the physical prior and data likelihood gives the effective potential sampled during the simulation:

$$
E_{\mathrm{total}}(\mathbf{x})
=
E_{\mathrm{Amber}}(\mathbf{x})
-
\beta^{-1}\ln p(D\mid\mathbf{x}).
$$

The force field and the data therefore play complementary roles. The restraints guide sampling toward regions that agree with the available information, while the physical model determines which conformations within those regions are energetically reasonable [1].

If every observation were assumed to be correct and independent, the likelihood term would reduce to a sum over all \(M\) restraint energies:

$$
E_{\mathrm{total}}(\mathbf{x})
=
E_{\mathrm{Amber}}(\mathbf{x})
+
\sum_{i=1}^{M}E_i^{\mathrm{rest}}(\mathbf{x}).
$$

This is the usual restrained-molecular-dynamics model, in which every restraint is enforced simultaneously. MELD instead allows the data to be sparse, ambiguous, or uncertain and does not require every supplied restraint to be correct [1].

### 2.1 Selecting a mutually compatible subset

Suppose a dataset contains \(M=100\) predicted contacts and previous experience suggests that approximately \(N=65\) are reliable. The appropriate assumption is not that all 100 contacts are correct, but that **about 65 of them are correct, although their identities are unknown**.

MELD represents this uncertainty by evaluating and ranking the restraint energies for the current conformation:

$$
E_{(1)}^{\mathrm{rest}}(\mathbf{x})
\leq
E_{(2)}^{\mathrm{rest}}(\mathbf{x})
\leq
\cdots
\leq
E_{(M)}^{\mathrm{rest}}(\mathbf{x}).
$$

Only the \(N\) lowest-energy restraints are then included in the effective potential:

$$
E_{\mathrm{total}}(\mathbf{x})
=
E_{\mathrm{Amber}}(\mathbf{x})
+
\sum_{i=1}^{N}E_{(i)}^{\mathrm{rest}}(\mathbf{x}),
\qquad N\leq M.
$$

The remaining \(M-N\) restraints contribute neither energy nor force during that evaluation. This selection is dynamic: as the protein changes conformation, the restraint energies are recalculated and a different subset may become active.

The number of restraints treated as reliable must be chosen by the user. Setting this number too high may force the simulation to account for incorrect information, whereas setting it too low may discard useful structural information.

This treatment has three important consequences.

#### Uncertain information can be tolerated

Incorrect restraints are not identified permanently at the beginning of the simulation. Instead, MELD favors conformations that satisfy a large, mutually compatible subset of the data while remaining physically reasonable. Restraints that are inconsistent with low-free-energy structures are less likely to remain active. In this sense, MELD samples both the protein structure and compatible interpretations of the available data [1].

#### Ambiguity can be represented explicitly

Consider a contact indicating that two side chains are close without specifying which atoms form the interaction. The possible atom-pair distances can be represented as alternatives, and MELD can select the interpretation that best matches the current conformation. The ambiguity therefore does not need to be resolved before the simulation begins.

#### Restraints guide rather than prescribe the structure

MELD restraints commonly contain a broad, flat region in which the energetic penalty is zero or small [1]. Once a restraint is satisfied, the force field determines the detailed geometry. The data narrow the conformational search without replacing the underlying physical model.

Because different conformations can activate different subsets of restraints, the resulting energy landscape may contain several competing funnels. MELD samples among these alternatives and favors the structural interpretations associated with the lowest free energy.

### 2.2 Restraints, groups, and collections

Real structural information often contains uncertainty at more than one level. MELD therefore organizes the input into a hierarchy of **restraints**, **groups**, and **collections**.

| Level | Interpretation | Selection parameter |
| :-- | :-- | :-- |
| **Restraint** | One geometrical condition, such as a distance or torsion | — |
| **Group** | One structural claim containing one or more geometrical interpretations | `KEEP` |
| **Collection** | A set of related structural claims with a shared expected reliability | `ACTIVE` |

The hierarchy is evaluated in two stages.

#### Selection within a group

Suppose group \(g\) contains \(m_g\) restraints with ordered energies

$$
E_{g,(1)}
\leq
E_{g,(2)}
\leq
\cdots
\leq
E_{g,(m_g)}.
$$

If the group has `KEEP` \(k_g\), its energy is the sum of its \(k_g\) lowest-energy restraints:

$$
E_g^{\mathrm{group}}(\mathbf{x})
=
\sum_{j=1}^{k_g}E_{g,(j)}^{\mathrm{rest}}(\mathbf{x}).
$$

The `KEEP` setting therefore describes how much of the information within one structural claim must be satisfied:

| Setting | Meaning |
| :-- | :-- |
| `KEEP 1` | Select the best-satisfied restraint in the group |
| `KEEP n` | Select the \(n\) best-satisfied restraints |
| `KEEP all` | Require every restraint in the group to contribute |

For example, an ambiguous contact between two side chains may be represented by many possible atom-pair distances placed in a single group with `KEEP 1`. The structural claim is that the two residues are in contact—not that every atom in one side chain must contact every atom in the other.

By contrast, several distances may be required together to define a particular β-strand register. These distances are complementary parts of the same geometry rather than alternatives. They would normally be placed in one group with `KEEP all`, or with another value that reflects how many distances must be satisfied for the pairing to be meaningful.

#### Selection within a collection

A collection contains groups representing comparable structural claims. After the energy of each group has been calculated, MELD ranks the groups within the collection:

$$
E_{(1)}^{\mathrm{group}}
\leq
E_{(2)}^{\mathrm{group}}
\leq
\cdots
\leq
E_{(G)}^{\mathrm{group}}.
$$

If the collection specifies `ACTIVE` \(A\), only the \(A\) lowest-energy groups contribute:

$$
E^{\mathrm{collection}}(\mathbf{x})
=
\sum_{g=1}^{A}E_{(g)}^{\mathrm{group}}(\mathbf{x}),
\qquad A\leq G.
$$

The remaining groups are inactive during that force evaluation. Thus, `ACTIVE` expresses how many group-level claims are expected to be compatible with the structure, while `KEEP` expresses how much of each selected claim must be satisfied.

The complete hierarchy can be summarized as

$$
\text{collection}
\;\xrightarrow{\;\texttt{ACTIVE}\;}
\text{selected groups}
\;\xrightarrow{\;\texttt{KEEP}\;}
\text{selected restraints}.
$$

For example, a collection containing 121 groups might select 21:

$$
121\ \text{groups}
\;\xrightarrow{\;\texttt{ACTIVE}\ 21\;}
21\ \text{active groups}.
$$

Within each selected group, one restraint might be chosen from 88 alternative distances:

$$
88\ \text{alternative distances}
\;\xrightarrow{\;\texttt{KEEP}\ 1\;}
1\ \text{active distance}.
$$

The scientific meaning of the hierarchy can therefore be summarized as follows:

- A **restraint** defines one geometrical condition.
- A **group** defines one structural claim and its internal ambiguity.
- A **collection** defines a set of related claims and how many are expected to be reliable.

This distinction is important. Writing 88 alternative atom-pair distances as 88 independent groups would treat them as separate structural claims. Placing them in one group with `KEEP 1` instead states that they are alternative atomic interpretations of a single residue-level contact.

### 2.3 Always-active restraints

MELD can also define restraints that bypass both levels of selection. These restraints remain active throughout the simulation and are not ranked against alternative restraints or groups.

Always-active restraints are appropriate only when the underlying information is considered effectively certain, such as a known covalent linkage or another chemically required interaction. Sparse, ambiguous, or uncertain predictions should normally remain within selectable groups and collections so that MELD can evaluate alternative interpretations.

An incorrect always-active restraint cannot be rejected by the selection procedure and may therefore bias the entire sampled ensemble.

---

## 3. Sampling the posterior with H,T-REMD

The posterior distribution explored by MELD contains many local minima separated by energetic barriers. A conventional molecular-dynamics trajectory can become trapped in one of these regions and may not sample alternative conformations within an accessible simulation time.

To improve sampling, MELD uses **Hamiltonian and temperature replica-exchange molecular dynamics (H,T-REMD)**. Multiple copies of the system are simulated in parallel along a replica ladder. The replicas differ in both temperature and restraint strength: replicas near the bottom of the ladder remain close to room temperature and experience strong restraints, whereas those near the top use higher temperatures and weaker restraints. This arrangement allows the upper replicas to explore conformational space broadly while the lower replicas concentrate sampling in low-free-energy regions that agree with the supplied information [1].

### 3.1 The replica ladder

A dimensionless parameter, \(\alpha\), identifies each position along the ladder. For replica \(n\) in a ladder containing \(N\) replicas,

$$
\alpha_n = \frac{n-1}{N-1},
\qquad n=1,\ldots,N.
$$

The two ends of the ladder are therefore

$$
\alpha_1=0
\qquad\text{and}\qquad
\alpha_N=1.
$$

Temperature and restraint strength are defined as functions of \(\alpha\):

$$
T_n=T(\alpha_n),
$$

$$
k_n=s(\alpha_n)\,k_0,
$$

where \(k_0\) is the base force constant and \(s(\alpha)\) is a scaling function applied along the ladder.

| Ladder position | Temperature | Restraint strength | Main role |
| :-- | :-- | :-- | :-- |
| \(\alpha=0\), the lowest replica | Typically near 300 K | Full or nearly full strength | Samples physically realistic, data-compatible conformations |
| \(\alpha=1\), the highest replica | Higher temperature | Weak or nearly absent | Crosses barriers and explores alternative conformations |

The exact temperature range and scaling profile are modeling choices rather than universal MELD constants. In the original applications, MELD used between 24 and 48 replicas, with exchange attempts every 20 or 50 ps. Some calculations increased the temperature from 300 to 450 K while reducing nonlocal distance-restraint strength toward the upper part of the ladder [1]. A particular tutorial may use a different upper temperature, replica count, or exchange interval.

### 3.2 Exchange between replicas

At regular intervals, neighboring replicas attempt to exchange configurations. The exchange criterion accounts for both the temperatures and Hamiltonians of the two replicas and is constructed to preserve the desired equilibrium distributions.

A conformation trapped near the bottom of the ladder can therefore move upward through successful exchanges. At higher temperature and weaker restraint strength, it can cross barriers, partially unfold, or adopt a different topology. It may later return to the lower replicas in a new region of conformational space.

The ladder can be viewed as dividing the sampling problem into two complementary tasks:

- the **upper replicas** promote exploration by weakening restraints and increasing thermal motion;
- the **lower replicas** evaluate the resulting conformations under physically relevant temperatures and stronger agreement with the data.

This exchange process does not guarantee that every relevant state will be sampled. Convergence still depends on factors such as simulation length, exchange efficiency, system size, force-field accuracy, and the design of the restraint hierarchy. It does, however, provide a systematic way for conformations to escape kinetic traps that would be difficult to overcome in a single conventional trajectory.

### 3.3 Temperature scalers, restraint scalers, and ramps

Three controls are commonly used to define how the simulation changes across the ladder and over time:

| Control | Purpose |
| :-- | :-- |
| `TSCALE` | Defines the temperature as a function of \(\alpha\) |
| `SCALER` | Defines how restraint strength changes with \(\alpha\) |
| `RAMP` | Defines how restraint strength changes with simulation time |

A geometric temperature schedule provides a smooth increase from the lower to the upper replicas. A restraint scaler can preserve full restraint strength near the bottom of the ladder and progressively weaken it toward the top.

Local and nonlocal information do not always need the same scaling behavior. Secondary-structure restraints describe local backbone geometry and may remain active throughout the ladder. By contrast, nonlocal distance restraints can strongly restrict the global fold and are often weakened in the upper replicas. This allows the protein to escape an incorrect topology before returning to the strongly restrained region. The original MELD study used this distinction in several applications: secondary-structure restraints remained at full strength, while imposed distance restraints weakened toward the upper replicas [1].

A time-dependent `RAMP` serves a different purpose. Rather than varying restraint strength between replicas, it introduces the restraints gradually during the initial part of the simulation. This is particularly helpful when all replicas begin from an extended chain. Applying strong nonlocal restraints immediately could produce abrupt collapse, poor local geometry, or large initial forces. A gradual ramp allows the system to relax while the data-derived potential is introduced.

### 3.4 Interpreting the sampled ensemble

When replica exchange is implemented correctly and sampling is converged, each ladder position samples the equilibrium distribution associated with its own temperature and Hamiltonian. The lowest replica therefore samples the restrained distribution defined by the physical model and the supplied data—the MELD posterior discussed in the previous section.

This distinction is important:

> The lowest-replica ensemble represents the force field **conditioned on the supplied information**. It is not the unbiased ensemble of the force field alone.

MELD is therefore intended to produce an ensemble rather than simply generate a single structure. Relative populations can provide information about the stability of competing conformational states, provided that sampling is sufficiently converged. In practice, low-free-energy regions are commonly identified by clustering the sampled conformations. The original MELD work emphasized that individual potential-energy values should not be used directly as free-energy scores because they do not include the entropic contribution [1].

A folded structure appearing once is consequently weaker evidence than a folded basin that is repeatedly visited and well populated. Analysis should focus on the distribution of conformations, transitions between states, and the stability of structural clusters rather than only on the single frame with the lowest energy or RMSD.

---

## 4. Sources of structural restraints

MELD can incorporate many forms of structural information, provided that the information can be translated into a suitable restraint potential and organized according to its uncertainty. The source of the information determines the geometry of the restraints, while its reliability and ambiguity determine how those restraints should be arranged into groups and collections.

### 4.1 Experimental and predicted information

Potential sources of structural information include NMR-derived contacts and torsions, EPR/DEER distance distributions, chemical cross-links, mutagenesis-derived contacts, and predicted residue–residue interactions. A typical MELD calculation may combine the protein sequence, a predicted secondary structure, and externally supplied residue-contact or distance information [1].

The conversion from an observation to a MELD restraint requires a scientific interpretation. For example:

- an assigned NOE may be represented by an interatomic distance;
- chemical shifts may be converted into backbone torsion ranges;
- a cross-link may define an upper-bound distance between candidate atoms;
- a residue-level contact prediction may correspond to several possible atom-pair distances;
- an uncertain contact list may contain separate groups, only a chosen fraction of which are active.

The restraint hierarchy should preserve the meaning of the original observation. Alternative atomic interpretations belong within a group, where `KEEP` controls how many alternatives are selected. Comparable observations with shared reliability belong within a collection, where `ACTIVE` controls how many groups contribute at a given force evaluation.

The [custom-collection tutorial](2_Protein_Folding/additional_tutorials/README.md) demonstrates this process using a predicted contact map and a chemical-shift-derived torsion table.

### 4.2 Coarse Physical Insights derived from sequence

MELD can also guide protein folding when no protein-specific experimental restraints are available. In this setting, it uses **Coarse Physical Insights (CPI)**: broad expectations about protein structure that are informative but not precise enough to determine the fold on their own.

These insights are useful because they reduce the conformational search without specifying a complete native structure. Their statistical and ambiguous nature is handled through MELD’s group-and-collection selection scheme.

The CPI restraints used in these tutorials are organized into three collections:

| Collection | Physical idea | Source of uncertainty |
| :-- | :-- | :-- |
| **`SS` — secondary structure** | Local regions tend to adopt the secondary structure predicted from the sequence | Some predicted residues or windows may be incorrect |
| **`SP` — strand pairing** | Predicted β-strands are likely to pair with other strands | The partner strand, orientation, and register may be unknown |
| **`HY` — hydrophobic contacts** | Hydrophobic residues tend to form a compact core in globular proteins | Only a small subset of all possible hydrophobic pairs should contact |

#### Secondary-structure collection (`SS`)

Secondary-structure predictions are translated into local restraints, commonly acting on overlapping sequence windows. These compound restraints may include both backbone torsions and local distances. The original MELD study used overlapping five-residue fragments and activated 75% of the resulting secondary-structure groups, reflecting the expectation that some predictions would be incorrect [1].

The exact active fraction is a modeling parameter and may differ between implementations or tutorials. Its purpose is not to weaken every secondary-structure restraint uniformly, but to allow MELD to leave out local assignments that are incompatible with the sampled conformation.

#### Strand-pairing collection (`SP`)

The strand-pairing collection describes possible interactions between predicted β-strands. The sequence may indicate which regions have β-strand character, but it does not uniquely determine which strands pair, whether they align in parallel or antiparallel orientations, or which residue register is adopted.

MELD can represent these possibilities as competing groups. Each group describes one proposed strand-pairing pattern, and only a limited number are activated. This allows the simulation to explore alternative β-sheet topologies without enforcing every possible pairing simultaneously.

#### Hydrophobic-core collection (`HY`)

The hydrophobic collection represents the general tendency of nonpolar side chains to become buried and form a compact core. It does not assume that every hydrophobic residue contacts every other hydrophobic residue. Such a requirement would overcompact the chain and would not reflect the packing geometry of real proteins.

Instead, possible hydrophobic contacts are generated as an intentionally overcomplete set, and MELD activates only a small subset. The force field then determines which compatible contacts produce physically realistic side-chain packing.

Together, the three collections provide complementary information:

- `SS` promotes plausible local backbone geometry;
- `SP` organizes β-strands into possible sheet topologies;
- `HY` favors the nonlocal collapse and packing expected in a globular fold.

The protein sequence and its secondary-structure prediction are therefore used to construct broad structural hypotheses rather than a predetermined native structure. MELD combines these hypotheses with the atomistic force field and evaluates competing interpretations through equilibrium sampling.

The complete CPI folding workflow is demonstrated in [2_Protein_Folding](2_Protein_Folding/README.md).

---

## 5. Glossary

| Term | Meaning |
| :-- | :-- |
| **alpha** | Ladder coordinate, 0 (cold, fully restrained) to 1 (hot, free). Sets temperature *and* restraint strength. |
| **Restraint** | One flat-bottom distance or torsion. |
| **Group** | One structural claim; `KEEP n` of its restraints are summed. |
| **Collection** | Comparable claims ranked against each other; `ACTIVE N` groups are enforced. |
| **`KEEP`** | How many restraints in a group are enforced. Resolves *internal* ambiguity. |
| **`ACTIVE`** | How many groups in a collection are enforced. Encodes *how much of the data you believe*. |
| **`ALWAYS`** | Restraints exempt from all selection. Chemistry only. |
| **`TSCALE`** | Temperature-vs-alpha curve. Overrides `temp0`. |
| **`SCALER`** | Force-constant-vs-alpha curve. |
| **`RAMP`** | Force-constant-vs-time curve; fades restraints in at the start. |
| **CPI** | Coarse Physical Insights — generic knowledge (SS, SP, HY) requiring no experimental data. |
| **`DISANG`** | AMBER restraint geometry. A numbered list. |
| **`INDXF`** | MELD index file. The hierarchy over those numbers, plus the ladder. |

---

## References

1. MacCallum, J. L.; Perez, A.; Dill, K. A. Determining Protein Structures by
   Combining Semireliable Data with Atomistic Physical Models by Bayesian
   Inference. *Proc. Natl. Acad. Sci. U.S.A.* **2015**, *112* (22), 6985–6990.
   [doi:10.1073/pnas.1506788112](https://doi.org/10.1073/pnas.1506788112) ·
   [PMC4460504](https://pmc.ncbi.nlm.nih.gov/articles/PMC4460504/)
   — *the MELD theory paper; §2 of this README follows it.*
2. Perez, A.; MacCallum, J. L.; Dill, K. A. Accelerating Molecular Simulations of
   Proteins Using Bayesian Inference on Weak Information. *Proc. Natl. Acad. Sci.
   U.S.A.* **2015**, *112* (38), 11846–11851.
   [doi:10.1073/pnas.1515561112](https://doi.org/10.1073/pnas.1515561112) ·
   [PMC4586851](https://pmc.ncbi.nlm.nih.gov/articles/PMC4586851/)
   — *the CPI paper; §4.2.*
3. Perez, A.; Morrone, J. A.; Dill, K. A. Accelerating Physical Simulations of
   Proteins by Leveraging External Knowledge. *WIREs Comput. Mol. Sci.* **2017**,
   *7*, e1309. [doi:10.1002/wcms.1309](https://doi.org/10.1002/wcms.1309) ·
   [PMC5612641](https://pmc.ncbi.nlm.nih.gov/articles/PMC5612641/)
   — *review; a good second read.*
4. Hopkins, C. W.; Le Grand, S.; Walker, R. C.; Roitberg, A. E. Long-Time-Step
   Molecular Dynamics through Hydrogen Mass Repartitioning. *J. Chem. Theory
   Comput.* **2015**, *11* (4), 1864–1874.
   [doi:10.1021/ct5010406](https://doi.org/10.1021/ct5010406)
   — *the 4 fs timestep.*

