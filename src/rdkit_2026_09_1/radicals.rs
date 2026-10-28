//! Immutable parsed-graph admission and versioned radical assignment.

use alloc::vec::Vec;

use thiserror::Error;

use super::element_data::ELEMENTS;
use crate::{
    NitrogenOxideCleanupError, Smiles, SmilesErrorWithSpan, TetrahedralNeighbor,
    atom::{Atom, AtomSyntax, bracketed::chirality::Chirality},
    bond::Bond,
};

/// Whether the versioned assignment rule assigns or preserves the count.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RadicalDecision {
    /// The original bracket atom has NoImplicit set and receives a count.
    Assigned,
    /// The original organic-subset atom retains its parser count of zero.
    Preserved,
}

/// A native valence cache whose initialized range is 0..=127.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NativeValenceCache {
    /// Bond contributions plus nonvertex explicit hydrogens.
    Explicit,
    /// Implicit valence, equal to implicit H in this admitted domain.
    Implicit,
}

/// A conservatively excluded cleanup family.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UnsupportedPreparation {
    /// An N-triple-N edge can participate in the separate nitrogen cleanup.
    NitrogenTripleNitrogen,
    /// A neutral P-double-O environment can require phosphorus cleanup.
    PhosphorusOxide,
    /// A neutral Cl/Br/I-double-O environment can require halogen cleanup.
    HalogenOxide,
}

/// An input for which no versioned parsed capability was constructed.
///
/// Unsupported states are not claims of native chemical invalidity. Errors
/// never substitute a zero count or modify an existing graph.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum RadicalPreparationError {
    /// The ordinary concrete SMILES parser rejected the input.
    #[error("SMILES parse error: {0}")]
    Parse(SmilesErrorWithSpan),
    /// Aromatic atoms or bonds need a separately qualified preparation path.
    #[error("aromatic state at atom {atom} is outside radical preparation")]
    UnsupportedAromaticState {
        /// Current/source atom index.
        atom: usize,
    },
    /// The acyclic constructor encountered a topological cycle.
    #[error("cycle through edge {atom_a}-{atom_b} is outside radical preparation")]
    UnsupportedCycle {
        /// One edge endpoint.
        atom_a: usize,
        /// The other endpoint.
        atom_b: usize,
    },
    /// Every atom in an admitted cyclic component must be neutral carbon.
    #[error("atom {atom} in a cyclic component must be neutral carbon")]
    UnsupportedCyclicAtom {
        /// Current/source atom index, including any attached tail.
        atom: usize,
    },
    /// Every bond in an admitted cyclic component must be an ordinary single
    /// bond.
    #[error("bond {atom_a}-{atom_b} in a cyclic component must be a single bond")]
    UnsupportedCyclicBond {
        /// One current/source edge endpoint.
        atom_a: usize,
        /// The other current/source edge endpoint.
        atom_b: usize,
    },
    /// Atom stereochemistry is outside saturated cyclic-component admission.
    #[error("stereochemistry at atom {atom} in a cyclic component is unsupported")]
    UnsupportedCyclicStereo {
        /// Current/source atom index, including any attached tail.
        atom: usize,
    },
    /// Only complete tetrahedral carbon slots in cyclic components are checked.
    #[error("tetrahedral stereochemistry at atom {atom} is outside checked admission")]
    UnsupportedTetrahedralStereo {
        /// Current/source atom index.
        atom: usize,
    },
    /// Stereo-neutral graph refinement did not distinguish all four ligands.
    /// This is an exclusion, not a claim that the atom cannot be stereogenic.
    #[error("tetrahedral ligands at atom {atom} are not independently distinguished")]
    UnresolvedTetrahedralLigands {
        /// Current/source atom index.
        atom: usize,
    },
    /// Checked H and radical state did not establish saturated cyclic carbon.
    #[error(
        "unsaturated atom {atom} in a cyclic component: total connectivity \
         {total_connectivity}, radical electrons {radical_electrons}"
    )]
    UnsaturatedCyclicAtom {
        /// Current/source atom index, including any attached tail.
        atom: usize,
        /// Graph degree plus captured explicit H plus checked implicit H.
        total_connectivity: usize,
        /// The count from the versioned radical policy.
        radical_electrons: u8,
    },
    /// Graph hydrogen removal/projection is outside this capability.
    #[error("hydrogen vertex {atom} is outside radical preparation")]
    UnsupportedHydrogenVertex {
        /// Current/source hydrogen atom index.
        atom: usize,
    },
    /// Only the declared nonmetal set may have incident bonds.
    #[error("bonded element {atomic_number} at atom {atom} is unsupported")]
    UnsupportedBondedElement {
        /// Current/source atom index.
        atom: usize,
        /// Element atomic number.
        atomic_number: u8,
    },
    /// A fixed conservative guard excluded a separate cleanup family.
    #[error("unsupported preparation {operation:?} at atom {atom}")]
    UnsupportedPreparation {
        /// Current/source atom index.
        atom: usize,
        /// Excluded cleanup family.
        operation: UnsupportedPreparation,
    },
    /// The owning N-oxide cleanup rejected its requested operation.
    #[error(transparent)]
    NitrogenOxideCleanup(#[from] NitrogenOxideCleanupError),
    /// Checked accumulation exceeded its wide representation.
    #[error("valence accumulation overflow at atom {atom}")]
    ArithmeticOverflow {
        /// Current/source atom index.
        atom: usize,
    },
    /// A value could not be represented by native signed integer arithmetic.
    #[error("native integer range exceeded at atom {atom}: {value}")]
    NativeIntOutOfRange {
        /// Current/source atom index.
        atom: usize,
        /// Unnarrowed value.
        value: i128,
    },
    /// An initialized signed-byte cache would be outside 0..=127.
    #[error("unsupported native {cache:?} cache at atom {atom}: {value}")]
    NativeCacheOutOfRange {
        /// Current/source atom index.
        atom: usize,
        /// The cache being initialized.
        cache: NativeValenceCache,
        /// Unnarrowed value; no native signed narrowing is reproduced.
        value: i128,
    },
    /// The declared strict nonaromatic native valence rule rejected an atom.
    #[error("explicit valence {valence} rejected at atom {atom}")]
    ExplicitValenceRejected {
        /// Current/source atom index.
        atom: usize,
        /// Actual explicit valence.
        valence: u8,
    },
    /// The graph's H cache differs from the versioned preparation result.
    #[error("hydrogen mismatch at atom {atom}: graph {graph}, expected {expected}")]
    HydrogenStateMismatch {
        /// Current/source atom index.
        atom: usize,
        /// Recorded graph implicit H.
        graph: u8,
        /// Versioned expected implicit H.
        expected: u8,
    },
    /// Cleanup changed an indexed field outside its permitted edit contract.
    #[error("cleanup changed {field} at atom {atom}")]
    CleanupStateMismatch {
        /// Current/source atom index (zero for a node-count mismatch).
        atom: usize,
        /// The violated indexed invariant.
        field: &'static str,
    },
    /// The final radical count cannot be represented without narrowing.
    #[error("radical count out of range at atom {atom}: {value}")]
    RadicalCountOutOfRange {
        /// Current/source atom index.
        atom: usize,
        /// Unnarrowed computed count.
        value: i128,
    },
}

impl From<SmilesErrorWithSpan> for RadicalPreparationError {
    fn from(error: SmilesErrorWithSpan) -> Self {
        Self::Parse(error)
    }
}

/// Versioned atom metadata belonging to one indexed owned graph.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AtomRadicalState {
    source_no_implicit: bool,
    explicit_hydrogens: u8,
    implicit_hydrogens: u8,
    explicit_valence: u8,
    radical_electrons: u8,
    decision: RadicalDecision,
}

impl AtomRadicalState {
    /// Returns the NoImplicit flag captured from the original parser syntax.
    #[must_use]
    pub const fn source_no_implicit(&self) -> bool {
        self.source_no_implicit
    }
    /// Returns the unchanged parsed nonvertex explicit-H count.
    #[must_use]
    pub const fn explicit_hydrogens(&self) -> u8 {
        self.explicit_hydrogens
    }
    /// Returns admitted implicit H, also the native implicit-valence cache.
    #[must_use]
    pub const fn implicit_hydrogens(&self) -> u8 {
        self.implicit_hydrogens
    }
    /// Returns the admitted native explicit-valence cache in 0..=127.
    #[must_use]
    pub const fn explicit_valence(&self) -> u8 {
        self.explicit_valence
    }
    /// Returns the versioned radical count for this exact owned graph.
    #[must_use]
    pub const fn radical_electrons(&self) -> u8 {
        self.radical_electrons
    }
    /// Distinguishes assignment from preservation even when counts are equal.
    #[must_use]
    pub const fn decision(&self) -> RadicalDecision {
        self.decision
    }
}

/// An immutable fresh-parse capability for a narrow versioned radical policy.
///
/// The graph must be nonaromatic, with no hydrogen vertices.
/// [`Self::try_parse`] admits only acyclic graphs.
/// [`Self::try_parse_with_saturated_carbocycles`] additionally admits
/// saturated, neutral carbon cyclic components without explicit stereo; its
/// restrictions cover their attached tails as well. The additive
/// [`Self::try_parse_with_tetrahedral_carbocycles`] checks a conservative
/// subset of tetrahedral centers within those components. Only
/// B/C/N/O/P/S/F/Cl/Br/I may be bonded; other concrete elements may be
/// isolated. The supported N-oxide cleanup runs before checked native valence,
/// hydrogen and radical rules. Other cleanup families are conservatively
/// rejected. Every initialized native valence cache must fit 0..=127.
///
/// This certifies indexed element/isotope/charge/bond/H/radical metadata in
/// that domain, not full sanitization, general stereochemistry, CIP or
/// normalization. Radical counts remain separate from the graph's legacy
/// radical policy. There is no admission constructor for arbitrary transformed
/// graphs.
#[derive(Debug, Clone)]
pub struct RadicalPreparedSmiles {
    molecule: Smiles,
    atoms: Vec<AtomRadicalState>,
    checked_tetrahedral: bool,
}

impl RadicalPreparedSmiles {
    /// Parses, prepares and admits an acyclic graph under the versioned policy.
    ///
    /// Origin NoImplicit flags are captured before cleanup can change syntax.
    /// No legacy radical helper, aromatic localization or native library is
    /// called. All metadata is sealed together with its immutable graph.
    ///
    /// # Errors
    /// Returns a typed error for parser exclusions, unsupported domains,
    /// cleanup refusals, valence/range failures or indexed H/state mismatches.
    pub fn try_parse(text: &str) -> Result<Self, RadicalPreparationError> {
        Self::parse_with_cycles(text, CycleAdmission::Reject)
    }

    /// Also admits connected components containing saturated carbon cycles.
    ///
    /// Every atom in each cyclic component, including attached tails, must be
    /// neutral carbon without atom stereo. Every edge must be a nonaromatic
    /// [`Bond::Single`]. Checked radical counts must be zero, and graph degree
    /// plus captured nonvertex explicit H plus checked implicit H must equal
    /// four at every carbon. Isotope and atom-map labels remain supported.
    ///
    /// Other disconnected acyclic components retain [`Self::try_parse`]'s
    /// policy, including its treatment of stereo, radicals and N-oxide cleanup.
    /// Fused, spiro, bridged and cage topology is not itself excluded. This
    /// constructor does not compute a cycle basis or certify later ring
    /// perception, aromaticity, stereochemistry or full sanitization.
    ///
    /// # Errors
    /// Returns the existing preparation errors or an indexed cyclic-component
    /// error when its structural, hydrogen or radical restrictions fail.
    pub fn try_parse_with_saturated_carbocycles(
        text: &str,
    ) -> Result<Self, RadicalPreparationError> {
        Self::parse_with_cycles(text, CycleAdmission::SaturatedCarbocycles)
    }

    /// Also admits checked tetrahedral tags in saturated carbon components.
    ///
    /// The whole cyclic-component restrictions of
    /// [`Self::try_parse_with_saturated_carbocycles`] remain in force. Tags
    /// must be `@`, `@@`, `@TH1` or `@TH2`, with four graph neighbors or three
    /// graph neighbors and one bracket H. Every tagged center must belong to a
    /// cyclic component, including when the center itself is on an attached
    /// tail.
    ///
    /// An unrooted, stereo-neutral graph refinement must distinguish the four
    /// ligands. Element/isotope seeds are refined with graph neighbors and
    /// checked total H, for at most atom-count / 2 + 1 rounds. Atom maps,
    /// source syntax, stereo tags and an individualized center never
    /// distinguish them. Ties are conservatively excluded, including symmetric
    /// or relative-stereo environments that this certificate cannot resolve. No
    /// tag is removed or reoriented. This does not establish CIP labels or
    /// general sanitization.
    ///
    /// # Errors
    /// Returns the existing checked preparation errors, or an indexed error
    /// for unsupported tags/slots, tagged acyclic components or unresolved
    /// ligand equivalence. Both matching chirality policies use this same
    /// domain.
    pub fn try_parse_with_tetrahedral_carbocycles(
        text: &str,
    ) -> Result<Self, RadicalPreparationError> {
        Self::parse_with_cycles(text, CycleAdmission::TetrahedralCarbocycles)
    }

    fn parse_with_cycles(
        text: &str,
        cycles: CycleAdmission,
    ) -> Result<Self, RadicalPreparationError> {
        let parsed: Smiles = text.parse()?;
        let origins: Vec<_> = parsed
            .nodes()
            .iter()
            .map(|atom| {
                Origin {
                    no_implicit: atom.syntax() == AtomSyntax::Bracket,
                    explicit_hydrogens: atom.hydrogen_count(),
                }
            })
            .collect();
        let cyclic_atoms = check_domain(&parsed, cycles)?;
        reject_other_cleanup(&parsed)?;
        let molecule = parsed.try_with_nitrogen_oxide_cleanup()?;
        check_cleanup(&parsed, &molecule)?;
        let mut atoms = Vec::with_capacity(origins.len());
        for (id, (node, origin)) in molecule.nodes().iter().zip(&origins).enumerate() {
            let z = number(node);
            let value = explicit_sum(&molecule, id)?;
            let (explicit_valence, implicit_hydrogens) =
                property_values(id, z, node.charge_value(), origin.no_implicit, value)?;
            let graph = molecule.implicit_hydrogen_count(id);
            if graph != implicit_hydrogens {
                return Err(RadicalPreparationError::HydrogenStateMismatch {
                    atom: id,
                    graph,
                    expected: implicit_hydrogens,
                });
            }
            let radical_electrons = if origin.no_implicit {
                assign_count(
                    id,
                    z,
                    node.charge_value(),
                    explicit_valence,
                    molecule.edge_count_for_node(id),
                )?
            } else {
                0
            };
            if cyclic_atoms[id] {
                let total_connectivity = molecule
                    .edge_count_for_node(id)
                    .checked_add(usize::from(origin.explicit_hydrogens))
                    .and_then(|value| value.checked_add(usize::from(implicit_hydrogens)))
                    .ok_or(RadicalPreparationError::ArithmeticOverflow { atom: id })?;
                if total_connectivity != 4 || radical_electrons != 0 {
                    return Err(RadicalPreparationError::UnsaturatedCyclicAtom {
                        atom: id,
                        total_connectivity,
                        radical_electrons,
                    });
                }
            }
            atoms.push(AtomRadicalState {
                source_no_implicit: origin.no_implicit,
                explicit_hydrogens: origin.explicit_hydrogens,
                implicit_hydrogens,
                explicit_valence,
                radical_electrons,
                decision: if origin.no_implicit {
                    RadicalDecision::Assigned
                } else {
                    RadicalDecision::Preserved
                },
            });
        }
        // Assigned atoms have NoImplicit=true; preserved atoms retain zero
        // radicals. Consequently the final property update has identical H
        // inputs. No metadata is installed into the legacy graph sidecar.
        let checked_tetrahedral = matches!(cycles, CycleAdmission::TetrahedralCarbocycles);
        if checked_tetrahedral {
            check_tetrahedral_ligands(&molecule, &atoms)?;
        }
        Ok(Self { molecule, atoms, checked_tetrahedral })
    }

    /// Borrows the exact graph whose indexed metadata this object owns.
    #[must_use]
    pub fn molecule(&self) -> &Smiles {
        &self.molecule
    }

    /// Returns one state per graph node, with entry `i` belonging to node `i`.
    #[must_use]
    pub fn atom_states(&self) -> &[AtomRadicalState] {
        &self.atoms
    }

    /// Returns the preserved source slots for one checked tetrahedral center.
    ///
    /// Only [`Self::try_parse_with_tetrahedral_carbocycles`] certifies these
    /// centers. Earlier constructors, untagged atoms and out-of-range indices
    /// return `None`. Orientation is still read from the exact owned graph via
    /// [`Smiles::tetrahedral_slots_anticlockwise`].
    #[must_use]
    pub fn checked_tetrahedral_neighbors(&self, atom: usize) -> Option<[TetrahedralNeighbor; 4]> {
        self.checked_tetrahedral.then(|| self.molecule.tetrahedral_stereo_neighbors(atom)).flatten()
    }

    /// Consumes the capability and returns an ordinary, unqualified graph.
    ///
    /// The separate native-policy counts are discarded. Existing graph APIs
    /// continue to use their original local radical and preparation policies.
    #[must_use]
    pub fn into_unqualified_smiles(self) -> Smiles {
        self.molecule
    }
}

struct Origin {
    no_implicit: bool,
    explicit_hydrogens: u8,
}

#[derive(Clone, Copy)]
enum CycleAdmission {
    Reject,
    SaturatedCarbocycles,
    TetrahedralCarbocycles,
}

fn number(atom: &Atom) -> u8 {
    u8::from(atom.element().expect("the concrete Smiles parser excludes wildcard atoms"))
}

fn bond_value(bond: Bond) -> u64 {
    match bond {
        Bond::Single | Bond::Up | Bond::Down => 1,
        Bond::Double => 2,
        Bond::Triple => 3,
        Bond::Quadruple => 4,
    }
}

fn explicit_sum(molecule: &Smiles, atom: usize) -> Result<i128, RadicalPreparationError> {
    let mut value = u64::from(molecule.nodes()[atom].hydrogen_count());
    for edge in molecule.edges_for_node(atom) {
        value = value
            .checked_add(bond_value(edge.bond()))
            .ok_or(RadicalPreparationError::ArithmeticOverflow { atom })?;
    }
    let value = i128::from(value);
    native_int(atom, value)?;
    Ok(value)
}

fn native_int(atom: usize, value: i128) -> Result<(), RadicalPreparationError> {
    i32::try_from(value)
        .map(|_| ())
        .map_err(|_| RadicalPreparationError::NativeIntOutOfRange { atom, value })
}

fn cache_value(
    atom: usize,
    cache: NativeValenceCache,
    value: i128,
) -> Result<u8, RadicalPreparationError> {
    if !(0..=127).contains(&value) {
        return Err(RadicalPreparationError::NativeCacheOutOfRange { atom, cache, value });
    }
    Ok(u8::try_from(value).expect("the initialized cache range fits u8"))
}

fn root(parent: &[usize], mut node: usize) -> usize {
    while parent[node] != node {
        node = parent[node];
    }
    node
}

fn check_domain(
    molecule: &Smiles,
    cycles: CycleAdmission,
) -> Result<Vec<bool>, RadicalPreparationError> {
    let mut parent: Vec<_> = (0..molecule.nodes().len()).collect();
    let mut cyclic = alloc::vec![false; parent.len()];
    for (id, atom) in molecule.nodes().iter().enumerate() {
        if atom.aromatic() {
            return Err(RadicalPreparationError::UnsupportedAromaticState { atom: id });
        }
        let z = number(atom);
        if z == 1 {
            return Err(RadicalPreparationError::UnsupportedHydrogenVertex { atom: id });
        }
        if molecule.edge_count_for_node(id) != 0 && !matches!(z, 5..=9 | 15..=17 | 35 | 53) {
            return Err(RadicalPreparationError::UnsupportedBondedElement {
                atom: id,
                atomic_number: z,
            });
        }
        for edge in molecule.edges_for_node(id) {
            if edge.is_aromatic() {
                return Err(RadicalPreparationError::UnsupportedAromaticState { atom: id });
            }
            let other = edge.target();
            if id < other {
                let a = root(&parent, id);
                let b = root(&parent, other);
                if a == b {
                    if matches!(cycles, CycleAdmission::Reject) {
                        return Err(RadicalPreparationError::UnsupportedCycle {
                            atom_a: id,
                            atom_b: other,
                        });
                    }
                    cyclic[a] = true;
                } else {
                    parent[b] = a;
                    cyclic[a] |= cyclic[b];
                }
            }
        }
        // Cleanup can initialize native caches before strict validation.
        cache_value(id, NativeValenceCache::Explicit, explicit_sum(molecule, id)?)?;
    }
    // Resolve membership after all unions: a later bridge can merge a cyclic
    // component into another root, and all attached atoms share its contract.
    for id in 0..parent.len() {
        cyclic[id] = cyclic[root(&parent, id)];
    }
    for (id, atom) in molecule.nodes().iter().enumerate() {
        if !cyclic[id] {
            if matches!(cycles, CycleAdmission::TetrahedralCarbocycles)
                && atom.chirality().is_some()
            {
                return Err(RadicalPreparationError::UnsupportedTetrahedralStereo { atom: id });
            }
            continue;
        }
        if number(atom) != 6 || atom.charge_value() != 0 {
            return Err(RadicalPreparationError::UnsupportedCyclicAtom { atom: id });
        }
        if atom.chirality().is_some() {
            if !matches!(cycles, CycleAdmission::TetrahedralCarbocycles) {
                return Err(RadicalPreparationError::UnsupportedCyclicStereo { atom: id });
            }
            if !matches!(
                atom.chirality(),
                Some(Chirality::At | Chirality::AtAt | Chirality::TH(1 | 2))
            ) || molecule.tetrahedral_stereo_neighbors(id).is_none()
            {
                return Err(RadicalPreparationError::UnsupportedTetrahedralStereo { atom: id });
            }
        }
        for edge in molecule.edges_for_node(id) {
            if edge.bond() != Bond::Single {
                return Err(RadicalPreparationError::UnsupportedCyclicBond {
                    atom_a: id,
                    atom_b: edge.target(),
                });
            }
        }
    }
    Ok(cyclic)
}

// Conservative equality partition for RDKit 2026.09.1's initial legacy CIP
// assignment (Chirality.cpp buildCIPInvariants/iterateCIPRanks; BSD-3-Clause).
// Within neutral carbon/single-bond components, each native round compares
// prior ranks, twice each graph neighbor's rank and one zero per total H.
// Comparing neighbor multisets and H separately has the same equality there.
// Element/isotope seeds omit native atom-map distinctions conservatively.
// We retain the native round bound and never individualize the center: rooted
// constitutional distinctions alone do not prove native tag retention.
fn check_tetrahedral_ligands(
    molecule: &Smiles,
    atoms: &[AtomRadicalState],
) -> Result<(), RadicalPreparationError> {
    if molecule.nodes().iter().all(|atom| atom.chirality().is_none()) {
        return Ok(());
    }
    let seeds: Vec<_> = molecule
        .nodes()
        .iter()
        .map(|atom| (number(atom), atom.isotope_mass_number().unwrap_or(0)))
        .collect();
    let mut colors = equivalence_colors(&seeds);
    for _ in 0..=(molecule.nodes().len() / 2) {
        let keys: Vec<_> = colors
            .iter()
            .zip(atoms)
            .enumerate()
            .map(|(id, (&color, state))| {
                let mut neighbors: Vec<_> = molecule
                    .edges_for_node(id)
                    .map(|edge| (bond_value(edge.bond()), colors[edge.target()]))
                    .collect();
                neighbors.sort_unstable();
                let hydrogens =
                    u16::from(state.explicit_hydrogens()) + u16::from(state.implicit_hydrogens());
                (color, neighbors, hydrogens)
            })
            .collect();
        let next = equivalence_colors(&keys);
        // Prior color is part of every key: classes only split. Equal class
        // counts therefore prove stabilization, regardless of numeric ranks.
        let stable = next.iter().max() == colors.iter().max();
        colors = next;
        if stable {
            break;
        }
    }
    for (center, atom) in molecule.nodes().iter().enumerate() {
        if atom.chirality().is_none() {
            continue;
        }
        let slots = molecule
            .tetrahedral_stereo_neighbors(center)
            .ok_or(RadicalPreparationError::UnsupportedTetrahedralStereo { atom: center })?;
        let ligands = slots.map(|slot| {
            match slot {
                TetrahedralNeighbor::Atom(atom) => Some(colors[atom]),
                TetrahedralNeighbor::BracketHydrogen => None,
                TetrahedralNeighbor::LonePair => {
                    unreachable!("neutral saturated carbon has no lone pair")
                }
            }
        });
        if ligands.iter().enumerate().any(|(id, ligand)| ligands[..id].contains(ligand)) {
            return Err(RadicalPreparationError::UnresolvedTetrahedralLigands { atom: center });
        }
    }
    Ok(())
}

fn equivalence_colors<T: Ord>(keys: &[T]) -> Vec<usize> {
    let mut order: Vec<_> = (0..keys.len()).collect();
    order.sort_unstable_by(|&left, &right| keys[left].cmp(&keys[right]));
    let mut colors = alloc::vec![0; keys.len()];
    let mut color = 0;
    for pair in order.windows(2) {
        color += usize::from(keys[pair[0]] != keys[pair[1]]);
        colors[pair[1]] = color;
    }
    colors
}

fn reject_other_cleanup(molecule: &Smiles) -> Result<(), RadicalPreparationError> {
    for (id, atom) in molecule.nodes().iter().enumerate() {
        let z = number(atom);
        for edge in molecule.edges_for_node(id) {
            let other = number(&molecule.nodes()[edge.target()]);
            let operation = match (z, atom.charge_value(), edge.bond(), other) {
                (7, _, Bond::Triple, 7) => Some(UnsupportedPreparation::NitrogenTripleNitrogen),
                (15, 0, Bond::Double, 8) => Some(UnsupportedPreparation::PhosphorusOxide),
                (17 | 35 | 53, 0, Bond::Double, 8) => Some(UnsupportedPreparation::HalogenOxide),
                _ => None,
            };
            if let Some(operation) = operation {
                return Err(RadicalPreparationError::UnsupportedPreparation {
                    atom: id,
                    operation,
                });
            }
        }
    }
    Ok(())
}

fn cleanup_mismatch(atom: usize, field: &'static str) -> RadicalPreparationError {
    RadicalPreparationError::CleanupStateMismatch { atom, field }
}

fn check_cleanup(before: &Smiles, after: &Smiles) -> Result<(), RadicalPreparationError> {
    if before.nodes().len() != after.nodes().len() {
        return Err(cleanup_mismatch(0, "node count"));
    }
    if before.number_of_bonds() != after.number_of_bonds() {
        return Err(cleanup_mismatch(0, "bond count"));
    }
    for (id, (old, new)) in before.nodes().iter().zip(after.nodes()).enumerate() {
        if old.element() != new.element()
            || old.isotope_mass_number() != new.isotope_mass_number()
            || old.class() != new.class()
            || old.aromatic() != new.aromatic()
            || old.chirality() != new.chirality()
        {
            return Err(cleanup_mismatch(id, "indexed atom identity"));
        }
        if old.hydrogen_count() != new.hydrogen_count()
            || before.implicit_hydrogen_count(id) != after.implicit_hydrogen_count(id)
        {
            return Err(cleanup_mismatch(id, "recorded hydrogens"));
        }
        let changed = old.charge_value() != new.charge_value();
        if changed
            && !matches!(
                (number(old), old.charge_value(), new.charge_value()),
                (7, 0, 1) | (8, 0, -1)
            )
        {
            return Err(cleanup_mismatch(id, "formal charge"));
        }
        if (changed && new.syntax() != AtomSyntax::Bracket)
            || (!changed && old.syntax() != new.syntax())
        {
            return Err(cleanup_mismatch(id, "atom syntax"));
        }
        let mut changed_bonds = 0;
        for edge in before.edges_for_node(id) {
            let other = edge.target();
            let Some(new_edge) = after.edge_for_node_pair((id, other)) else {
                return Err(cleanup_mismatch(id, "connectivity"));
            };
            if edge.ring_num() != new_edge.ring_num() {
                return Err(cleanup_mismatch(id, "bond provenance"));
            }
            if edge.descriptor() != new_edge.descriptor() {
                let neighbor_old = &before.nodes()[other];
                let neighbor_new = &after.nodes()[other];
                if edge.bond() != Bond::Double
                    || new_edge.bond() != Bond::Single
                    || edge.is_aromatic()
                    || new_edge.is_aromatic()
                    || !matches!((number(old), number(neighbor_old)), (7, 8) | (8, 7))
                    || !changed
                    || neighbor_old.charge_value() != 0
                    || neighbor_new.charge_value() != -new.charge_value()
                {
                    return Err(cleanup_mismatch(id, "bond edit"));
                }
                changed_bonds += 1;
            }
        }
        if changed_bonds != usize::from(changed) {
            return Err(cleanup_mismatch(id, "charge/bond correspondence"));
        }
    }
    Ok(())
}

fn effective_number(z: u8, charge: i8) -> usize {
    if ELEMENTS[usize::from(z)].valences == [-1] {
        usize::from(z)
    } else {
        usize::try_from((i16::from(z) - i16::from(charge)).clamp(0, 118)).unwrap()
    }
}

fn hypervalent(z: u8, effective: usize) -> bool {
    (matches!(z, 15 | 16) && effective > 16) || (matches!(z, 33 | 34) && effective > 34)
}

fn property_values(
    id: usize,
    z: u8,
    charge: i8,
    no_implicit: bool,
    value: i128,
) -> Result<(u8, u8), RadicalPreparationError> {
    native_int(id, value)?;
    let explicit = cache_value(id, NativeValenceCache::Explicit, value)?;
    let original = ELEMENTS[usize::from(z)].valences;
    let effective = effective_number(z, charge);
    let valences = ELEMENTS[effective].valences;
    let mut maximum = i128::from(*valences.last().unwrap());
    let mut offset = 0;
    if hypervalent(z, effective) {
        maximum = i128::from(*original.last().unwrap());
        offset = -i128::from(charge);
    }
    // Kept faithful to the arithmetic rule even though H vertices are not
    // admitted by the public constructor.
    if z == 1 && charge == -1 {
        maximum = 2;
    }
    native_int(id, value + offset)?;
    if maximum >= 0 && *original.last().unwrap() >= 0 && value + offset > maximum {
        return Err(RadicalPreparationError::ExplicitValenceRejected {
            atom: id,
            valence: explicit,
        });
    }
    if no_implicit || effective == 0 || valences[0] == -1 {
        return Ok((explicit, 0));
    }
    // Parser radical counts are zero. This property route is only used with
    // zero prior counts; no generic transformed-state inference is exposed.
    let (allowed, explicit_plus_radical) = if hypervalent(z, effective) {
        (original, value - i128::from(charge))
    } else {
        (valences, value)
    };
    native_int(id, explicit_plus_radical)?;
    let implicit = allowed.iter().take_while(|&&v| v >= 0).find_map(|&v| {
        (i128::from(v) >= explicit_plus_radical).then_some(i128::from(v) - explicit_plus_radical)
    });
    let implicit = match implicit {
        Some(value) => value,
        None if *allowed.last().unwrap() != -1 && *original.last().unwrap() > 0 => {
            return Err(RadicalPreparationError::ExplicitValenceRejected {
                atom: id,
                valence: explicit,
            });
        }
        None => 0,
    };
    Ok((explicit, cache_value(id, NativeValenceCache::Implicit, implicit)?))
}

fn assign_count(
    id: usize,
    z: u8,
    charge: i8,
    valence: u8,
    degree: usize,
) -> Result<u8, RadicalPreparationError> {
    let entry = &ELEMENTS[usize::from(z)];
    let outer = i128::from(entry.outer);
    let charge = i128::from(charge);
    let value = if entry.valences == [-1] {
        if degree != 0 { 0 } else { (outer - charge).max(0) % 2 }
    } else {
        let total = i128::from(valence);
        let base = if matches!(z, 1 | 2) { 2 } else { 8 };
        let mut count = base - outer - total + charge;
        if count < 0 {
            count = 0;
            if entry.valences.len() > 1 {
                for &allowed in entry.valences {
                    let candidate = i128::from(allowed) - total + charge;
                    if candidate >= 0 {
                        count = candidate;
                        break;
                    }
                }
            }
        }
        let early = outer - total - charge;
        if early >= 0 {
            count = count.min(early);
        }
        count
    };
    native_int(id, value)?;
    u8::try_from(value)
        .map_err(|_| RadicalPreparationError::RadicalCountOutOfRange { atom: id, value })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn native_capability_does_not_replace_the_legacy_graph_policy() {
        let original: Smiles = "[Fe+2]".parse().unwrap();
        let before = format!("{original:?}");
        assert_eq!(original.radical_electron_counts().as_ref(), [2]);
        let prepared = RadicalPreparedSmiles::try_parse("[Fe+2]").unwrap();
        assert_eq!(prepared.atom_states()[0].radical_electrons(), 0);
        assert_eq!(prepared.molecule().radical_electron_counts().as_ref(), [2]);
        assert_eq!(format!("{original:?}"), before);
        assert_eq!(prepared.into_unqualified_smiles().radical_electron_counts().as_ref(), [2]);
    }

    #[test]
    fn initialized_cache_bounds_do_not_reproduce_signed_narrowing() {
        for cache in [NativeValenceCache::Explicit, NativeValenceCache::Implicit] {
            for value in [0, 1, 126, 127] {
                assert_eq!(i128::from(cache_value(0, cache, value).unwrap()), value);
            }
            for value in [-1, 128, 255, 256, 2_147_483_647, 2_147_483_648] {
                assert!(matches!(
                    cache_value(0, cache, value),
                    Err(RadicalPreparationError::NativeCacheOutOfRange { .. })
                ));
            }
        }
    }

    #[test]
    fn native_integer_guard_is_distinct_from_cache_capacity() {
        assert!(native_int(0, i128::from(i32::MAX)).is_ok());
        assert!(matches!(
            native_int(0, i128::from(i32::MAX) + 1),
            Err(RadicalPreparationError::NativeIntOutOfRange { .. })
        ));
        assert!(native_int(0, i128::from(i32::MIN)).is_ok());
        assert!(native_int(0, i128::from(i32::MIN) - 1).is_err());
    }

    #[test]
    fn native_effective_zero_and_unpreferred_elements_use_distinct_rules() {
        assert_eq!(effective_number(6, 15), 0);
        assert_eq!(effective_number(26, 15), 26);
        assert_eq!(ELEMENTS[0].valences, [-1]);
        assert_eq!(ELEMENTS[0].outer, 0);
        assert_eq!(assign_count(0, 26, 2, 0, 0).unwrap(), 0);
        assert_eq!(assign_count(0, 26, 3, 0, 0).unwrap(), 1);
        assert_eq!(assign_count(0, 26, 3, 1, 1).unwrap(), 0);
    }
}
