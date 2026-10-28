//! Immutable parsed-graph admission and versioned radical assignment.

use alloc::vec::Vec;

use thiserror::Error;

use super::element_data::ELEMENTS;
use crate::{
    NitrogenOxideCleanupError, Smiles, SmilesErrorWithSpan,
    atom::{Atom, AtomSyntax},
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
    /// The graph contains a cycle; lexical ring labels alone do not imply one.
    #[error("cycle through edge {atom_a}-{atom_b} is outside radical preparation")]
    UnsupportedCycle {
        /// One edge endpoint.
        atom_a: usize,
        /// The other endpoint.
        atom_b: usize,
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
/// The graph must be acyclic and nonaromatic, with no hydrogen vertices.
/// Only B/C/N/O/P/S/F/Cl/Br/I may be bonded; other concrete elements may be
/// isolated. The supported N-oxide cleanup runs before checked native valence,
/// hydrogen and radical rules. Other cleanup families are conservatively
/// rejected. Every initialized native valence cache must fit 0..=127.
///
/// This certifies indexed element/isotope/charge/bond/H/radical metadata in
/// that domain, not full sanitization, stereochemistry, CIP or normalization.
/// Radical counts remain separate from the graph's legacy radical policy.
/// There is no admission constructor for arbitrary transformed graphs.
#[derive(Debug, Clone)]
pub struct RadicalPreparedSmiles {
    molecule: Smiles,
    atoms: Vec<AtomRadicalState>,
}

impl RadicalPreparedSmiles {
    /// Parses, prepares and admits an owned graph under the versioned policy.
    ///
    /// Origin NoImplicit flags are captured before cleanup can change syntax.
    /// No legacy radical helper, aromatic localization or native library is
    /// called. All metadata is sealed together with its immutable graph.
    ///
    /// # Errors
    /// Returns a typed error for parser exclusions, unsupported domains,
    /// cleanup refusals, valence/range failures or indexed H/state mismatches.
    pub fn try_parse(text: &str) -> Result<Self, RadicalPreparationError> {
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
        check_domain(&parsed)?;
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
        Ok(Self { molecule, atoms })
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

fn check_domain(molecule: &Smiles) -> Result<(), RadicalPreparationError> {
    let mut parent: Vec<_> = (0..molecule.nodes().len()).collect();
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
                    return Err(RadicalPreparationError::UnsupportedCycle {
                        atom_a: id,
                        atom_b: other,
                    });
                }
                parent[b] = a;
            }
        }
        // Cleanup can initialize native caches before strict validation.
        cache_value(id, NativeValenceCache::Explicit, explicit_sum(molecule, id)?)?;
    }
    Ok(())
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
