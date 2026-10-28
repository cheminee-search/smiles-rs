//! Optional weights of the current recorded graph using RDKit 2026.09.1 tables.
//!
//! This module projects concrete atoms in their stored order into the
//! separately versioned molecular-formulas calculator. It uses bracket H and
//! the existing recorded implicit-H cache; it does not sanitize, infer H,
//! prepare a molecule, change atom order, or convert through an aggregated
//! molecular formula. Existing scientific formula conversions and weight
//! methods are unchanged.

use alloc::vec::Vec;
use core::num::NonZeroU16;

use molecular_formulas_native::rdkit_2026_09_1::{AtomMassInput, molecular_weights};
pub use molecular_formulas_native::rdkit_2026_09_1::{MolecularWeightError, MolecularWeights};

use crate::{Smiles, atom::Atom};

/// A recorded graph cannot be projected or is outside the calculator's domain.
///
/// Errors refer to the current graph's zero-based atom order. No partial weight
/// result is returned. Projection errors precede calculator errors because all
/// records are collected before their weights are calculated.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum RecordedWeightError {
    /// The graph does not have one recorded implicit-H count per atom.
    #[error("recorded hydrogen cache has {hydrogens} entries for {atoms} atoms")]
    HydrogenCacheLength {
        /// Number of graph atoms.
        atoms: usize,
        /// Number of recorded implicit-H entries.
        hydrogens: usize,
    },
    /// An atom lacks the concrete element required by the mass calculator.
    #[error("atom {atom_index} has no concrete element")]
    MissingElement {
        /// Zero-based atom position in the current graph.
        atom_index: usize,
    },
    /// A recorded isotope number is zero, which is not natural abundance.
    #[error("atom {atom_index} records isotope zero")]
    ZeroIsotope {
        /// Zero-based atom position in the current graph.
        atom_index: usize,
    },
    /// The projected atom or H counts exceed the checked native-width domain.
    #[error("recorded-weight calculation failed: {0}")]
    Calculation(#[source] MolecularWeightError),
}

/// Computes native-vintage weights from this graph's ordered recorded atoms.
///
/// Every graph vertex contributes one record, including ordinary, isotopic and
/// charged hydrogen atoms. Only bracket H plus the current recorded implicit H
/// contribute additional nonvertex hydrogens. Graph-H neighbors are not added
/// a second time. Isotope numbers and formal charges are preserved.
///
/// This borrows the graph without changing it and allocates a temporary record
/// vector. It performs no sanitization, aromaticity assignment, radical
/// inference, H normalization, fragment selection or formula reparse. A
/// successful calculation does not certify a prepared descriptor contract.
/// Evaluate transformed graphs separately: H representation and atom order can
/// change floating-point rounding even when their compositions agree.
///
/// Average weight is conventionally expressed in g/mol and exact weight in Da;
/// the calculator's isotope fallback and electron correction are documented on
/// [`MolecularWeights`]. Existing scientific mass conventions are unaffected.
///
/// # Errors
/// Returns [`RecordedWeightError`] for a cache-length mismatch, a nonconcrete
/// atom, isotope zero, or a native-width calculation bound. Shape is checked
/// first, then the `u32` atom-count bound, then each element/isotope in graph
/// order. All records are projected before the calculator checks H bounds.
/// Neither unprojectable atoms nor failed calculations are silently omitted.
///
/// # Examples
/// ```
/// use smiles_rs::{
///     Smiles,
///     rdkit_2026_09_1::{RecordedWeightError, recorded_weights},
/// };
/// let molecule: Smiles = "CCO".parse().expect("valid declared SMILES");
/// let weights = recorded_weights(&molecule)?;
/// assert!(weights.average() > 0.0);
/// assert!(weights.exact() > 0.0);
/// # Ok::<(), RecordedWeightError>(())
/// ```
pub fn recorded_weights(molecule: &Smiles) -> Result<MolecularWeights, RecordedWeightError> {
    let records = project_records(molecule.nodes(), molecule.implicit_hydrogen_counts())?;
    molecular_weights(&records).map_err(RecordedWeightError::Calculation)
}

fn validate_shape(atoms: usize, hydrogens: usize) -> Result<(), RecordedWeightError> {
    if atoms != hydrogens {
        return Err(RecordedWeightError::HydrogenCacheLength { atoms, hydrogens });
    }
    u32::try_from(atoms)
        .map(|_| ())
        .map_err(|_| RecordedWeightError::Calculation(MolecularWeightError::AtomCountOverflow))
}

fn project_records(
    atoms: &[Atom],
    hydrogens: &[u8],
) -> Result<Vec<AtomMassInput>, RecordedWeightError> {
    validate_shape(atoms.len(), hydrogens.len())?;
    let mut records = Vec::with_capacity(atoms.len());
    for (atom_index, atom) in atoms.iter().enumerate() {
        let element = atom.element().ok_or(RecordedWeightError::MissingElement { atom_index })?;
        let isotope = match atom.isotope_mass_number() {
            Some(number) => {
                Some(
                    NonZeroU16::new(number)
                        .ok_or(RecordedWeightError::ZeroIsotope { atom_index })?,
                )
            }
            None => None,
        };
        records.push(AtomMassInput {
            element,
            isotope,
            formal_charge: atom.charge_value(),
            nonvertex_hydrogens: u16::from(atom.hydrogen_count())
                + u16::from(hydrogens[atom_index]),
        });
    }
    Ok(records)
}

#[cfg(test)]
mod tests;
