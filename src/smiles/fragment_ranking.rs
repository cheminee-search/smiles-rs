//! Numeric candidates for default fragment-parent selection.

use alloc::vec::Vec;

use super::{Smiles, SmilesAtomPolicy};
use crate::atom::Atom;

/// One connected component tied for the maximal fragment-parent numeric key.
///
/// Multiple candidates remain unresolved: no canonical-string tie breaker has
/// been applied. The pseudo-weight is a ranking score, not molecular mass.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FragmentCandidate {
    component_id: usize,
    atom_count: u32,
    pseudo_weight: f64,
}

impl FragmentCandidate {
    /// Returns the identifier used by [`Smiles::connected_components`].
    #[inline]
    #[must_use]
    pub const fn component_id(&self) -> usize {
        self.component_id
    }

    /// Returns graph vertices plus bracket and cached inferred hydrogens.
    /// Graph hydrogen atoms are already counted as vertices.
    #[inline]
    #[must_use]
    pub const fn atom_count_including_hydrogens(&self) -> u32 {
        self.atom_count
    }

    /// Returns the secondary ranking score, not a physical molecular mass.
    /// See [`Smiles::fragment_parent_candidates`] for its arithmetic policy.
    #[inline]
    #[must_use]
    pub const fn pseudo_weight(&self) -> f64 {
        self.pseudo_weight
    }

    fn add_atom(
        &mut self,
        atom_id: usize,
        atom: Atom,
        inferred_hydrogens: u8,
    ) -> Result<(), FragmentRankingError> {
        let element = atom.element().ok_or(FragmentRankingError::UnknownElement { atom_id })?;
        let hydrogens = u32::from(atom.hydrogen_count()) + u32::from(inferred_hydrogens);
        self.atom_count = self
            .atom_count
            .checked_add(1 + hydrogens)
            .ok_or(FragmentRankingError::AtomCountOverflow { component_id: self.component_id })?;

        // Preserve the native source's two separately rounded additions. Do
        // not combine H with the atom term or accumulate charge separately.
        let atom_number_term =
            100 * u32::from(u8::from(element)) + u32::from(atom.isotope_mass_number().unwrap_or(0));
        self.pseudo_weight += f64::from(atom_number_term) - f64::from(atom.charge_value()) * 0.1;
        self.pseudo_weight += f64::from(hydrogens);
        Ok(())
    }
}

/// A graph cannot be ranked within the concrete-element, checked-count policy.
/// No partial candidate set is returned on error.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum FragmentRankingError {
    /// An atom has no concrete element and hence no supported numeric score.
    #[error("Atom {atom_id} has no concrete element for fragment ranking")]
    UnknownElement {
        /// Original parent atom identifier.
        atom_id: usize,
    },
    /// A component's inclusive-H atom count exceeds the native `u32` range.
    /// Native unsigned wraparound is deliberately unsupported.
    #[error("Inclusive hydrogen atom count exceeds u32 for component {component_id}")]
    AtomCountOverflow {
        /// Connected-component identifier.
        component_id: usize,
    },
}

impl<AtomPolicy: SmilesAtomPolicy> Smiles<AtomPolicy> {
    /// Returns every component with the maximal default fragment-parent
    /// numeric key, in connected-component order.
    ///
    /// This follows the numeric policy of the
    /// [RDKit 2026.09.1 chooser](https://github.com/rdkit/rdkit/blob/fece8caa860bdf6c9c82bdb9253f22a0cf87150c/Code/GraphMol/MolStandardize/Fragment.cpp):
    /// no organic preference, atom counting enabled, and hydrogens included.
    /// First maximize graph vertices plus bracket and cached inferred H.
    /// Then maximize a pseudo-weight accumulated in increasing parent atom ID
    /// order. For each atom, separately add `100 * atomic_number + isotope -
    /// 0.1 * charge` and its nonvertex H count using ordinary `f64` arithmetic.
    /// Exact floating-point comparisons determine numeric ties; this score
    /// is not molecular mass and can depend on atom order and H representation.
    ///
    /// The graph is used as supplied, without sanitization, H removal,
    /// canonicalization or extraction. Multiple results are an unresolved
    /// native canonical-string tie, including identical-looking components.
    /// No native or local canonical-string tie breaker is applied. An empty
    /// graph returns no candidates; a single component still has its numeric
    /// key computed. This is not the complete native fragment-parent operation.
    ///
    /// # Errors
    /// Returns [`FragmentRankingError::UnknownElement`] for an unknown element
    /// anywhere in the graph, or [`FragmentRankingError::AtomCountOverflow`]
    /// rather than wrapping a component's inclusive-H count beyond `u32`.
    ///
    /// # Examples
    /// ```
    /// use smiles_rs::Smiles;
    /// let parent: Smiles = "O.N".parse().expect("valid SMILES");
    /// let candidates = parent.fragment_parent_candidates()?;
    /// assert_eq!(candidates.len(), 1);
    /// assert_eq!(candidates[0].component_id(), 1);
    /// assert_eq!(candidates[0].atom_count_including_hydrogens(), 4);
    /// # Ok::<(), smiles_rs::FragmentRankingError>(())
    /// ```
    #[allow(clippy::float_cmp)] // Exact equality defines the source policy's numeric ties.
    pub fn fragment_parent_candidates(
        &self,
    ) -> Result<Vec<FragmentCandidate>, FragmentRankingError> {
        let components = self.connected_components();
        let mut candidates = (0..components.number_of_components())
            .map(|component_id| {
                FragmentCandidate { component_id, atom_count: 0, pseudo_weight: 0.0 }
            })
            .collect::<Vec<_>>();
        // One global ordered pass preserves each component's ordered atom
        // subsequence, without extracting or repeatedly scanning components.
        for (atom_id, &atom) in self.nodes().iter().enumerate() {
            candidates[components.component_of_node(atom_id)].add_atom(
                atom_id,
                atom,
                self.implicit_hydrogen_count(atom_id),
            )?;
        }

        let Some(mut best) = candidates.first().copied() else {
            return Ok(candidates);
        };
        for &candidate in &candidates[1..] {
            if candidate.atom_count > best.atom_count
                || (candidate.atom_count == best.atom_count
                    && candidate.pseudo_weight > best.pseudo_weight)
            {
                best = candidate;
            }
        }
        candidates.retain(|candidate| {
            candidate.atom_count == best.atom_count && candidate.pseudo_weight == best.pseudo_weight
        });
        Ok(candidates)
    }
}

#[cfg(test)]
mod tests {
    use super::{FragmentCandidate, FragmentRankingError, Smiles};
    use crate::WildcardSmiles;

    #[test]
    fn unknown_elements_fail_the_entire_ranking() {
        let parent: WildcardSmiles = "*.C".parse().unwrap();
        assert_eq!(
            parent.inner.fragment_parent_candidates(),
            Err(FragmentRankingError::UnknownElement { atom_id: 0 })
        );
    }

    #[test]
    fn count_bound_is_checked_instead_of_wrapping() {
        let parent: Smiles = "[NH4+].C".parse().unwrap();
        let atom = parent.nodes()[0];
        let mut score =
            FragmentCandidate { component_id: 7, atom_count: u32::MAX - 5, pseudo_weight: 0.0 };
        score.add_atom(0, atom, 0).unwrap();
        assert_eq!(score.atom_count, u32::MAX);
        assert_eq!(score.pseudo_weight.to_bits(), 703.9_f64.to_bits());
        assert_eq!(
            score.add_atom(0, atom, 0),
            Err(FragmentRankingError::AtomCountOverflow { component_id: 7 })
        );
        assert_eq!(score.atom_count, u32::MAX);
        assert_eq!(score.pseudo_weight.to_bits(), 703.9_f64.to_bits());
    }
}
