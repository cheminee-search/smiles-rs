//! Typed neighbors for mapped tetrahedral orientation.

use elements_rs::Element;

use super::{Smiles, SmilesAtomPolicy, StereoNeighbor, WildcardSmiles};
use crate::atom::bracketed::chirality::Chirality;

/// Identity of one position around a tetrahedral stereocenter.
///
/// A graph hydrogen, including an isotope or charged hydrogen, is an
/// [`Atom`](Self::Atom). A lone pair is never represented as a hydrogen.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TetrahedralNeighbor {
    /// A neighboring graph atom, named by its node index.
    Atom(usize),
    /// A hydrogen written inside the center's brackets, without a graph index.
    BracketHydrogen,
    /// The virtual fourth position of a supported three-neighbor sulfur atom.
    LonePair,
}

impl<AtomPolicy: SmilesAtomPolicy> Smiles<AtomPolicy> {
    /// Returns four typed tetrahedral neighbors in parsed reference order.
    ///
    /// Supports `@`, `@@`, `@TH1`, and `@TH2` with four represented neighbors.
    /// A terminal lone-pair position is also supported for nonaromatic sulfur
    /// with three graph neighbors, nonaromatic bonds of total order four, and
    /// no bracket hydrogen. Formal charge is not an additional restriction.
    /// Valence-three sulfonium, other lone-pair elements, and other
    /// coordination classes are outside this lone-pair extension.
    ///
    /// This reports parsed slots, not canonical ranking, stereogenicity, or
    /// chemical validity. It does not sanitize the molecule. Returns `None`
    /// for absent/unsupported markup, invalid centers, or incomplete or
    /// non-distinct slot identities.
    ///
    /// # Examples
    ///
    /// ```
    /// use smiles_rs::{
    ///     Smiles,
    ///     TetrahedralNeighbor::{Atom, LonePair},
    /// };
    /// let molecule: Smiles = "CC[S@](=O)CCC".parse()?;
    /// assert_eq!(
    ///     molecule.tetrahedral_stereo_neighbors(2),
    ///     Some([Atom(1), Atom(3), Atom(4), LonePair])
    /// );
    /// # Ok::<(), smiles_rs::SmilesErrorWithSpan>(())
    /// ```
    #[must_use]
    pub fn tetrahedral_stereo_neighbors(&self, atom_id: usize) -> Option<[TetrahedralNeighbor; 4]> {
        if !matches!(
            self.node_by_id(atom_id)?.chirality()?,
            Chirality::At | Chirality::AtAt | Chirality::TH(1 | 2)
        ) {
            return None;
        }
        let parsed = self.parsed_stereo_neighbors_row(atom_id);
        if parsed.len() != 4 && !(parsed.len() == 3 && self.has_tetrahedral_lone_pair(atom_id)) {
            return None;
        }
        let mut neighbors = [TetrahedralNeighbor::LonePair; 4];
        for (index, neighbor) in parsed.iter().enumerate() {
            neighbors[index] = match *neighbor {
                StereoNeighbor::Atom(atom) => TetrahedralNeighbor::Atom(atom),
                StereoNeighbor::ExplicitHydrogen => TetrahedralNeighbor::BracketHydrogen,
            };
        }
        if neighbors
            .iter()
            .enumerate()
            .any(|(index, neighbor)| neighbors[..index].contains(neighbor))
        {
            return None;
        }
        Some(neighbors)
    }

    /// Reads tetrahedral orientation relative to caller-named typed slots.
    ///
    /// Looking from `neighbors[0]`, `Some(true)` means anticlockwise (`@` or
    /// `@TH1`); `Some(false)` means clockwise. All four identities returned by
    /// [`Self::tetrahedral_stereo_neighbors`] must appear exactly once.
    /// A bracket hydrogen and a lone pair are distinct identities.
    ///
    /// Returns `None` for unsupported centers or a nonmatching slot set. This
    /// reports parsed orientation without establishing stereogenicity.
    ///
    /// # Examples
    ///
    /// ```
    /// use smiles_rs::{
    ///     Smiles,
    ///     TetrahedralNeighbor::{Atom, LonePair},
    /// };
    /// let molecule: Smiles = "CC[S@](=O)CCC".parse()?;
    /// let neighbors = [Atom(1), Atom(3), Atom(4), LonePair];
    /// assert_eq!(molecule.tetrahedral_slots_anticlockwise(2, neighbors), Some(true));
    /// let mut swapped = neighbors;
    /// swapped.swap(1, 2);
    /// assert_eq!(molecule.tetrahedral_slots_anticlockwise(2, swapped), Some(false));
    /// # Ok::<(), smiles_rs::SmilesErrorWithSpan>(())
    /// ```
    #[must_use]
    pub fn tetrahedral_slots_anticlockwise(
        &self,
        atom_id: usize,
        neighbors: [TetrahedralNeighbor; 4],
    ) -> Option<bool> {
        let parsed = self.tetrahedral_stereo_neighbors(atom_id)?;
        let mut anticlockwise =
            matches!(self.node_by_id(atom_id)?.chirality()?, Chirality::At | Chirality::TH(1));
        let mut order = [0; 4];
        let mut seen = 0_u8;
        for (slot, neighbor) in neighbors.into_iter().enumerate() {
            let index = parsed.iter().position(|&parsed| parsed == neighbor)?;
            if seen & (1 << index) != 0 {
                return None;
            }
            seen |= 1 << index;
            for &previous in &order[..slot] {
                anticlockwise ^= previous > index;
            }
            order[slot] = index;
        }
        Some(anticlockwise)
    }

    // The valence-four sulfur subset of RDKit 2026.09.1's tetrahedral-center
    // eligibility. This does not generalize the ordering rules of sulfonium,
    // phosphorus, selenium, or nitrogen centers.
    pub(super) fn has_tetrahedral_lone_pair(&self, atom_id: usize) -> bool {
        let Some(atom) = self.node_by_id(atom_id) else {
            return false;
        };
        let parsed = self.parsed_stereo_neighbors_row(atom_id);
        atom.element() == Some(Element::S)
            && !atom.aromatic()
            && atom.hydrogen_count() == 0
            && self.edge_count_for_node(atom_id) == 3
            && parsed.len() == 3
            && parsed.iter().all(|neighbor| matches!(neighbor, StereoNeighbor::Atom(_)))
            && self.edges_for_node(atom_id).all(|edge| !edge.is_aromatic())
            && super::implicit_hydrogens::explicit_valence(self, atom_id) == 4
    }
}

impl WildcardSmiles {
    /// See [`Smiles::tetrahedral_stereo_neighbors`].
    #[inline]
    #[must_use]
    pub fn tetrahedral_stereo_neighbors(&self, atom_id: usize) -> Option<[TetrahedralNeighbor; 4]> {
        self.inner.tetrahedral_stereo_neighbors(atom_id)
    }

    /// See [`Smiles::tetrahedral_slots_anticlockwise`].
    #[inline]
    #[must_use]
    pub fn tetrahedral_slots_anticlockwise(
        &self,
        atom_id: usize,
        neighbors: [TetrahedralNeighbor; 4],
    ) -> Option<bool> {
        self.inner.tetrahedral_slots_anticlockwise(atom_id, neighbors)
    }
}
