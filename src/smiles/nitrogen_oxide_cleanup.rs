//! An opt-in, bounded nitrogen-oxide charge-separation operation.
//!
//! The chemical rule follows RDKit 2026.09.1 `MolOps::cleanUp`'s neutral
//! explicit-valence-five nitrogen/double-oxygen branch. This is an independent
//! implementation, not the other cleanup rules or a sanitization pipeline.

use alloc::{boxed::Box, vec::Vec};

use elements_rs::Element;
use geometric_traits::traits::{SparseValuedMatrix2DRef, SparseValuedMatrixRef};
use thiserror::Error;

use super::{BondEntry, BondMatrix, Smiles};
use crate::{atom::Atom, bond::Bond};

/// An unsupported state encountered by nitrogen-oxide charge separation.
///
/// Errors apply to the whole operation; the borrowed input remains unchanged.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
#[non_exhaustive]
pub enum NitrogenOxideCleanupError {
    /// A potentially affected nitrogen has aromatic or unsupported bond typing.
    #[error("unsupported nitrogen-oxide bond environment at atom {nitrogen}")]
    UnsupportedEnvironment {
        /// Current graph atom ID of the nitrogen.
        nitrogen: usize,
    },
    /// The first eligible oxygen cannot be established from stored bond order.
    #[error("ambiguous nitrogen-oxide oxygen order at atom {nitrogen}")]
    AmbiguousOxygenOrder {
        /// Current graph atom ID of the nitrogen.
        nitrogen: usize,
    },
    /// A double-bonded oxygen is nonterminal, aromatic or bears nonvertex H.
    #[error("unsupported nitrogen-oxide oxygen at atom {oxygen}")]
    UnsupportedOxygen {
        /// Current graph atom ID of the oxygen.
        oxygen: usize,
    },
    /// An affected nitrogen has a nonzero inferred-H override.
    #[error("unsupported nitrogen-oxide hydrogen state at atom {atom}")]
    UnsupportedHydrogens {
        /// Current graph atom ID of the affected atom.
        atom: usize,
    },
    /// An affected atom carries a stereo tag or parsed stereo-neighbor row.
    #[error("unsupported nitrogen-oxide stereochemistry at atom {atom}")]
    UnsupportedStereochemistry {
        /// Current graph atom ID of the affected atom.
        atom: usize,
    },
    /// A present radical sidecar assigns nonzero electrons to an affected atom.
    #[error("stored nitrogen-oxide radical electrons at atom {atom}")]
    StoredRadicalElectrons {
        /// Current graph atom ID of the affected atom.
        atom: usize,
    },
    /// A saved source cannot receive the same indexed edits without ambiguity.
    #[error("nitrogen-oxide saved-source conflict at depth {depth}")]
    SavedSourceConflict {
        /// Source depth, with one denoting the immediate saved source.
        depth: usize,
    },
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct OxideEdit {
    nitrogen: usize,
    oxygen: usize,
}

impl Smiles {
    /// Returns an owned graph with supported neutral nitrogen oxides separated
    /// into `N+`/`O-` joined by a single bond.
    ///
    /// The rule requires nonaromatic neutral nitrogen with explicit valence
    /// five (graph bonds plus bracket H, excluding inferred H). Its first
    /// neutral double-oxygen neighbor must be terminal and have no nonvertex
    /// H. Other incident bonds must be single or neutral terminal double-O.
    /// Stored order selects among non-ring candidates; a sole closure candidate
    /// is supported, but multiple closure-only candidates are ambiguous.
    /// Selection precedes oxygen support checks; an unsupported first oxygen
    /// is never skipped. This representation-sensitive rule does not certify
    /// native adjacency order for arbitrary transformed graphs.
    ///
    /// Atom IDs, isotopes, maps, per-atom H, remote graph-relative stereo and
    /// existing radical sidecars are preserved. An absent radical sidecar stays
    /// absent. Every saved Kekule source receives the same validated edits.
    /// Affected atoms become bracket atoms so their charges survive rendering.
    /// Source copying is iterative, including the unchanged no-match case.
    ///
    /// This is **not sanitization**: no H removal, aromaticity, radical or CIP
    /// assignment, or other nitrogen/phosphorus/halogen cleanup is performed.
    /// Unrelated chemistry is unchanged and uncertified. Prepared targets and
    /// fingerprints must be computed afresh from the returned graph.
    ///
    /// # Errors
    /// Returns [`NitrogenOxideCleanupError`] for an unsupported potential
    /// mutation or inconsistent saved source. The input is never mutated.
    ///
    /// # Examples
    /// ```
    /// use smiles_rs::{NitrogenOxideCleanupError, Smiles};
    /// let original: Smiles = "O=N(C)(C)C".parse().expect("valid syntax");
    /// let cleaned = original.try_with_nitrogen_oxide_cleanup()?;
    /// assert_eq!(original.nodes()[0].charge_value(), 0);
    /// assert_eq!(cleaned.nodes()[0].charge_value(), -1);
    /// assert_eq!(cleaned.nodes()[1].charge_value(), 1);
    /// # Ok::<(), NitrogenOxideCleanupError>(())
    /// ```
    pub fn try_with_nitrogen_oxide_cleanup(&self) -> Result<Self, NitrogenOxideCleanupError> {
        let mut edits = Vec::new();
        for nitrogen in 0..self.atom_nodes.len() {
            if let Some(edit) = self.nitrogen_oxide_edit(nitrogen)? {
                edits.push(edit);
            }
        }

        let mut layers = Vec::new();
        let mut next = self.kekulization_source.as_deref();
        while let Some(layer) = next {
            layers.push(layer);
            next = layer.kekulization_source.as_deref();
        }
        for (index, layer) in layers.iter().enumerate() {
            if !self.cleanup_source_agrees(layer, &edits) {
                return Err(NitrogenOxideCleanupError::SavedSourceConflict { depth: index + 1 });
            }
        }

        let mut source = None;
        for layer in layers.into_iter().rev() {
            let mut cleaned = layer.apply_nitrogen_oxide_edits(&edits);
            cleaned.kekulization_source = source;
            source = Some(Box::new(cleaned));
        }
        let mut cleaned = self.apply_nitrogen_oxide_edits(&edits);
        cleaned.kekulization_source = source;
        Ok(cleaned)
    }

    fn nitrogen_oxide_edit(
        &self,
        nitrogen: usize,
    ) -> Result<Option<OxideEdit>, NitrogenOxideCleanupError> {
        let atom = self.atom_nodes[nitrogen];
        if atom.element() != Some(Element::N) || atom.charge_value() != 0 {
            return Ok(None);
        }
        let neighbors = || self.bond_matrix.sparse_row_entries(nitrogen);
        let candidates = neighbors()
            .filter(|(neighbor, entry)| {
                entry.bond() == Bond::Double
                    && self.atom_nodes[*neighbor].element() == Some(Element::O)
                    && self.atom_nodes[*neighbor].charge_value() == 0
            })
            .collect::<Vec<_>>();
        if candidates.is_empty() {
            return Ok(None);
        }
        if atom.aromatic() || neighbors().any(|(_, entry)| entry.aromatic()) {
            return Err(NitrogenOxideCleanupError::UnsupportedEnvironment { nitrogen });
        }
        // Only equality with five matters; saturation cannot turn overflow
        // into a false five-valent match.
        let valence = neighbors().fold(atom.hydrogen_count(), |total, (_, entry)| {
            total.saturating_add(match entry.bond() {
                Bond::Single | Bond::Up | Bond::Down => 1,
                Bond::Double => 2,
                Bond::Triple => 3,
                Bond::Quadruple => 4,
            })
        });
        if valence != 5 {
            return Ok(None);
        }
        let oxygen = first_oxygen(nitrogen, &candidates)?;
        self.validate_cleanup_oxygen(oxygen)?;
        for (neighbor, entry) in neighbors() {
            match entry.bond() {
                Bond::Single | Bond::Up | Bond::Down => {}
                Bond::Double
                    if self.atom_nodes[neighbor].element() == Some(Element::O)
                        && self.atom_nodes[neighbor].charge_value() == 0 =>
                {
                    self.validate_cleanup_oxygen(neighbor)?;
                }
                _ => return Err(NitrogenOxideCleanupError::UnsupportedEnvironment { nitrogen }),
            }
        }
        if self.implicit_hydrogen_cache[nitrogen] != 0 {
            return Err(NitrogenOxideCleanupError::UnsupportedHydrogens { atom: nitrogen });
        }
        for affected in [nitrogen, oxygen] {
            if self.atom_nodes[affected].chirality().is_some()
                || !self.parsed_stereo_neighbors[affected].is_empty()
            {
                return Err(NitrogenOxideCleanupError::UnsupportedStereochemistry {
                    atom: affected,
                });
            }
            if self.radical_electrons.as_ref().is_some_and(|counts| counts[affected] != 0) {
                return Err(NitrogenOxideCleanupError::StoredRadicalElectrons { atom: affected });
            }
        }
        Ok(Some(OxideEdit { nitrogen, oxygen }))
    }

    fn validate_cleanup_oxygen(&self, oxygen: usize) -> Result<(), NitrogenOxideCleanupError> {
        let atom = self.atom_nodes[oxygen];
        if atom.aromatic()
            || atom.hydrogen_count() != 0
            || self.implicit_hydrogen_cache[oxygen] != 0
            || self.bond_matrix.sparse_row_entries(oxygen).count() != 1
        {
            return Err(NitrogenOxideCleanupError::UnsupportedOxygen { oxygen });
        }
        Ok(())
    }

    fn cleanup_source_agrees(&self, source: &Self, edits: &[OxideEdit]) -> bool {
        if edits.is_empty() {
            return true;
        }
        if self.atom_nodes.len() != source.atom_nodes.len()
            || !self
                .bond_matrix
                .sparse_entries()
                .map(|(edge, _)| edge)
                .eq(source.bond_matrix.sparse_entries().map(|(edge, _)| edge))
        {
            return false;
        }
        edits.iter().all(|edit| {
            source.nitrogen_oxide_edit(edit.nitrogen) == Ok(Some(*edit))
                && [edit.nitrogen, edit.oxygen].into_iter().all(|id| {
                    self.atom_nodes[id] == source.atom_nodes[id]
                        && self.implicit_hydrogen_cache[id] == source.implicit_hydrogen_cache[id]
                        && self.parsed_stereo_neighbors[id] == source.parsed_stereo_neighbors[id]
                })
                && self.bond_matrix.sparse_row_entries(edit.nitrogen).all(|(neighbor, entry)| {
                    source.bond_entry_for_node_pair((edit.nitrogen, neighbor)).is_some_and(
                        |other| {
                            entry.descriptor() == other.descriptor()
                                && entry.order() == other.order()
                                && entry.ring_num() == other.ring_num()
                        },
                    )
                })
        })
    }

    fn apply_nitrogen_oxide_edits(&self, edits: &[OxideEdit]) -> Self {
        if edits.is_empty() {
            return self.clone_without_kekulization_source();
        }
        let mut atoms = self.atom_nodes.clone();
        for edit in edits {
            atoms[edit.nitrogen] = charged_bracket_atom(atoms[edit.nitrogen], 1);
            atoms[edit.oxygen] = charged_bracket_atom(atoms[edit.oxygen], -1);
        }
        // Validated affected inferred counts are already zero. Preserve the
        // whole cache, rather than re-inferring any unrelated atom's state.
        let bonds = BondMatrix::from_sorted_upper_triangular_entries(
            atoms.len(),
            self.bond_matrix.sparse_entries().filter(|((row, column), _)| row < column).map(
                |((row, column), entry)| {
                    let selected = [(row, column), (column, row)].into_iter().any(|(n, o)| {
                        edits
                            .binary_search_by_key(&n, |edit| edit.nitrogen)
                            .is_ok_and(|index| edits[index].oxygen == o)
                    });
                    (row, column, if selected { entry.with_bond(Bond::Single) } else { *entry })
                },
            ),
        )
        .unwrap_or_else(|_| unreachable!("existing indexed bond topology is unchanged"));
        Self::from_bond_matrix_parts_with_sidecars(
            atoms,
            bonds,
            self.parsed_stereo_neighbors.clone(),
            self.implicit_hydrogen_cache.clone(),
            self.radical_electrons.clone(),
            None,
        )
    }
}

fn first_oxygen(
    nitrogen: usize,
    candidates: &[(usize, &BondEntry)],
) -> Result<usize, NitrogenOxideCleanupError> {
    let mut non_ring =
        candidates.iter().filter(|(_, entry)| entry.ring_num().is_none()).collect::<Vec<_>>();
    non_ring.sort_unstable_by_key(|(_, entry)| entry.order());
    if non_ring.windows(2).any(|pair| pair[0].1.order() == pair[1].1.order()) {
        return Err(NitrogenOxideCleanupError::AmbiguousOxygenOrder { nitrogen });
    }
    if let Some((oxygen, _)) = non_ring.first() {
        return Ok(*oxygen);
    }
    if let [(oxygen, _)] = candidates {
        return Ok(*oxygen);
    }
    Err(NitrogenOxideCleanupError::AmbiguousOxygenOrder { nitrogen })
}

fn charged_bracket_atom(atom: Atom, charge: i8) -> Atom {
    let mut builder = Atom::builder()
        .with_symbol(atom.symbol())
        .with_hydrogens(atom.hydrogen_count())
        .with_class(atom.class());
    if let Some(isotope) = atom.isotope_mass_number() {
        builder = builder.with_isotope(isotope);
    }
    // Affected aromatic/stereo states were refused before constructing a plan.
    builder.build().with_charge_value(charge)
}

#[cfg(test)]
mod tests {
    use alloc::{boxed::Box, format};

    use geometric_traits::traits::SparseValuedMatrixRef;

    use super::{BondMatrix, NitrogenOxideCleanupError, Smiles};

    fn change_order(molecule: &mut Smiles, edge: (usize, usize), order: usize) {
        molecule.bond_matrix = BondMatrix::from_sorted_upper_triangular_entries(
            molecule.atom_nodes.len(),
            molecule.bond_matrix.sparse_entries().filter_map(|((row, column), entry)| {
                (row < column).then_some((
                    row,
                    column,
                    if (row, column) == edge { entry.with_order(order) } else { *entry },
                ))
            }),
        )
        .unwrap();
    }

    #[test]
    fn preserves_radical_absence_and_exact_present_sidecars() {
        let mut molecule: Smiles = "O=N(C)(C)C".parse().unwrap();
        assert!(molecule.radical_electrons.is_none());
        assert!(molecule.try_with_nitrogen_oxide_cleanup().unwrap().radical_electrons.is_none());
        molecule.radical_electrons = Some(vec![0, 0, 1, 0, 0]);
        assert_eq!(
            molecule.try_with_nitrogen_oxide_cleanup().unwrap().radical_electrons,
            molecule.radical_electrons
        );
        for affected in [0, 1] {
            molecule.radical_electrons.as_mut().unwrap()[affected] = 1;
            let original = format!("{molecule:?}");
            assert_eq!(
                molecule.try_with_nitrogen_oxide_cleanup().unwrap_err(),
                NitrogenOxideCleanupError::StoredRadicalElectrons { atom: affected }
            );
            assert_eq!(format!("{molecule:?}"), original);
            molecule.radical_electrons.as_mut().unwrap()[affected] = 0;
        }
    }

    #[test]
    fn rejects_inferred_hydrogen_override_without_materializing_it() {
        let mut molecule: Smiles = "O=N(C)(C)C".parse().unwrap();
        molecule.implicit_hydrogen_cache[1] = 1;
        assert_eq!(
            molecule.try_with_nitrogen_oxide_cleanup().unwrap_err(),
            NitrogenOxideCleanupError::UnsupportedHydrogens { atom: 1 }
        );
    }

    #[test]
    fn iterative_copy_preserves_nested_sources_on_match_and_no_match() {
        for input in ["O=N(C)(C)C", "C[N+](C)(C)[O-]"] {
            let mut molecule: Smiles = input.parse().unwrap();
            let mut source = molecule.clone_without_kekulization_source();
            source.radical_electrons = Some(vec![0; source.atom_nodes.len()]);
            source.kekulization_source =
                Some(Box::new(molecule.clone_without_kekulization_source()));
            molecule.kekulization_source = Some(Box::new(source));
            let original = format!("{molecule:?}");
            let cleaned = molecule.try_with_nitrogen_oxide_cleanup().unwrap();
            let mut before = Some(&molecule);
            let mut after = Some(&cleaned);
            let mut count = 0;
            while let (Some(old), Some(new)) = (before, after) {
                assert_eq!(old.radical_electrons, new.radical_electrons);
                assert_eq!(old.implicit_hydrogen_cache, new.implicit_hydrogen_cache);
                assert_eq!(old.parsed_stereo_neighbors, new.parsed_stereo_neighbors);
                assert_eq!(
                    format!(
                        "{:?}",
                        old.clone_without_kekulization_source()
                            .try_with_nitrogen_oxide_cleanup()
                            .unwrap()
                    ),
                    format!("{:?}", new.clone_without_kekulization_source())
                );
                before = old.kekulization_source.as_deref();
                after = new.kekulization_source.as_deref();
                count += 1;
            }
            assert_eq!(count, 3);
            assert!(before.is_none() && after.is_none());
            assert_eq!(format!("{molecule:?}"), original);
        }
    }

    #[test]
    fn rejects_source_charge_hydrogen_bond_and_order_conflicts_atomically() {
        for conflict in 0..4 {
            let mut molecule: Smiles = "O=N(C)(C)C".parse().unwrap();
            let mut source = molecule.clone_without_kekulization_source();
            match conflict {
                0 => source.atom_nodes[0] = source.atom_nodes[0].with_charge_value(-1),
                1 => source.implicit_hydrogen_cache[1] = 1,
                2 => source = source.try_with_nitrogen_oxide_cleanup().unwrap(),
                3 => change_order(&mut source, (0, 1), 50),
                _ => unreachable!(),
            }
            molecule.kekulization_source = Some(Box::new(source));
            let original = format!("{molecule:?}");
            assert_eq!(
                molecule.try_with_nitrogen_oxide_cleanup().unwrap_err(),
                NitrogenOxideCleanupError::SavedSourceConflict { depth: 1 }
            );
            assert_eq!(format!("{molecule:?}"), original);
        }
    }

    #[test]
    fn source_cannot_choose_a_different_oxygen_or_duplicate_competing_order() {
        let mut molecule: Smiles = "O=N(=O)C".parse().unwrap();
        let mut source = molecule.clone_without_kekulization_source();
        change_order(&mut source, (0, 1), 50);
        molecule.kekulization_source = Some(Box::new(source));
        assert_eq!(
            molecule.try_with_nitrogen_oxide_cleanup().unwrap_err(),
            NitrogenOxideCleanupError::SavedSourceConflict { depth: 1 }
        );
        molecule.kekulization_source = None;
        let order = molecule.bond_entry_for_node_pair((0, 1)).unwrap().order();
        change_order(&mut molecule, (1, 2), order);
        assert_eq!(
            molecule.try_with_nitrogen_oxide_cleanup().unwrap_err(),
            NitrogenOxideCleanupError::AmbiguousOxygenOrder { nitrogen: 1 }
        );
    }

    #[test]
    fn preserves_remote_stereo_rows_and_all_bond_metadata() {
        for input in ["N[C@@H](C)CN(=O)(C)C", "C/C=C/CN(=O)(C)C", "[H]/C(F)=C(/Cl)CN(=O)(C)C"] {
            let molecule: Smiles = input.parse().unwrap();
            let cleaned = molecule.try_with_nitrogen_oxide_cleanup().unwrap();
            assert_eq!(molecule.parsed_stereo_neighbors, cleaned.parsed_stereo_neighbors);
            for ((pair, old), (new_pair, new)) in
                molecule.bond_matrix.sparse_entries().zip(cleaned.bond_matrix.sparse_entries())
            {
                assert_eq!(pair, new_pair);
                assert_eq!(old.order(), new.order());
                assert_eq!(old.ring_num(), new.ring_num());
                assert_eq!(old.aromatic(), new.aromatic());
                if old.bond() != new.bond() {
                    assert_eq!(old.bond(), crate::bond::Bond::Double);
                    assert_eq!(new.bond(), crate::bond::Bond::Single);
                    assert_ne!(cleaned.atom_nodes[pair.0].charge_value(), 0);
                    assert_ne!(cleaned.atom_nodes[pair.1].charge_value(), 0);
                }
            }
        }
    }

    #[test]
    fn aromatic_perception_source_retains_the_same_indexed_cleanup() {
        let molecule: Smiles = "c1ccccc1N(=O)(C)C".parse().unwrap();
        let perception = molecule.perceive_aromaticity().unwrap();
        let aromatic = perception.aromaticized();
        assert!(aromatic.kekulization_source.is_some());
        let cleaned = aromatic.try_with_nitrogen_oxide_cleanup().unwrap();
        let restored = cleaned.kekulize().unwrap();
        assert_eq!(cleaned.nodes()[6].charge_value(), 1);
        assert_eq!(cleaned.nodes()[7].charge_value(), -1);
        assert_eq!(restored.nodes()[6].charge_value(), 1);
        assert_eq!(restored.nodes()[7].charge_value(), -1);
        assert_eq!(cleaned.implicit_hydrogen_cache, aromatic.implicit_hydrogen_cache);
    }
}
