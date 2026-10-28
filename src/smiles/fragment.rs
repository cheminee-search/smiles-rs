//! Standalone fragments carved out of a parent [`Smiles`] graph.
//!
//! A [`Fragment`] keeps the parent's atom and bond typing (including aromatic
//! flags) and the parent's radical electron counts. Arbitrary atom/bond
//! fragments recompute implicit hydrogen counts for their own connectivity;
//! whole-component fragments preserve the parent's H cache, stereo ordering
//! and saved localization source. Atom ids
//! are remapped to a compact local numbering, with the parent correspondence
//! kept private behind [`Fragment::local_id`] and [`Fragment::parent_id`].

use alloc::{boxed::Box, string::String, vec::Vec};

use geometric_traits::traits::SparseValuedMatrixRef;

use super::{BondMatrixBuilder, ConcreteAtoms, Smiles, SmilesAtomPolicy};
use crate::{
    bond::bond_edge::BondEdge,
    errors::{RootError, SubgraphError},
};

/// A standalone molecule carved out of a parent graph, with a private mapping
/// back to the parent's atom ids.
#[derive(Debug, Clone)]
pub struct Fragment<AtomPolicy = ConcreteAtoms> {
    smiles: Smiles<AtomPolicy>,
    parent_of_local: Vec<usize>,
    local_of_parent: Vec<usize>,
}

impl<AtomPolicy: SmilesAtomPolicy> Fragment<AtomPolicy> {
    /// Returns the fragment as a standalone [`Smiles`].
    #[inline]
    #[must_use]
    pub fn smiles(&self) -> &Smiles<AtomPolicy> {
        &self.smiles
    }

    /// Consumes the fragment and returns the standalone [`Smiles`].
    #[inline]
    #[must_use]
    pub fn into_smiles(self) -> Smiles<AtomPolicy> {
        self.smiles
    }

    /// Returns the number of atoms in the fragment.
    #[inline]
    #[must_use]
    pub fn atom_count(&self) -> usize {
        self.parent_of_local.len()
    }

    /// Returns the fragment-local id of a parent atom, if it is in the
    /// fragment.
    #[inline]
    #[must_use]
    pub fn local_id(&self, parent_atom: usize) -> Option<usize> {
        self.local_of_parent.get(parent_atom).copied().filter(|&local| local != usize::MAX)
    }

    /// Returns the parent id of a fragment atom.
    ///
    /// # Panics
    ///
    /// Panics if `local_atom` is not a valid fragment atom id.
    #[inline]
    #[must_use]
    pub fn parent_id(&self, local_atom: usize) -> usize {
        self.parent_of_local[local_atom]
    }

    /// Renders the fragment anchored at a parent atom. When `isomeric` is
    /// false, isotopes and stereochemistry are dropped (RDKit
    /// `isomericSmiles=False`).
    ///
    /// # Errors
    ///
    /// Returns [`RootError::AtomNotInFragment`] if `parent_atom` is not part of
    /// the fragment.
    pub fn render_rooted(&self, parent_atom: usize, isomeric: bool) -> Result<String, RootError> {
        let local = self.local_id(parent_atom).ok_or(RootError::AtomNotInFragment(parent_atom))?;
        if isomeric {
            Ok(self.smiles.render_rooted(local))
        } else {
            Ok(self.smiles.non_isomeric().render_rooted(local))
        }
    }
}

impl<AtomPolicy: SmilesAtomPolicy> Smiles<AtomPolicy> {
    /// Extracts one complete connected component without changing its atom,
    /// bond, hydrogen, radical or stereo state.
    ///
    /// Component identifiers use [`Self::connected_components`] ordering. The
    /// fragment retains parent/local atom mappings and projects any saved
    /// localization source onto the same parent atom IDs. No sanitization,
    /// canonicalization or hydrogen re-inference is performed.
    ///
    /// # Errors
    /// Returns [`SubgraphError::ComponentOutOfRange`] if `component_id` is not
    /// a connected component of this graph.
    ///
    /// # Examples
    /// ```
    /// use smiles_rs::Smiles;
    /// let molecule: Smiles = "CC.O".parse().expect("valid SMILES");
    /// let fragment = molecule.fragment_from_component(0).expect("first component");
    /// assert_eq!(fragment.atom_count(), 2);
    /// assert_eq!(fragment.parent_id(0), 0);
    /// assert_eq!(fragment.smiles().to_string(), "CC");
    /// ```
    pub fn fragment_from_component(
        &self,
        component_id: usize,
    ) -> Result<Fragment<AtomPolicy>, SubgraphError> {
        let components = self.connected_components();
        if component_id >= components.number_of_components() {
            return Err(SubgraphError::ComponentOutOfRange(component_id));
        }
        let parent_of_local = components.node_ids_of_component(component_id).collect::<Vec<_>>();
        let mut local_of_parent = vec![usize::MAX; self.nodes().len()];
        for (local, &parent) in parent_of_local.iter().enumerate() {
            local_of_parent[parent] = local;
        }

        let mut smiles = self.induced_component_with_mapping(&parent_of_local, &local_of_parent);
        let mut output = &mut smiles;
        let mut source = self.kekulization_source.as_deref();
        while let Some(saved) = source {
            // Aromaticity assignment and the operations preserving its source
            // retain node IDs and connectivity. Orders, aromatic labels and H
            // spelling can differ, so project the saved graph's own state.
            debug_assert_eq!(saved.nodes().len(), self.nodes().len());
            debug_assert!(
                saved
                    .bond_matrix()
                    .sparse_entries()
                    .map(|(edge, _)| edge)
                    .eq(self.bond_matrix().sparse_entries().map(|(edge, _)| edge))
            );
            output.kekulization_source = Some(Box::new(
                saved.induced_component_with_mapping(&parent_of_local, &local_of_parent),
            ));
            output = output
                .kekulization_source
                .as_deref_mut()
                .unwrap_or_else(|| unreachable!("source was just inserted"));
            source = saved.kekulization_source.as_deref();
        }

        Ok(Fragment { smiles, parent_of_local, local_of_parent })
    }

    /// Builds a fragment from an atom set: every bond among the chosen atoms is
    /// included (an induced subgraph). Implicit hydrogen counts are recomputed
    /// for the fragment, while radical electron counts are kept from the
    /// parent.
    ///
    /// # Errors
    ///
    /// Returns [`SubgraphError::AtomOutOfRange`] if any id is not a valid atom.
    pub fn fragment_from_atoms(
        &self,
        atoms: impl IntoIterator<Item = usize>,
    ) -> Result<Fragment<AtomPolicy>, SubgraphError> {
        let node_count = self.nodes().len();
        let mut parent_of_local = Vec::new();
        let mut local_of_parent = vec![usize::MAX; node_count];
        for atom in atoms {
            if atom >= node_count {
                return Err(SubgraphError::AtomOutOfRange(atom));
            }
            register_local(&mut parent_of_local, &mut local_of_parent, atom);
        }

        let mut builder = BondMatrixBuilder::default();
        for ((row, column), entry) in self.bond_matrix.sparse_entries() {
            if row >= column {
                continue;
            }
            let (local_row, local_column) = (local_of_parent[row], local_of_parent[column]);
            if local_row == usize::MAX || local_column == usize::MAX {
                continue;
            }
            builder
                .push_edge_with_descriptor(local_row, local_column, entry.descriptor(), None)
                .unwrap_or_else(|_| unreachable!("induced fragment preserves a simple graph"));
        }

        Ok(self.finish_fragment(parent_of_local, local_of_parent, builder))
    }

    /// Builds a fragment from a bond set (RDKit `PathToSubmol`): exactly the
    /// listed bonds are kept and their endpoints become the fragment's atoms.
    /// Implicit hydrogen counts are recomputed for the fragment, while radical
    /// electron counts are kept from the parent.
    ///
    /// # Errors
    ///
    /// Returns [`SubgraphError::AtomOutOfRange`] if an endpoint is not a valid
    /// atom, or [`SubgraphError::BondReferencesUnknownAtom`] if a listed bond
    /// is not an edge of the parent graph.
    pub fn fragment_from_bonds(
        &self,
        bonds: impl IntoIterator<Item = BondEdge>,
    ) -> Result<Fragment<AtomPolicy>, SubgraphError> {
        let node_count = self.nodes().len();
        let mut parent_of_local = Vec::new();
        let mut local_of_parent = vec![usize::MAX; node_count];
        let mut builder = BondMatrixBuilder::default();
        for edge in bonds {
            let [source, target] = edge.endpoints();
            if source >= node_count {
                return Err(SubgraphError::AtomOutOfRange(source));
            }
            if target >= node_count {
                return Err(SubgraphError::AtomOutOfRange(target));
            }
            let descriptor = self
                .edge_for_node_pair((source, target))
                .ok_or(SubgraphError::BondReferencesUnknownAtom(target))?
                .descriptor();
            let local_source = register_local(&mut parent_of_local, &mut local_of_parent, source);
            let local_target = register_local(&mut parent_of_local, &mut local_of_parent, target);
            builder
                .push_edge_with_descriptor(local_source, local_target, descriptor, None)
                .map_err(|_| SubgraphError::BondReferencesUnknownAtom(target))?;
        }

        Ok(self.finish_fragment(parent_of_local, local_of_parent, builder))
    }

    fn finish_fragment(
        &self,
        parent_of_local: Vec<usize>,
        local_of_parent: Vec<usize>,
        builder: BondMatrixBuilder,
    ) -> Fragment<AtomPolicy> {
        let atom_nodes: Vec<_> =
            parent_of_local.iter().map(|&parent| self.atom_nodes[parent]).collect();
        let atom_count = atom_nodes.len();
        let parsed_stereo_neighbors = vec![Vec::new(); atom_count];
        let parent_radical_electrons = self.radical_electron_counts();
        let radical_electrons = Some(
            parent_of_local
                .iter()
                .copied()
                .map(|parent| parent_radical_electrons[parent])
                .collect(),
        );
        let smiles = Self::from_bond_matrix_parts_with_parsed_stereo_and_radicals(
            atom_nodes,
            builder.finish(atom_count),
            parsed_stereo_neighbors,
            radical_electrons,
        );
        Fragment { smiles, parent_of_local, local_of_parent }
    }
}

/// Assigns (or returns) the local id of a parent atom, recording it on first
/// sight.
fn register_local(
    parent_of_local: &mut Vec<usize>,
    local_of_parent: &mut [usize],
    parent_atom: usize,
) -> usize {
    let local = &mut local_of_parent[parent_atom];
    if *local == usize::MAX {
        *local = parent_of_local.len();
        parent_of_local.push(parent_atom);
    }
    *local
}

#[cfg(test)]
mod tests {
    use alloc::{boxed::Box, vec::Vec};

    use super::{Fragment, Smiles};
    use crate::{SubgraphError, smiles::StereoNeighbor};

    fn assert_mapped_state(parent: &Smiles, selected: &Smiles, fragment: &Fragment) {
        assert_eq!(selected.nodes().len(), fragment.atom_count());
        let parent_radicals = parent.radical_electron_counts();
        let selected_radicals = selected.radical_electron_counts();
        for local in 0..fragment.atom_count() {
            let original = fragment.parent_id(local);
            assert_eq!(fragment.local_id(original), Some(local));
            assert_eq!(selected.nodes()[local], parent.nodes()[original]);
            assert_eq!(
                selected.implicit_hydrogen_count(local),
                parent.implicit_hydrogen_count(original)
            );
            assert_eq!(selected_radicals[local], parent_radicals[original]);
            let expected_neighbors = parent
                .parsed_stereo_neighbors_row(original)
                .iter()
                .map(|neighbor| {
                    match *neighbor {
                        StereoNeighbor::Atom(atom) => {
                            StereoNeighbor::Atom(fragment.local_id(atom).unwrap())
                        }
                        StereoNeighbor::ExplicitHydrogen => StereoNeighbor::ExplicitHydrogen,
                    }
                })
                .collect::<Vec<_>>();
            assert_eq!(selected.parsed_stereo_neighbors_row(local), expected_neighbors);
            assert_eq!(
                selected.edges_for_node(local).count(),
                parent.edges_for_node(original).count()
            );
            for edge in parent.edges_for_node(original) {
                let neighbor = fragment.local_id(edge.other(original).unwrap()).unwrap();
                assert_eq!(
                    selected.edge_for_node_pair((local, neighbor)).unwrap().descriptor(),
                    edge.descriptor()
                );
            }
        }
    }

    fn assert_component(parent: &Smiles, id: usize) {
        let before = format!("{parent:?}");
        let components = parent.connected_components();
        let expected_atoms = components.node_ids_of_component(id).collect::<Vec<_>>();
        let fragment = parent.fragment_from_component(id).unwrap();
        assert_eq!(
            (0..fragment.atom_count()).map(|i| fragment.parent_id(i)).collect::<Vec<_>>(),
            expected_atoms
        );
        for atom in 0..parent.nodes().len() {
            assert_eq!(fragment.local_id(atom).is_some(), components.component_of_node(atom) == id);
        }
        assert_eq!(fragment.local_id(parent.nodes().len()), None);
        let mut original = Some(parent);
        let mut selected = Some(fragment.smiles());
        while let Some(source) = original {
            let graph = selected.expect("every saved source must be retained");
            assert_mapped_state(source, graph, &fragment);
            original = source.kekulization_source.as_deref();
            selected = graph.kekulization_source.as_deref();
        }
        assert!(selected.is_none());
        assert_eq!(format!("{parent:?}"), before);
    }

    #[test]
    fn whole_component_keeps_parsed_stereo_rows() {
        let smiles: Smiles = "[Na+].N[C@@H](C)C(=O)O.[Cl-]".parse().unwrap();
        let components = smiles.connected_components();
        let component = components.component_of_node(1);
        assert_component(&smiles, component);
    }

    #[test]
    fn whole_components_preserve_ring_stereo_h_and_radicals() {
        for input in [
            "[Na+].F[C@]1(Cl)CCC(C)C1",
            "N[C@@H]1CC[C@H](O)C1.[Na+]",
            "[H][C@]1(F)CCC(Cl)C1.[Na+]",
            "N[C@H]1.F.C(O)CC1",
            r"Br/C([H])=C(\[H])F.[Na+]",
            r"F/C([H])=C(/[H])C([H])=C([H])C(/[H])=C(\[H])Cl.[Na+]",
            "C[S@](=O)CC.[Na+]",
            "[CH3].O",
            "[13CH4].C",
            "[H:7]C.C",
        ] {
            let parent: Smiles = input.parse().unwrap();
            for component in 0..parent.connected_components().number_of_components() {
                assert_component(&parent, component);
            }
        }
    }

    #[test]
    fn whole_components_project_saved_localization_sources() {
        for input in ["C1=CC=CC=C1.CC", "[13CH]1=CC=CC=C1.[Na+]"] {
            let original: Smiles = input.parse().unwrap();
            let assignment = original.aromaticity_assignment();
            let aromatic = original.try_with_aromaticity_assignment(&assignment).unwrap();
            assert!(aromatic.kekulization_source.is_some());
            for id in 0..aromatic.connected_components().number_of_components() {
                assert_component(&aromatic, id);
                let fragment = aromatic.fragment_from_component(id).unwrap();
                let restored = fragment.smiles().kekulize().unwrap();
                assert_mapped_state(&original, &restored, &fragment);
                let reperceived = restored.aromaticity_assignment();
                let expected = fragment.smiles().aromaticity_assignment();
                assert_eq!(reperceived.status(), expected.status());
                assert_eq!(reperceived.atom_ids(), expected.atom_ids());
                assert_eq!(reperceived.bond_edges(), expected.bond_edges());
            }
        }
    }

    #[test]
    fn whole_component_projects_each_saved_source_with_the_same_mapping() {
        let original: Smiles = "[13CH]1=CC=CC=C1.[Na+]".parse().unwrap();
        let assignment = original.aromaticity_assignment();
        let mut parent = original.try_with_aromaticity_assignment(&assignment).unwrap();
        // Public aromaticity operations normally keep one source. Exercise an
        // internal multi-link chain to ensure projection never truncates it.
        parent.kekulization_source.as_deref_mut().unwrap().kekulization_source =
            Some(Box::new(original.clone()));
        for id in 0..parent.connected_components().number_of_components() {
            assert_component(&parent, id);
        }
    }

    #[test]
    fn whole_component_rejects_invalid_identifiers() {
        let parent: Smiles = "CCO".parse().unwrap();
        for id in [1, usize::MAX] {
            assert_eq!(
                parent.fragment_from_component(id).unwrap_err(),
                SubgraphError::ComponentOutOfRange(id)
            );
        }
        let empty = parent.fragment_from_atoms([]).unwrap().into_smiles();
        assert_eq!(
            empty.fragment_from_component(0).unwrap_err(),
            SubgraphError::ComponentOutOfRange(0)
        );
    }

    /// Regression test for ClusterFuzzLite run 38211063533: an environment
    /// fragment with two chemically-equivalent explicit `[H]` atoms (a
    /// genuine graph automorphism) rendered a different ring-closure digit
    /// depending only on which order its atoms were handed to
    /// `fragment_from_atoms`.
    #[test]
    fn render_rooted_is_independent_of_fragment_atom_order() {
        let smiles: Smiles = "c1[H]N0cc1c[H]0nc1c[H]0N#1c[H]0nccc".parse().unwrap();
        let center = 0;
        let radius = 3;
        let environment = smiles.atom_environment(center, radius).unwrap();
        let atoms: Vec<usize> = environment.atoms().collect();
        let reversed: Vec<usize> = atoms.iter().rev().copied().collect();

        let forward =
            smiles.fragment_from_atoms(atoms).unwrap().render_rooted(center, false).unwrap();
        let backward =
            smiles.fragment_from_atoms(reversed).unwrap().render_rooted(center, false).unwrap();

        assert_eq!(forward, backward);
    }
}
