//! Atom statistics of a concrete graph and its recorded hydrogen counts.

use elements_rs::Element;

use super::Smiles;

/// An exact atom count could not be represented as a `u64`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
#[error("atom statistics exceed u64 count capacity")]
pub struct AtomCountOverflow;

/// Immutable counts from a concrete graph and its recorded hydrogen metadata.
///
/// These statistics do not sanitize, normalize, or prepare a molecule. Inferred
/// hydrogens retain the raw graph policy of
/// [`Smiles::implicit_hydrogen_counts`]. Graph hydrogen atoms, bracket
/// hydrogens, and inferred hydrogens are counted separately. Isotope, charge,
/// and radical annotations are not reinterpreted.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct AtomStatistics {
    graph_atoms: u64,
    bracket_hydrogens: u64,
    inferred_hydrogens: u64,
    atoms_including_hydrogens: u64,
    heavy_atoms: u64,
    carbon_atoms: u64,
    carbon_total_connectivity_four: u64,
    nitrogen_oxygen_atoms: u64,
    hydrogens_on_nitrogen_oxygen: u64,
}

impl AtomStatistics {
    /// Number of graph vertices, including graph hydrogen atoms.
    #[must_use]
    pub const fn graph_atoms(&self) -> u64 {
        self.graph_atoms
    }

    /// Hydrogens recorded within atom brackets, without separate graph
    /// vertices.
    #[must_use]
    pub const fn bracket_hydrogens(&self) -> u64 {
        self.bracket_hydrogens
    }

    /// Sum of the graph's recorded inferred hydrogen counts.
    #[must_use]
    pub const fn inferred_hydrogens(&self) -> u64 {
        self.inferred_hydrogens
    }

    /// Graph vertices plus bracket and inferred hydrogens.
    ///
    /// Graph hydrogen atoms are already vertices and are not added again.
    #[must_use]
    pub const fn atoms_including_hydrogens(&self) -> u64 {
        self.atoms_including_hydrogens
    }

    /// Graph atoms with atomic number greater than one, excluding every H
    /// isotope.
    #[must_use]
    pub const fn heavy_atoms(&self) -> u64 {
        self.heavy_atoms
    }

    /// Carbon graph atoms, independent of isotope, charge, or radical
    /// annotation.
    #[must_use]
    pub const fn carbon_atoms(&self) -> u64 {
        self.carbon_atoms
    }

    /// Carbon atoms whose graph degree plus bracket and inferred H equals four.
    ///
    /// This tests total connectivity, not a perceived hybridization state.
    #[must_use]
    pub const fn carbon_total_connectivity_four(&self) -> u64 {
        self.carbon_total_connectivity_four
    }

    /// Fraction of carbon atoms with total connectivity four, or zero if none.
    ///
    /// Computed by ordinary `f64` conversion and division of the exact counts.
    /// Counts above the exact integer range of `f64` can lose precision here;
    /// the integer accessors retain their exact values.
    #[must_use]
    #[allow(clippy::cast_precision_loss)]
    pub fn carbon_total_connectivity_four_fraction(&self) -> f64 {
        if self.carbon_atoms == 0 {
            0.0
        } else {
            self.carbon_total_connectivity_four as f64 / self.carbon_atoms as f64
        }
    }

    /// Nitrogen and oxygen graph atoms, regardless of donor or acceptor
    /// chemistry.
    #[must_use]
    pub const fn nitrogen_oxygen_atoms(&self) -> u64 {
        self.nitrogen_oxygen_atoms
    }

    /// Bracket, inferred, and adjacent graph hydrogens on nitrogen or oxygen.
    ///
    /// Isotopic and charged graph H count as hydrogen. This counts attachment
    /// incidences: a graph H adjacent to two N/O atoms contributes twice.
    #[must_use]
    pub const fn hydrogens_on_nitrogen_oxygen(&self) -> u64 {
        self.hydrogens_on_nitrogen_oxygen
    }
}

impl Smiles {
    /// Counts atoms and recorded hydrogens without changing the supplied graph.
    ///
    /// This uses the raw recorded hydrogen policy of
    /// [`Self::implicit_hydrogen_counts`]. It does not sanitize, remove graph
    /// H, perceive hybridization, normalize charge/radicals, select fragments,
    /// or perform a preparation round trip. The result describes this graph
    /// state, not a native descriptor or a chemically validated molecule.
    ///
    /// The calculation borrows existing graph and hydrogen storage and needs no
    /// molecule clone or intermediate collection. It is available for concrete
    /// [`Smiles`]; wildcard and SMARTS query graphs are outside this API.
    ///
    /// # Errors
    /// Returns [`AtomCountOverflow`] if an exact count exceeds `u64` capacity.
    ///
    /// # Examples
    ///
    /// ```
    /// use smiles_rs::Smiles;
    /// let molecule: Smiles = "CCO".parse().expect("valid SMILES");
    /// let counts = molecule.atom_statistics()?;
    /// assert_eq!(counts.graph_atoms(), 3);
    /// assert_eq!(counts.inferred_hydrogens(), 6);
    /// assert_eq!(counts.atoms_including_hydrogens(), 9);
    /// assert_eq!(counts.heavy_atoms(), 3);
    /// assert_eq!(counts.hydrogens_on_nitrogen_oxygen(), 1);
    /// assert_eq!(counts.carbon_total_connectivity_four_fraction(), 1.0);
    /// # Ok::<(), smiles_rs::AtomCountOverflow>(())
    /// ```
    pub fn atom_statistics(&self) -> Result<AtomStatistics, AtomCountOverflow> {
        let mut result = AtomStatistics {
            graph_atoms: count_u64(self.nodes().len())?,
            ..AtomStatistics::default()
        };
        for (atom_id, atom) in self.nodes().iter().enumerate() {
            let element = atom.element().expect("concrete Smiles atoms have an element");
            let bracket_hydrogens = u64::from(atom.hydrogen_count());
            let inferred_hydrogens = u64::from(self.implicit_hydrogen_count(atom_id));
            let non_graph_hydrogens = bracket_hydrogens + inferred_hydrogens;
            add_count(&mut result.bracket_hydrogens, bracket_hydrogens)?;
            add_count(&mut result.inferred_hydrogens, inferred_hydrogens)?;
            if element != Element::H {
                add_count(&mut result.heavy_atoms, 1)?;
            }
            if element == Element::C {
                add_count(&mut result.carbon_atoms, 1)?;
                let connectivity = count_u64(self.edge_count_for_node(atom_id))?
                    .checked_add(non_graph_hydrogens)
                    .ok_or(AtomCountOverflow)?;
                if connectivity == 4 {
                    add_count(&mut result.carbon_total_connectivity_four, 1)?;
                }
            }
            if matches!(element, Element::N | Element::O) {
                add_count(&mut result.nitrogen_oxygen_atoms, 1)?;
                add_count(&mut result.hydrogens_on_nitrogen_oxygen, non_graph_hydrogens)?;
                let graph_hydrogens = self
                    .edges_for_node(atom_id)
                    .filter_map(|edge| edge.other(atom_id))
                    .filter(|&neighbor| self.nodes()[neighbor].element() == Some(Element::H))
                    .count();
                add_count(&mut result.hydrogens_on_nitrogen_oxygen, count_u64(graph_hydrogens)?)?;
            }
        }
        result.atoms_including_hydrogens = result
            .graph_atoms
            .checked_add(result.bracket_hydrogens)
            .and_then(|count| count.checked_add(result.inferred_hydrogens))
            .ok_or(AtomCountOverflow)?;
        Ok(result)
    }
}

fn count_u64(count: usize) -> Result<u64, AtomCountOverflow> {
    u64::try_from(count).map_err(|_| AtomCountOverflow)
}

fn add_count(count: &mut u64, addition: u64) -> Result<(), AtomCountOverflow> {
    *count = count.checked_add(addition).ok_or(AtomCountOverflow)?;
    Ok(())
}
