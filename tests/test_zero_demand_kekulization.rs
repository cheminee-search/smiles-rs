//! Standalone localization of aromatic components with no new pi-bond demand.

use smiles_rs::{AromaticityStatus, KekulizationError, Smiles, bond::Bond};

fn assert_all_single_ring_localization(source: &str) {
    let molecule = source.parse::<Smiles>().expect("reference SMILES should parse");
    let localized = molecule
        .kekulize_standalone()
        .unwrap_or_else(|error| panic!("zero-demand component {source} did not localize: {error}"));

    assert_eq!(localized.nodes().len(), molecule.nodes().len());
    assert_eq!(localized.number_of_bonds(), molecule.number_of_bonds());
    assert_eq!(localized.implicit_hydrogen_counts(), molecule.implicit_hydrogen_counts());
    for (original, actual) in molecule.nodes().iter().zip(localized.nodes()) {
        assert!(!actual.aromatic());
        assert_eq!(original.symbol(), actual.symbol());
        assert_eq!(original.isotope_mass_number(), actual.isotope_mass_number());
        assert_eq!(original.charge(), actual.charge());
        assert_eq!(original.hydrogen_count(), actual.hydrogen_count());
        assert_eq!(original.chirality(), actual.chirality());
        assert_eq!(original.class(), actual.class());
    }
    for atom_id in 0..molecule.nodes().len() {
        for edge in molecule.edges_for_node(atom_id) {
            let actual = localized.edge_for_node_pair((edge.source(), edge.target())).unwrap();
            assert!(!actual.is_aromatic());
            assert_eq!(actual.bond(), if edge.is_aromatic() { Bond::Single } else { edge.bond() });
        }
    }

    let reperceived = localized.perceive_aromaticity().expect("reference aromaticity is supported");
    assert_eq!(reperceived.status(), AromaticityStatus::Complete);
    for (atom_id, atom) in molecule.nodes().iter().enumerate() {
        assert_eq!(reperceived.assignment().contains_atom(atom_id), atom.aromatic());
        for edge in molecule.edges_for_node(atom_id) {
            assert_eq!(
                reperceived.assignment().contains_edge(edge.source(), edge.target()),
                edge.is_aromatic()
            );
        }
    }
}

#[test]
fn standalone_localizes_imino_thione_ring_without_new_double_bonds() {
    assert_all_single_ring_localization("N=c1[nH]c(=S)ss1");
}

#[test]
fn standalone_localizes_dioxo_ring_without_new_double_bonds() {
    assert_all_single_ring_localization("O=c1[nH]c(=O)ss1");
}

#[test]
fn empty_matching_supports_multiple_lone_pair_donor_elements() {
    for source in ["O=c1[nH]c(=O)os1", "O=c1[nH]c(=O)oo1", "s1ssss1", "o1oooo1", "[nH]1[nH][nH]1"] {
        assert_all_single_ring_localization(source);
    }
}

#[test]
fn empty_matching_preserves_isotopes_maps_hydrogens_charge_and_stereo() {
    for source in [
        "[18O]=c1[nH]c(=[34S])ss1",
        "[O:1]=[c:2]1[nH:3][c:4](=[O:5])[s:6][s:7]1",
        "O=c1n([H])c(=O)ss1",
        "O=c1n([2H])c(=O)ss1",
        "O=c1[n-]c(=O)ss1.[Na+]",
        "O=c1n(C[C@H](F)Cl)c(=O)ss1",
        "O=c1n(C[C@@H](F)Cl)c(=O)ss1",
        "F/C=C/F.O=c1[nH]c(=O)ss1",
        "F/C=C\\F.O=c1[nH]c(=O)ss1",
    ] {
        assert_all_single_ring_localization(source);
    }
}

#[test]
fn zero_demand_components_can_coexist_with_ordinary_matching_components() {
    for source in ["O=c1[nH]c(=O)ss1.c1ccccc1", "c1ccccc1.O=c1[nH]c(=O)ss1"] {
        let molecule = source.parse::<Smiles>().unwrap();
        let localized = molecule.kekulize_standalone().expect("each component should localize");
        assert_eq!(localized.implicit_hydrogen_counts(), molecule.implicit_hydrogen_counts());
        let reperceived = localized.perceive_aromaticity().unwrap();
        assert_eq!(reperceived.assignment(), &molecule.aromaticity_assignment());
    }
}

#[test]
fn empty_matching_is_independent_of_input_atom_order() {
    for source in ["N=c1[nH]c(=S)ss1", "O=c1[nH]c(=O)ss1"] {
        let molecule = source.parse::<Smiles>().unwrap();
        for root in 0..molecule.nodes().len() {
            assert_all_single_ring_localization(&molecule.render_rooted(root));
        }
    }
}

#[test]
fn empty_matching_rejects_inconsistent_aromatic_support() {
    for source in [
        "s1sss1",
        "o1ooo1",
        "[nH]1[nH][nH][nH]1",
        "C=c1[nH]c(=O)ss1",
        "O=c1[nH2+]c(=O)ss1",
        // RDKit accepts these directional spellings, but the parser does not
        // flag their directional ring bond as aromatic. Preserve the explicit
        // unsupported result until that separate representation gap is fixed.
        "C/N=c1/[nH]c(=O)ss1",
        "C/N=c1\\[nH]c(=O)ss1",
    ] {
        assert_eq!(
            source.parse::<Smiles>().unwrap().kekulize_standalone(),
            Err(KekulizationError::NoPerfectMatching { candidate_atom_count: 0 }),
            "unsupported zero-demand component {source}"
        );
    }
    for source in ["c1cccc1", "n1cccc1", "O=c1[nH]c(=O)ss1.c1cccc1"] {
        assert!(matches!(
            source.parse::<Smiles>().unwrap().kekulize_standalone(),
            Err(KekulizationError::NoPerfectMatching { .. })
        ));
    }
}
