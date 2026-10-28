//! Structural admission regressions; these are not a native differential panel.

use smiles_rs::{
    Smiles, TetrahedralNeighbor,
    rdkit_2026_09_1::{RadicalPreparationError, RadicalPreparedSmiles},
};

fn prepare(text: &str) -> Result<RadicalPreparedSmiles, RadicalPreparationError> {
    RadicalPreparedSmiles::try_parse_with_tetrahedral_carbocycles(text)
}

#[test]
fn checked_centers_preserve_source_slots_graph_and_native_h_state() {
    for text in [
        "C[C@H]1CC(C)CCC1",
        "C[C@@H]1CC(C)CCC1",
        "C[C@TH1H]1CC(C)CCC1",
        "C[C@TH2H]1CC(C)CCC1",
        "C[C@@H](CC(C)CCC1)1",
        "[C@H]1(C)CC(C)CCC1",
        "C[C@]1(CC)CC(C)CCC1",
        "C1CC1[C@H](C)CC",
        "C[C@H]1C[13CH2]CCC1",
    ] {
        let raw: Smiles = text.parse().unwrap();
        let prepared = prepare(text).unwrap();
        assert_eq!(raw.nodes(), prepared.molecule().nodes(), "{text}");
        assert_eq!(raw.implicit_hydrogen_counts(), prepared.molecule().implicit_hydrogen_counts());
        for (atom, state) in prepared.atom_states().iter().enumerate() {
            assert_eq!(state.radical_electrons(), 0, "{text}");
            assert_eq!(
                raw.edges_for_node(atom).collect::<Vec<_>>(),
                prepared.molecule().edges_for_node(atom).collect::<Vec<_>>()
            );
            assert_eq!(
                prepared.checked_tetrahedral_neighbors(atom),
                raw.tetrahedral_stereo_neighbors(atom)
            );
        }
        let clone = prepared.clone();
        drop(prepared);
        for atom in 0..raw.nodes().len() {
            if let Some(slots) = clone.checked_tetrahedral_neighbors(atom) {
                assert_eq!(
                    raw.tetrahedral_slots_anticlockwise(atom, slots),
                    clone.molecule().tetrahedral_slots_anticlockwise(atom, slots)
                );
                assert!(!slots.contains(&TetrahedralNeighbor::LonePair));
            }
        }
        assert_eq!(clone.checked_tetrahedral_neighbors(raw.nodes().len()), None);
    }
}

#[test]
fn map_and_source_h_spelling_cannot_resolve_equivalent_ligands() {
    for text in [
        "C[C@H]1CCCCC1",
        "C[C@H]1[CH2:7]CCCC1",
        "C[C@H]1[CH2]CCCC1",
        "C[C@H]1CCC(C)CC1",
        "C1CC1[C@H](C)C",
        "[C@H]12C3C4C1C5C2C3C45",
    ] {
        assert!(
            matches!(
                prepare(text),
                Err(RadicalPreparationError::UnresolvedTetrahedralLigands { .. })
            ),
            "{text}"
        );
    }
}

#[test]
fn unsupported_slots_classes_and_acyclic_tags_remain_errors() {
    for text in [
        "[C@H]1CC1",
        "C[C@SP1H]1CC(C)CCC1",
        "C[C@TB1H]1CC(C)CCC1",
        "C[C@OH1H]1CC(C)CCC1",
        "C1CC1.F[C@H](Cl)Br",
    ] {
        assert!(
            matches!(
                prepare(text),
                Err(RadicalPreparationError::UnsupportedTetrahedralStereo { .. })
            ),
            "{text}"
        );
    }
}

#[test]
fn original_constructors_do_not_gain_stereo_certificates() {
    let acyclic = RadicalPreparedSmiles::try_parse("F[C@H](Cl)Br").unwrap();
    assert_eq!(acyclic.checked_tetrahedral_neighbors(1), None);
    assert!(matches!(
        RadicalPreparedSmiles::try_parse_with_saturated_carbocycles("C[C@H]1CC(C)CCC1"),
        Err(RadicalPreparationError::UnsupportedCyclicStereo { atom: 1 })
    ));
    assert!(matches!(
        RadicalPreparedSmiles::try_parse("C[C@H]1CC(C)CCC1"),
        Err(RadicalPreparationError::UnsupportedCycle { .. })
    ));
}

#[test]
fn cyclic_whole_component_exclusions_are_preserved() {
    for text in
        ["C1CC1O", "C1CC1[CH2+]", "C1=CC1", "C1CC1[CH2]", "[CH3]1CC1", "[H]C1CC1", "c1ccccc1"]
    {
        assert_eq!(
            format!("{:?}", prepare(text).unwrap_err()),
            format!(
                "{:?}",
                RadicalPreparedSmiles::try_parse_with_saturated_carbocycles(text).unwrap_err()
            ),
            "{text}"
        );
    }
}

#[test]
fn disconnected_nonstereo_radicals_keep_their_sealed_state() {
    let prepared = prepare("[Fe+3].[CH3].C[C@H]1CC(C)CCC1").unwrap();
    assert_eq!(prepared.atom_states()[0].radical_electrons(), 1);
    assert_eq!(prepared.atom_states()[1].radical_electrons(), 1);
    assert!(prepared.checked_tetrahedral_neighbors(3).is_some());
}
