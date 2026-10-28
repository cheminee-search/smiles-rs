//! Declared structural regressions, not native differential qualification.

use smiles_rs::{
    Smiles,
    rdkit_2026_09_1::{RadicalDecision, RadicalPreparationError, RadicalPreparedSmiles},
};

fn prepare(text: &str) -> Result<RadicalPreparedSmiles, RadicalPreparationError> {
    RadicalPreparedSmiles::try_parse_with_saturated_carbocycles(text)
}

fn assert_same_indexed_graph(before: &Smiles, after: &Smiles) {
    assert_eq!(before.nodes(), after.nodes());
    assert_eq!(before.implicit_hydrogen_counts(), after.implicit_hydrogen_counts());
    assert_eq!(before.number_of_bonds(), after.number_of_bonds());
    for id in 0..before.nodes().len() {
        assert_eq!(
            before.edges_for_node(id).collect::<Vec<_>>(),
            after.edges_for_node(id).collect::<Vec<_>>()
        );
    }
}

#[test]
fn named_saturated_topologies_preserve_indexed_graphs_and_checked_hydrogens() {
    // These spellings come from existing parser/ring/canonicalization fixtures.
    // The expected H vectors are declared structural expectations.
    let cases: &[(&str, &[u8])] = &[
        ("C1CC1", &[2, 2, 2]),
        ("C1CCC1", &[2, 2, 2, 2]),
        ("C1CCCCC1", &[2, 2, 2, 2, 2, 2]),
        ("CC1CCCCC1", &[3, 1, 2, 2, 2, 2, 2]),
        ("C1CCCCC1CC2CCCCC2", &[2, 2, 2, 2, 2, 1, 2, 1, 2, 2, 2, 2, 2]),
        ("C1CCC2(CC1)CCC2", &[2, 2, 2, 0, 2, 2, 2, 2, 2]),
        ("C1CC2CCC1C2", &[2, 2, 1, 2, 2, 1, 2]),
        ("C1CCC2CCCCC2C1", &[2, 2, 2, 1, 2, 2, 2, 2, 1, 2]),
        ("C12C3C4C1C5C2C3C45", &[1, 1, 1, 1, 1, 1, 1, 1]),
    ];
    for &(text, expected_h) in cases {
        let raw: Smiles = text.parse().unwrap();
        let prepared = prepare(text).unwrap();
        assert_same_indexed_graph(&raw, prepared.molecule());
        assert_eq!(prepared.atom_states().len(), expected_h.len(), "{text}");
        for (state, &hydrogens) in prepared.atom_states().iter().zip(expected_h) {
            assert_eq!(state.implicit_hydrogens(), hydrogens, "{text}");
            assert_eq!(state.explicit_hydrogens(), 0, "{text}");
            assert_eq!(state.explicit_valence(), 4 - hydrogens, "{text}");
            assert_eq!(state.radical_electrons(), 0, "{text}");
            assert_eq!(state.decision(), RadicalDecision::Preserved, "{text}");
        }
        let cloned = prepared.clone();
        assert_eq!(cloned.atom_states(), prepared.atom_states());
        assert_same_indexed_graph(prepared.molecule(), cloned.molecule());
        assert_same_indexed_graph(&raw, &cloned.into_unqualified_smiles());
        assert!(matches!(
            RadicalPreparedSmiles::try_parse(text),
            Err(RadicalPreparationError::UnsupportedCycle { .. })
        ));
    }
}

#[test]
fn bracket_hydrogens_preserve_origin_even_when_total_hydrogens_agree() {
    let implicit = prepare("C1CC1").unwrap();
    let explicit = prepare("[CH2]1[CH2][CH2]1").unwrap();
    for (left, right) in implicit.atom_states().iter().zip(explicit.atom_states()) {
        assert!(!left.source_no_implicit());
        assert!(right.source_no_implicit());
        assert_eq!(left.decision(), RadicalDecision::Preserved);
        assert_eq!(right.decision(), RadicalDecision::Assigned);
        assert_eq!(left.implicit_hydrogens(), 2);
        assert_eq!(right.explicit_hydrogens(), 2);
        assert_eq!(right.implicit_hydrogens(), 0);
        assert_eq!(right.explicit_valence(), 4);
        assert_eq!(right.radical_electrons(), 0);
    }
    let text = "[13CH2:7]1CC1";
    let labelled = prepare(text).unwrap();
    assert_eq!(labelled.molecule().nodes()[0].isotope_mass_number(), Some(13));
    assert_eq!(labelled.molecule().nodes()[0].class(), 7);
    assert_same_indexed_graph(&text.parse().unwrap(), labelled.molecule());
}

#[test]
fn cyclic_restrictions_cover_atoms_and_bonds_on_attached_tails() {
    for (text, expected_atom) in
        [("N1CC1", 0), ("OC1CC1", 0), ("C1CC1O", 3), ("[CH+]1CC1", 0), ("C1CC1[CH2+]", 3)]
    {
        assert!(
            matches!(
                prepare(text),
                Err(RadicalPreparationError::UnsupportedCyclicAtom { atom })
                    if atom == expected_atom
            ),
            "{text}"
        );
    }
    for text in ["C1=CC1", "C1CC1C=C", "C1CC1C=O", "C1C/C1", "C1=CC=CC=C1"] {
        assert!(
            matches!(prepare(text), Err(RadicalPreparationError::UnsupportedCyclicBond { .. })),
            "{text}"
        );
    }
    for (text, expected_atom) in [("[C@H]1CC1", 0), ("C1CC1[C@H](C)C", 3)] {
        assert!(
            matches!(
                prepare(text),
                Err(RadicalPreparationError::UnsupportedCyclicStereo { atom })
                    if atom == expected_atom
            ),
            "{text}"
        );
    }
}

#[test]
fn unsaturated_and_overvalent_cyclic_carbons_remain_typed_failures() {
    for (text, expected_atom, total, radicals) in
        [("[C]1CC1", 0, 2, 2), ("[CH]1CC1", 0, 3, 1), ("C1CC1[CH2]", 3, 3, 1)]
    {
        assert!(
            matches!(
                prepare(text),
                Err(RadicalPreparationError::UnsaturatedCyclicAtom {
                    atom, total_connectivity, radical_electrons,
                }) if atom == expected_atom && total_connectivity == total
                    && radical_electrons == radicals
            ),
            "{text}"
        );
    }
    for text in ["[CH3]1CC1", "C1CC1(C)(C)C"] {
        assert!(
            matches!(
                prepare(text),
                Err(RadicalPreparationError::ExplicitValenceRejected { valence: 5, .. })
            ),
            "{text}"
        );
    }
}

#[test]
fn component_flags_follow_late_unions_and_lexical_closures() {
    // The final bridge merges the already cyclic component into another root.
    // Closure labels also join lexical fragments without creating extra cycles.
    let accepted = "C1.C2C3.C1C23";
    let prepared = prepare(accepted).unwrap();
    assert_eq!(prepared.molecule().implicit_hydrogen_counts(), &[3, 2, 2, 2, 1]);
    assert_same_indexed_graph(&accepted.parse().unwrap(), prepared.molecule());
    assert!(matches!(
        prepare("C1.C2C3.O1C23"),
        Err(RadicalPreparationError::UnsupportedCyclicAtom { atom: 3 })
    ));
    let lexical_only = "C1.O.CC1";
    let old = RadicalPreparedSmiles::try_parse(lexical_only).unwrap();
    let new = prepare(lexical_only).unwrap();
    assert_eq!(old.atom_states(), new.atom_states());
    assert_same_indexed_graph(old.molecule(), new.molecule());
}

#[test]
fn disconnected_acyclic_companions_retain_existing_preparation() {
    let companions = "[Fe+3].O=N(C)(C)C.F/C=C/F";
    let acyclic = RadicalPreparedSmiles::try_parse(companions).unwrap();
    let mixed = prepare(&format!("C1CC1.{companions}")).unwrap();
    assert_eq!(&mixed.atom_states()[3..], acyclic.atom_states());
    assert_eq!(&mixed.molecule().nodes()[3..], acyclic.molecule().nodes());
    assert_eq!(mixed.atom_states()[3].radical_electrons(), 1);
    assert_eq!(mixed.molecule().nodes()[4].charge_value(), -1);
    assert_eq!(mixed.molecule().nodes()[5].charge_value(), 1);
    assert!(!mixed.atom_states()[4].source_no_implicit());
    assert!(!mixed.atom_states()[5].source_no_implicit());
    let mixed_reversed = prepare(&format!("{companions}.C1CC1")).unwrap();
    assert_eq!(&mixed_reversed.atom_states()[..acyclic.atom_states().len()], acyclic.atom_states());
}

#[test]
fn acyclic_acceptance_and_errors_are_identical_between_constructors() {
    for text in ["CC", "[CH3]", "[13CH3:7]", "[Fe+2]", "O=N(C)(C)C", "F/C=C/F"] {
        let old = RadicalPreparedSmiles::try_parse(text).unwrap();
        let new = prepare(text).unwrap();
        assert_eq!(old.atom_states(), new.atom_states(), "{text}");
        assert_same_indexed_graph(old.molecule(), new.molecule());
    }
    for text in ["*", "[H]C", "[Fe]N", "[CH5]", "CN=N#N", "C=P(=O)C", "O=Cl(=O)(=O)O"] {
        assert_eq!(
            format!("{:?}", RadicalPreparedSmiles::try_parse(text).unwrap_err()),
            format!("{:?}", prepare(text).unwrap_err()),
            "{text}"
        );
    }
    assert!(matches!(
        prepare("c1ccccc1"),
        Err(RadicalPreparationError::UnsupportedAromaticState { atom: 0 })
    ));
    assert!(matches!(
        prepare("[H]C1CC1"),
        Err(RadicalPreparationError::UnsupportedHydrogenVertex { atom: 0 })
    ));
}
