//! Regression contracts for the opt-in, versioned parsed radical capability.

use smiles_rs::{
    atom::AtomSyntax,
    rdkit_2026_09_1::{
        RadicalDecision, RadicalPreparationError, RadicalPreparedSmiles, UnsupportedPreparation,
    },
};

#[test]
fn native_unpreferred_valence_rule_handles_the_observed_iron_gap() {
    // These native counts are retained prior observations, not a new oracle
    // run.
    let divalent = RadicalPreparedSmiles::try_parse("[Fe+2]").unwrap();
    assert_eq!(divalent.atom_states()[0].radical_electrons(), 0);
    assert_eq!(divalent.atom_states()[0].decision(), RadicalDecision::Assigned);
    let trivalent = RadicalPreparedSmiles::try_parse("[Fe+3].CC[O-]").unwrap();
    assert_eq!(trivalent.atom_states()[0].radical_electrons(), 1);
    assert_eq!(trivalent.atom_states().len(), trivalent.molecule().nodes().len());
}

#[test]
fn parser_origin_distinguishes_preservation_from_assignment() {
    let ordinary = RadicalPreparedSmiles::try_parse("CC").unwrap();
    for atom in ordinary.atom_states() {
        assert!(!atom.source_no_implicit());
        assert_eq!(atom.decision(), RadicalDecision::Preserved);
        assert_eq!(atom.radical_electrons(), 0);
        assert_eq!(atom.explicit_hydrogens(), 0);
        assert_eq!(atom.implicit_hydrogens(), 3);
        assert_eq!(atom.explicit_valence(), 1);
    }
    let bracket = RadicalPreparedSmiles::try_parse("[CH3]").unwrap();
    let atom = &bracket.atom_states()[0];
    assert!(atom.source_no_implicit());
    assert_eq!(atom.decision(), RadicalDecision::Assigned);
    assert_eq!(atom.radical_electrons(), 1);
    assert_eq!(atom.explicit_hydrogens(), 3);
    assert_eq!(atom.implicit_hydrogens(), 0);
}

#[test]
fn cleanup_keeps_the_original_flags_when_it_creates_bracket_atoms() {
    let prepared = RadicalPreparedSmiles::try_parse("O=N(C)(C)C").unwrap();
    assert_eq!(prepared.molecule().nodes()[0].charge_value(), -1);
    assert_eq!(prepared.molecule().nodes()[1].charge_value(), 1);
    for id in [0, 1] {
        assert_eq!(prepared.molecule().nodes()[id].syntax(), AtomSyntax::Bracket);
        assert!(!prepared.atom_states()[id].source_no_implicit());
        assert_eq!(prepared.atom_states()[id].decision(), RadicalDecision::Preserved);
    }
    let explicitly_bracketed = RadicalPreparedSmiles::try_parse("[O-][N+](C)(C)C").unwrap();
    assert!(explicitly_bracketed.atom_states()[0].source_no_implicit());
    assert!(explicitly_bracketed.atom_states()[1].source_no_implicit());
    assert_eq!(explicitly_bracketed.atom_states()[0].decision(), RadicalDecision::Assigned);
}

#[test]
fn capability_clone_keeps_indexed_metadata_and_consumption_is_explicit() {
    let prepared = RadicalPreparedSmiles::try_parse("[13CH3]").unwrap();
    let cloned = prepared.clone();
    assert_eq!(prepared.atom_states(), cloned.atom_states());
    assert_eq!(prepared.molecule().nodes(), cloned.molecule().nodes());
    let ordinary_graph = cloned.into_unqualified_smiles();
    assert_eq!(ordinary_graph.nodes()[0].isotope_mass_number(), Some(13));
    let mapped = RadicalPreparedSmiles::try_parse("[CH3:7]").unwrap();
    assert_eq!(mapped.molecule().nodes()[0].class(), 7);
}

#[test]
fn unsupported_states_remain_typed_instead_of_becoming_zero_counts() {
    assert!(matches!(
        RadicalPreparedSmiles::try_parse("c1ccccc1"),
        Err(RadicalPreparationError::UnsupportedAromaticState { .. })
    ));
    assert!(matches!(
        RadicalPreparedSmiles::try_parse("C1=CC=CC=C1"),
        Err(RadicalPreparationError::UnsupportedCycle { .. })
    ));
    assert!(matches!(
        RadicalPreparedSmiles::try_parse("[H]C"),
        Err(RadicalPreparationError::UnsupportedHydrogenVertex { .. })
    ));
    assert!(matches!(
        RadicalPreparedSmiles::try_parse("[Fe]N"),
        Err(RadicalPreparationError::UnsupportedBondedElement { .. })
    ));
    assert!(matches!(
        RadicalPreparedSmiles::try_parse("*"),
        Err(RadicalPreparationError::Parse(_))
    ));
    assert!(matches!(
        RadicalPreparedSmiles::try_parse("[Fe+16]"),
        Err(RadicalPreparationError::Parse(_))
    ));
}

#[test]
fn fixed_cleanup_guards_reject_before_any_oxide_mutation() {
    for (text, expected) in [
        ("CN=N#N", UnsupportedPreparation::NitrogenTripleNitrogen),
        ("C=P(=O)C", UnsupportedPreparation::PhosphorusOxide),
        ("O=Cl(=O)(=O)O", UnsupportedPreparation::HalogenOxide),
    ] {
        assert!(matches!(RadicalPreparedSmiles::try_parse(text),
            Err(RadicalPreparationError::UnsupportedPreparation { operation, .. })
                if operation == expected));
    }
}

#[test]
fn geometry_is_not_a_blanket_preparation_exclusion() {
    assert!(matches!(
        RadicalPreparedSmiles::try_parse("O=[N@](C)(F)Cl"),
        Err(RadicalPreparationError::NitrogenOxideCleanup(_))
    ));
    let remote = RadicalPreparedSmiles::try_parse("O=N(C)(C)C.F/C=C/F").unwrap();
    assert_eq!(remote.molecule().nodes()[0].charge_value(), -1);
    assert_eq!(remote.molecule().nodes()[1].charge_value(), 1);
}

#[test]
fn a_lexical_closure_label_does_not_imply_a_topological_cycle() {
    let prepared = RadicalPreparedSmiles::try_parse("C1.O.CC1").unwrap();
    assert!(prepared.molecule().edge_for_node_pair((0, 3)).unwrap().ring_num().is_some());
    assert_eq!(prepared.atom_states().len(), 4);
}
