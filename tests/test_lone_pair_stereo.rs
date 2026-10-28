//! Tetrahedral sulfur with a distinct virtual lone-pair neighbor.

use TetrahedralNeighbor::{Atom, BracketHydrogen, LonePair};
use serde::Deserialize;
use smiles_rs::{
    Smiles, TetrahedralNeighbor, WildcardSmiles, atom::bracketed::chirality::Chirality,
};

#[test]
fn sulfoxide_exposes_semantic_tetrahedral_chirality() {
    let molecule: Smiles = "CC[S@](=O)CCC".parse().expect("valid sulfoxide");
    assert!(matches!(molecule.smarts_tetrahedral_chirality(2), Some(Chirality::TH(1 | 2))));
}

#[test]
fn legacy_none_still_names_a_bracket_hydrogen_not_a_lone_pair() {
    let molecule: Smiles = "CC[S@](=O)CCC".parse().expect("valid sulfoxide");
    assert_eq!(
        molecule.tetrahedral_neighbors_anticlockwise(2, [Some(1), Some(3), Some(4), None]),
        None
    );
}

#[test]
fn sulfur_orientation_matches_native_neighbor_order_including_rings_and_roots() {
    #[derive(Deserialize)]
    struct Case {
        smiles: String,
        center: usize,
        neighbors: [usize; 3],
        anticlockwise: bool,
    }
    #[derive(Deserialize)]
    struct Fixture {
        cases: Vec<Case>,
    }
    let fixture: Fixture =
        serde_json::from_str(include_str!("fixtures/sulfur_tetrahedral_rdkit_2026_09_1.json"))
            .unwrap();
    assert_eq!(fixture.cases.len(), 24);
    for case in fixture.cases {
        let molecule: WildcardSmiles = case.smiles.parse().expect("declared SMILES");
        let [first, second, third] = case.neighbors;
        let slots = [Atom(first), Atom(second), Atom(third), LonePair];
        assert_eq!(
            molecule.tetrahedral_slots_anticlockwise(case.center, slots),
            Some(case.anticlockwise),
            "{}",
            case.smiles
        );
        let mut swapped = slots;
        swapped.swap(0, 3);
        assert_eq!(
            molecule.tetrahedral_slots_anticlockwise(case.center, swapped),
            Some(!case.anticlockwise)
        );
        if !case.smiles.contains('*') {
            let concrete: Smiles = case.smiles.parse().unwrap();
            assert_eq!(
                concrete.tetrahedral_slots_anticlockwise(case.center, slots),
                Some(case.anticlockwise)
            );
        }
    }
}

#[test]
fn lone_pairs_and_hydrogens_have_distinct_slot_identities() {
    let molecule: Smiles = "CC[S@](=O)CCC".parse().unwrap();
    let slots = [Atom(1), Atom(3), Atom(4), LonePair];
    assert_eq!(molecule.tetrahedral_stereo_neighbors(2), Some(slots));
    for wrong in [
        [Atom(1), Atom(3), Atom(4), BracketHydrogen],
        [Atom(1), Atom(3), LonePair, LonePair],
        [Atom(1), Atom(3), Atom(3), LonePair],
        [Atom(1), Atom(3), Atom(2), LonePair],
        [Atom(1), Atom(3), Atom(usize::MAX), LonePair],
    ] {
        assert_eq!(molecule.tetrahedral_slots_anticlockwise(2, wrong), None);
    }
    assert_eq!(molecule.tetrahedral_stereo_neighbors(usize::MAX), None);
    assert_eq!(molecule.tetrahedral_slots_anticlockwise(usize::MAX, slots), None);

    let hydrogen: Smiles = "C[S@H](=O)CC".parse().unwrap();
    let bracket_slots = hydrogen.tetrahedral_stereo_neighbors(1).unwrap();
    assert!(bracket_slots.contains(&BracketHydrogen));
    assert!(!bracket_slots.contains(&LonePair));
    let expanded = hydrogen.with_explicit_hydrogens();
    let graph_slots = expanded.tetrahedral_stereo_neighbors(1).unwrap();
    assert!(graph_slots.iter().all(|slot| matches!(slot, Atom(_))));
    assert_eq!(expanded.tetrahedral_slots_anticlockwise(1, graph_slots), Some(true));
}

#[test]
fn unsupported_lone_pair_families_and_coordination_classes_are_not_inferred() {
    for (input, center) in [
        ("CCS(=O)CCC", 2),
        ("C[S@](CC)CCC", 1),
        ("[S@+](C)(CC)CCC", 0),
        ("C[Se@](=O)CC", 1),
        ("C[P@](CC)CCC", 1),
        ("C[N@](CC)CCC", 1),
        ("CC[S@AL1](=O)CCC", 2),
        ("CC[S@SP1](=O)CCC", 2),
        ("F[Pt@SP1](Cl)(Br)I", 1),
        ("F[P@TB1](Cl)(Br)(I)N", 1),
        ("F[Co@OH1](Cl)(Br)(I)(N)O", 1),
    ] {
        let molecule: Smiles = input.parse().unwrap();
        assert_eq!(molecule.tetrahedral_stereo_neighbors(center), None, "{input}");
    }
    let sulfone: Smiles = "C[S@](=O)(=O)CC".parse().unwrap();
    assert!(!sulfone.tetrahedral_stereo_neighbors(1).unwrap().contains(&LonePair));
}

#[test]
fn parsed_slots_do_not_certify_stereogenicity() {
    for input in ["C[S@](=O)C", "O=[S@]1CCCC1"] {
        let molecule: Smiles = input.parse().unwrap();
        assert!(molecule.tetrahedral_stereo_neighbors(1).unwrap().contains(&LonePair));
        assert_eq!(molecule.smarts_tetrahedral_chirality(1), None);
    }
}
