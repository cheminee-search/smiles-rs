//! Opt-in nitrogen-oxide charge separation on immutable parsed graphs.

use elements_rs::Element;
use smiles_rs::{NitrogenOxideCleanupError, Smiles, bond::Bond};

#[test]
fn charge_separation_preserves_indexed_hydrogens_and_original_graph() {
    let molecule: Smiles = "O=N(C)(C)C".parse().unwrap();
    let original = format!("{molecule:?}");
    let cleaned: Smiles = molecule.try_with_nitrogen_oxide_cleanup().unwrap();

    assert_eq!(cleaned.nodes()[0].charge_value(), -1);
    assert_eq!(cleaned.nodes()[1].charge_value(), 1);
    assert_eq!(cleaned.edge_for_node_pair((0, 1)).unwrap().bond(), Bond::Single);
    assert_eq!(cleaned.implicit_hydrogen_counts(), &[0, 0, 3, 3, 3]);
    assert_eq!(format!("{molecule:?}"), original);
    assert_eq!(
        format!("{:?}", cleaned.try_with_nitrogen_oxide_cleanup().unwrap()),
        format!("{cleaned:?}")
    );

    let reparsed: Smiles = cleaned.to_string().parse().unwrap();
    for (element, charge) in [(Element::N, 1), (Element::O, -1)] {
        let atom = reparsed.nodes().iter().find(|atom| atom.element() == Some(element)).unwrap();
        assert_eq!(atom.charge_value(), charge);
        assert_eq!(atom.hydrogen_count(), 0);
    }
}

#[test]
fn true_nonmatches_remain_exact_owned_copies() {
    for input in [
        "C[N+](C)(C)[O-]",
        "[O-]=N(C)(C)C",
        "O=[N+](C)(C)C",
        "O=[N-](C)(C)C",
        "[O+]=N(C)(C)C",
        "CN=O",
        "O=N(C)C",
        "O=N(=O)(C)C",
        "CN=N#N",
        "CP(=O)(C)C",
        "O=Cl(=O)(=O)O",
        "[CH5]",
    ] {
        let molecule: Smiles = input.parse().unwrap();
        let cleaned = molecule.try_with_nitrogen_oxide_cleanup().unwrap();
        assert_eq!(format!("{cleaned:?}"), format!("{molecule:?}"), "{input}");
    }
}

#[test]
fn preserves_isotopes_maps_bracket_and_vertex_hydrogens_and_remote_stereo() {
    for input in [
        "[18O:8]=[N:7]([CH3:1])([CH3:2])[CH3:3]",
        "O=[15N](C)(C)C",
        "O=[NH](C)C",
        "[H]N(=O)(C)C",
        "[2H]N(=O)(C)C",
        "O=[NH2]C",
        "O=[NH3]",
        "N[C@@H](C)CN(=O)(C)C",
        "C/C=C/CN(=O)(C)C",
        "[H]/C(F)=C(/Cl)CN(=O)(C)C",
        "c1ccccc1N(=O)(C)C",
        "O=N(C)(C)C.O=N(CC)(C)C",
        "O=N(C)(C)CN(=O)(C)C",
    ] {
        let molecule: Smiles = input.parse().unwrap();
        let cleaned = molecule.try_with_nitrogen_oxide_cleanup().unwrap();
        assert_eq!(cleaned.nodes().len(), molecule.nodes().len());
        assert_eq!(cleaned.implicit_hydrogen_counts(), molecule.implicit_hydrogen_counts());
        let mut changed = 0;
        for (old, new) in molecule.nodes().iter().zip(cleaned.nodes()) {
            assert_eq!(old.symbol(), new.symbol());
            assert_eq!(old.isotope_mass_number(), new.isotope_mass_number());
            assert_eq!(old.class(), new.class());
            assert_eq!(old.hydrogen_count(), new.hydrogen_count());
            assert_eq!(old.aromatic(), new.aromatic());
            assert_eq!(old.chirality(), new.chirality());
            if old.charge_value() == new.charge_value() {
                assert_eq!(old, new);
            } else {
                changed += 1;
            }
        }
        assert!(changed >= 2, "{input}");
        assert_eq!(
            cleaned.nodes().iter().map(|a| i32::from(a.charge_value())).sum::<i32>(),
            molecule.nodes().iter().map(|a| i32::from(a.charge_value())).sum::<i32>()
        );
    }
}

#[test]
fn chooses_first_eligible_oxygen_without_skipping_unsupported_environments() {
    for input in [
        "[OH:8]=[N:7](C)=[O:9]",
        "[O:9]=[N:7](C)=[OH:8]",
        "[OH]=N(C)(C)C",
        "[H]O=N(C)(C)C",
        "[2H]O=N(C)(C)C",
        "CO=N(C)(C)C",
    ] {
        let molecule: Smiles = input.parse().unwrap();
        assert!(
            matches!(
                molecule.try_with_nitrogen_oxide_cleanup(),
                Err(NitrogenOxideCleanupError::UnsupportedOxygen { .. })
            ),
            "{input}"
        );
    }
    for input in ["[O:8]=[N:7](=[O:9])[CH3:1]", "[O:9]=[N:7](=[O:8])[CH3:1]"] {
        let molecule: Smiles = input.parse().unwrap();
        let cleaned = molecule.try_with_nitrogen_oxide_cleanup().unwrap();
        assert_eq!(cleaned.nodes()[0].charge_value(), -1);
        assert_eq!(cleaned.nodes()[2].charge_value(), 0);
    }
}

#[test]
fn ordinary_bond_precedes_closure_and_multiple_closure_candidates_are_explicit() {
    for input in ["O1.N(=O)(C)=1", "O1.N(C)(=O)=1", "[O:8]1.[N:7](=[O:9])(C)=1"] {
        let molecule: Smiles = input.parse().unwrap();
        let cleaned = molecule.try_with_nitrogen_oxide_cleanup().unwrap();
        assert_eq!(cleaned.nodes()[0].charge_value(), 0, "{input}");
        assert_eq!(cleaned.nodes().iter().filter(|a| a.charge_value() == -1).count(), 1);
    }
    for input in ["O2.O1.N(C)=2=1", "[O:8]2.[O:9]1.[N:7](C)=2=1", "[O:8]1.[O:9]2.[N:7](C)=2=1"] {
        let molecule: Smiles = input.parse().unwrap();
        assert!(
            matches!(
                molecule.try_with_nitrogen_oxide_cleanup(),
                Err(NitrogenOxideCleanupError::AmbiguousOxygenOrder { .. })
            ),
            "{input}"
        );
    }
}

#[test]
fn sole_closure_candidate_is_supported() {
    let molecule: Smiles = "O1.N(C)(C)(C)=1".parse().unwrap();
    assert_eq!(molecule.try_with_nitrogen_oxide_cleanup().unwrap().nodes()[0].charge_value(), -1);
}

#[test]
fn parser_exclusions_remain_explicit_and_are_not_cleanup_results() {
    for input in
        ["", "O2.O1.N(C)(=2)=1", "[O:8]2.[O:9]1.[N:7](C)(=2)=1", "[O:8]1.[O:9]2.[N:7](C)(=2)=1"]
    {
        assert!(input.parse::<Smiles>().is_err(), "{input}");
    }
}

#[test]
fn unsupported_stereo_aromatic_and_multiple_bonds_fail_atomically() {
    for input in [
        "O=[N@](C)(F)Cl",
        "O=[N@@](C)(F)Cl",
        "O=n1ccccc1",
        "C=N(=O)C",
        "O=N(C)(C)C.[OH]=N(C)(C)C",
        "O=N(C)(C)C.O=n1ccccc1",
    ] {
        let molecule: Smiles = input.parse().unwrap();
        let before = format!("{molecule:?}");
        assert!(molecule.try_with_nitrogen_oxide_cleanup().is_err(), "{input}");
        assert_eq!(format!("{molecule:?}"), before);
    }
}
