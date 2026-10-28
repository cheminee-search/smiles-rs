//! Canonicalization must preserve alkene geometry while collapsing ordinary H.

use smiles_rs::{DoubleBondStereoConfig, Smiles, WildcardSmiles};

fn configurations(molecule: &Smiles) -> Vec<DoubleBondStereoConfig> {
    let mut result = Vec::new();
    for left in 0..molecule.nodes().len() {
        for right in left + 1..molecule.nodes().len() {
            if let Some(configuration) = molecule.double_bond_stereo_config(left, right) {
                result.push(configuration);
            }
        }
    }
    result
}

fn assert_preserved(source: &str, implicit: &str, expected: DoubleBondStereoConfig) {
    let parsed: Smiles = source.parse().unwrap();
    assert_eq!(configurations(&parsed), vec![expected]);
    let canonical = parsed.canonicalize();
    assert_eq!(canonical.nodes().len(), 4, "ordinary H vertices should collapse");
    assert_eq!(configurations(&canonical), vec![expected], "{source}");
    assert_eq!(
        canonical.to_string(),
        implicit.parse::<Smiles>().unwrap().canonicalize().to_string()
    );
    let reparsed: Smiles = canonical.to_string().parse().unwrap();
    assert_eq!(configurations(&reparsed), vec![expected]);
    assert_eq!(canonical.to_string(), reparsed.canonicalize().to_string());
}

#[test]
fn canonicalization_preserves_e_when_directional_hydrogen_is_removed() {
    assert_preserved(r"Br/C([H])=C(\[H])F", "F/C=C/Br", DoubleBondStereoConfig::E);
}

#[test]
fn canonicalization_preserves_z_when_directional_hydrogen_is_removed() {
    assert_preserved(r"C(=C(\[H])F)(/Br)[H]", r"F/C=C\Br", DoubleBondStereoConfig::Z);
}

#[test]
fn coupled_alkenes_keep_all_four_configurations() {
    let cases = [
        (r"[H]/C(F)=C(/C([H])=C(\[H])Cl)[H]", r"F/C=C\C=C\Cl"),
        (r"[H]/C(F)=C(/C([H])=C(/[H])Cl)[H]", r"F/C=C\C=C/Cl"),
        (r"[H]/C(F)=C(\C([H])=C(\[H])Cl)[H]", r"F/C=C/C=C\Cl"),
        (r"[H]/C(F)=C(\C([H])=C(/[H])Cl)[H]", "F/C=C/C=C/Cl"),
    ];
    let mut canonical_strings = std::collections::BTreeSet::new();
    for (source, implicit) in cases {
        let parsed: Smiles = source.parse().unwrap();
        assert_eq!(configurations(&parsed).len(), 2);
        let canonical = parsed.canonicalize();
        assert_eq!(configurations(&canonical).len(), 2, "{source}");
        assert_eq!(
            canonical.to_string(),
            implicit.parse::<Smiles>().unwrap().canonicalize().to_string(),
            "{source}"
        );
        canonical_strings.insert(canonical.to_string());
    }
    assert_eq!(canonical_strings.len(), 4);
}

#[test]
fn replacing_a_removed_reference_reverses_its_direction() {
    for (source, implicit) in
        [("[H]/C(*)=C(/F)Cl", r"*/C=C(\F)Cl"), (r"[H]/C(*)=C(\F)Cl", r"*/C=C(/F)Cl")]
    {
        let parsed: WildcardSmiles = source.parse().unwrap();
        // H outranks the wildcard in the source graph, so this exercises the
        // reference replacement, not just the loss of a directional H bond.
        assert_eq!(parsed.double_bond_stereo_atoms(1, 3), Some((0, 5)));
        let canonical = parsed.canonicalize();
        assert_eq!(canonical.nodes().len(), parsed.nodes().len() - 1);
        assert_eq!(
            canonical.to_string(),
            implicit.parse::<WildcardSmiles>().unwrap().canonicalize().to_string()
        );
        assert_eq!(
            canonical.to_string(),
            canonical.to_string().parse::<WildcardSmiles>().unwrap().canonicalize().to_string()
        );
    }
}

#[test]
fn special_hydrogens_keep_their_metadata() {
    for source in [
        r"[2H]/C(Br)=C(\[H])F",
        r"[2H]/C(Br)=C(\[3H])F",
        r"[H:7]/C(Br)=C(\[H])F",
        r"[H:7]/C(Br)=C(\[H:8])F",
        r"[H+]/C(Br)=C(\[H])F",
        r"[H-]/C(Br)=C(\[H])F",
    ] {
        let parsed: Smiles = source.parse().unwrap();
        let metadata = |molecule: &Smiles| {
            let mut result = molecule
                .nodes()
                .iter()
                .filter(|atom| {
                    atom.element().is_some_and(|element| u8::from(element) == 1)
                        && (atom.isotope_mass_number().is_some()
                            || atom.charge_value() != 0
                            || atom.class() != 0)
                })
                .map(|atom| (atom.isotope_mass_number(), atom.charge_value(), atom.class()))
                .collect::<Vec<_>>();
            result.sort_unstable();
            result
        };
        let canonical = parsed.canonicalize();
        assert_eq!(metadata(&parsed), metadata(&canonical), "{source}");
        assert_eq!(configurations(&canonical).len(), 1, "{source}");
    }
}

#[test]
fn collapse_does_not_invent_unspecified_or_small_ring_stereo() {
    for source in [
        "[H]C(Br)=C([H])F",
        "[H]/C(Br)=C(F)[H]",
        r"[H]/C([H])=C(\[H])F",
        r"[H]/C1=C(\F)CCCC1",
        "[H]/CBr",
    ] {
        let parsed: Smiles = source.parse().unwrap();
        assert!(configurations(&parsed.canonicalize()).is_empty(), "{source}");
    }
}

#[test]
fn partially_specified_chains_preserve_unknown_geometry() {
    for (source, expected) in [
        (r"F/C([H])=C(/[H])C([H])=C(\[2H])Cl", 1),
        (r"F/C([H])=C(/[H])C([H])=C([H])C(/[H])=C(\[H])Cl", 2),
        (r"F/C([2H])=C(/[2H])C([2H])=C([2H])C(/[2H])=C(\[2H])Cl", 2),
    ] {
        let parsed: Smiles = source.parse().unwrap();
        assert_eq!(configurations(&parsed).len(), expected);
        let canonical = parsed.canonicalize();
        assert_eq!(configurations(&canonical).len(), expected, "{source}");
        let reparsed: Smiles = canonical.to_string().parse().unwrap();
        assert_eq!(configurations(&reparsed).len(), expected, "{source}");
        assert_eq!(canonical.to_string(), reparsed.canonicalize().to_string());
    }
}

#[test]
fn partial_geometry_retains_h_only_in_the_affected_component() {
    let source = r"F/C([H])=C(/[H])C([H])=C([H])C(/[H])=C(\[H])Cl.[H]C";
    let parsed: Smiles = source.parse().unwrap();
    let canonical = parsed.canonicalize();
    let ordinary_h_vertices = |molecule: &Smiles| {
        molecule
            .nodes()
            .iter()
            .filter(|atom| {
                atom.element().is_some_and(|element| u8::from(element) == 1)
                    && atom.isotope_mass_number().is_none()
            })
            .count()
    };
    assert_eq!(ordinary_h_vertices(&parsed), 7);
    assert_eq!(ordinary_h_vertices(&canonical), 6);
    assert_eq!(configurations(&canonical).len(), 2);
    assert_eq!(canonical.to_string(), canonical.canonicalize().to_string());
}

#[test]
fn partial_geometry_is_invariant_across_native_rooted_spellings() {
    // Generated independently with RDKit 2026.09.1, H retained. The fixture
    // records every atom root, exact-string deduplication and generator hash.
    let fixture: serde_json::Value =
        serde_json::from_str(include_str!("fixtures/partial_ez_roots_rdkit_2026_09_1.json"))
            .unwrap();
    let seeds = fixture["seeds"].as_array().unwrap();
    let variants = fixture["variants"].as_array().unwrap();
    assert_eq!(variants.len(), 38);
    for variant in variants {
        let input = variant["input"].as_str().unwrap();
        let canonical = input.parse::<Smiles>().unwrap().canonicalize();
        for origin in variant["origins"].as_array().unwrap() {
            let seed = seeds.iter().find(|seed| seed["id"] == origin["seed_id"]).unwrap();
            let source: Smiles = seed["input"].as_str().unwrap().parse().unwrap();
            assert_eq!(canonical.to_string(), source.canonicalize().to_string(), "{input}");
            assert_eq!(configurations(&canonical).len(), configurations(&source).len());
        }
        assert_eq!(
            canonical.to_string(),
            canonical.to_string().parse::<Smiles>().unwrap().canonicalize().to_string(),
            "{input}"
        );
    }
}

#[test]
fn hydrogen_collapse_preserves_original_inferred_and_bracket_h() {
    let hydrogen_count = |molecule: &Smiles| {
        molecule
            .nodes()
            .iter()
            .enumerate()
            .map(|(id, atom)| {
                usize::from(atom.hydrogen_count())
                    + usize::from(molecule.implicit_hydrogen_count(id))
                    + usize::from(atom.element().is_some_and(|element| u8::from(element) == 1))
            })
            .sum::<usize>()
    };
    for source in [
        "[H]C",
        "[H]CC",
        "[H]/CBr",
        "[H][C]",
        "[H][CH]",
        "[H][C+]",
        "[H][N+](C)(C)C",
        "[H][O+]([H])C",
        "[H][O-]",
        "[H:7]C",
        "[2H]C",
        "[H]C([H])([H])[H]",
    ] {
        let parsed: Smiles = source.parse().unwrap();
        let canonical = parsed.canonicalize();
        assert_eq!(hydrogen_count(&parsed), hydrogen_count(&canonical), "{source}");
        let reparsed: Smiles = canonical.to_string().parse().unwrap();
        assert_eq!(hydrogen_count(&parsed), hydrogen_count(&reparsed), "{source}");
        assert_eq!(canonical.to_string(), reparsed.canonicalize().to_string());
    }
}
