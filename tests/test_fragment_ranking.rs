//! Default numeric fragment-parent candidate contract.

use smiles_rs::Smiles;

fn candidates(input: &str) -> Vec<(usize, u32, f64)> {
    input
        .parse::<Smiles>()
        .unwrap()
        .fragment_parent_candidates()
        .unwrap()
        .iter()
        .map(|candidate| {
            (
                candidate.component_id(),
                candidate.atom_count_including_hydrogens(),
                candidate.pseudo_weight(),
            )
        })
        .collect()
}

#[test]
fn primary_count_overrides_pseudo_weight_and_clears_earlier_ties() {
    assert_eq!(candidates("O.N"), vec![(1, 4, 703.0)]);
    assert_eq!(candidates("O.O.N"), vec![(2, 4, 703.0)]);
    assert_eq!(candidates("N.O.O"), vec![(0, 4, 703.0)]);
    assert_eq!(candidates("C.[Xe]"), vec![(0, 5, 604.0)]);
}

#[test]
fn numeric_ties_remain_all_candidates_in_component_order() {
    assert_eq!(candidates("CO.NN"), vec![(0, 6, 1404.0), (1, 6, 1404.0)]);
    assert_eq!(candidates("CC.CC.CC"), vec![(0, 8, 1206.0), (1, 8, 1206.0), (2, 8, 1206.0)]);
    assert_eq!(candidates("[CH3:1][OH:2].[CH3:3][OH:4]"), vec![(0, 6, 1404.0), (1, 6, 1404.0)]);
}

#[test]
fn explicit_h_representation_is_not_silently_normalized() {
    assert_eq!(candidates("[H]C.C"), vec![(0, 5, 703.0)]);
    assert_eq!(candidates("[H]C([H])([H])[H].C"), vec![(0, 5, 1000.0)]);
}

#[test]
fn isotope_number_and_charge_contribute_without_a_mass_table() {
    assert_eq!(candidates("[13CH4].C"), vec![(0, 5, 617.0)]);
    assert_eq!(candidates("[CH3+].[CH3-]"), vec![(1, 4, 603.1)]);
}

#[test]
fn ranking_does_not_prefer_organic_or_mutate_the_parent() {
    let parent: Smiles = "[NH4+].C".parse().unwrap();
    let before = format!("{parent:?}");
    let best = parent.fragment_parent_candidates().unwrap();
    assert_eq!(best.len(), 1);
    assert_eq!(best[0].component_id(), 0);
    assert_eq!(best[0].atom_count_including_hydrogens(), 5);
    assert_eq!(best[0].pseudo_weight().to_bits(), 703.9_f64.to_bits());
    assert_eq!(format!("{parent:?}"), before);
    let fragment = parent.fragment_from_component(best[0].component_id()).unwrap();
    assert_eq!(fragment.parent_id(0), 0);
}

#[test]
fn empty_and_single_component_graphs_have_numeric_results() {
    let parent: Smiles = "CCO".parse().unwrap();
    let empty = parent.fragment_from_atoms([]).unwrap().into_smiles();
    assert!(empty.fragment_parent_candidates().unwrap().is_empty());
    assert_eq!(candidates("CCO"), vec![(0, 9, 2006.0)]);
}

#[test]
fn repeated_charge_retains_source_rounding_and_unresolved_ties() {
    // Independent replay of the pinned native source operations is one ULP
    // below regrouping the integer terms and dividing once by ten. These are
    // source-policy bits, not a measurement of the native internal accumulator.
    let input = concat!(
        "[NH3+]C[NH2+]C[NH2+]C[NH2+]C[NH2+]C[NH2+]C[NH2+]C[NH3+].",
        "C([NH2+]C[NH2+]C[NH2+]C[NH3+])[NH2+]C[NH2+]C[NH2+]C[NH3+]"
    );
    let actual = candidates(input)
        .into_iter()
        .map(|(id, count, weight)| (id, count, weight.to_bits()))
        .collect::<Vec<_>>();
    assert_eq!(actual, vec![(0, 47, 0x40c3_3399_9999_9999), (1, 47, 0x40c3_3399_9999_9999)]);
}
