//! Statistics follow the supplied graph and its recorded hydrogen metadata.

use serde::Deserialize;
use smiles_rs::{AtomStatistics, Smiles};

#[derive(Deserialize)]
struct ExpectedCounts {
    graph_atoms: u64,
    bracket_hydrogens: u64,
    inferred_hydrogens: u64,
    atoms_including_hydrogens: u64,
    heavy_atoms: u64,
    carbon_atoms: u64,
    carbon_total_connectivity_four: u64,
    nitrogen_oxygen_atoms: u64,
    hydrogens_on_nitrogen_oxygen: u64,
    carbon_total_connectivity_four_fraction_bits: u64,
}

#[derive(Deserialize)]
struct Case {
    id: String,
    smiles: String,
    counts: ExpectedCounts,
    concrete_parser_supported: bool,
}

#[derive(Deserialize)]
struct Fixture {
    rows: Vec<Case>,
}

fn fixture() -> Fixture {
    serde_json::from_str(include_str!("fixtures/atom_statistics_raw_rdkit_2026_09_1.json")).unwrap()
}

fn counts_array(counts: AtomStatistics) -> [u64; 9] {
    [
        counts.graph_atoms(),
        counts.bracket_hydrogens(),
        counts.inferred_hydrogens(),
        counts.atoms_including_hydrogens(),
        counts.heavy_atoms(),
        counts.carbon_atoms(),
        counts.carbon_total_connectivity_four(),
        counts.nitrogen_oxygen_atoms(),
        counts.hydrogens_on_nitrogen_oxygen(),
    ]
}

fn assert_expected(actual: AtomStatistics, expected: &ExpectedCounts, id: &str) {
    assert_eq!(
        counts_array(actual),
        [
            expected.graph_atoms,
            expected.bracket_hydrogens,
            expected.inferred_hydrogens,
            expected.atoms_including_hydrogens,
            expected.heavy_atoms,
            expected.carbon_atoms,
            expected.carbon_total_connectivity_four,
            expected.nitrogen_oxygen_atoms,
            expected.hydrogens_on_nitrogen_oxygen,
        ],
        "{id}"
    );
    assert_eq!(
        actual.carbon_total_connectivity_four_fraction().to_bits(),
        expected.carbon_total_connectivity_four_fraction_bits,
        "{id}"
    );
}

#[test]
fn raw_counts_match_independently_recorded_metadata() {
    let fixture = fixture();
    assert_eq!(fixture.rows.len(), 41);
    let mut exclusions = Vec::new();
    for case in fixture.rows {
        let parsed = case.smiles.parse::<Smiles>();
        if !case.concrete_parser_supported {
            assert!(parsed.is_err(), "{}: exclusion changed", case.id);
            exclusions.push(case.smiles);
            continue;
        }
        let molecule = parsed.unwrap_or_else(|error| panic!("{}: {error}", case.id));
        let actual = molecule.atom_statistics().unwrap();
        assert_expected(actual, &case.counts, &case.id);
    }
    assert_eq!(exclusions, ["*", "[999C]"]);
}

#[test]
fn graph_hydrogen_incidence_includes_bridging_isotopic_and_charged_vertices() {
    let fixture: Fixture =
        serde_json::from_str(include_str!("fixtures/atom_statistics_graph_h_rdkit_2026_09_1.json"))
            .unwrap();
    assert_eq!(fixture.rows.len(), 4);
    for case in fixture.rows {
        let molecule: Smiles = case.smiles.parse().unwrap();
        let raw = molecule.atom_statistics().unwrap();
        assert_expected(raw, &case.counts, &case.id);
        assert_composition_unchanged(
            raw,
            molecule.with_explicit_hydrogens().atom_statistics().unwrap(),
        );
        assert_composition_unchanged(raw, molecule.canonicalize().atom_statistics().unwrap());
    }
}

fn assert_composition_unchanged(before: AtomStatistics, after: AtomStatistics) {
    assert_eq!(before.atoms_including_hydrogens(), after.atoms_including_hydrogens());
    assert_eq!(before.heavy_atoms(), after.heavy_atoms());
    assert_eq!(before.carbon_atoms(), after.carbon_atoms());
    assert_eq!(before.carbon_total_connectivity_four(), after.carbon_total_connectivity_four());
    assert_eq!(before.nitrogen_oxygen_atoms(), after.nitrogen_oxygen_atoms());
    assert_eq!(before.hydrogens_on_nitrogen_oxygen(), after.hydrogens_on_nitrogen_oxygen());
}

#[test]
fn materializing_hydrogens_moves_counts_into_vertices_without_double_counting() {
    for case in fixture().rows.into_iter().filter(|case| case.concrete_parser_supported) {
        let molecule: Smiles = case.smiles.parse().unwrap();
        let original = molecule.atom_statistics().unwrap();
        let expanded = molecule.with_explicit_hydrogens().atom_statistics().unwrap();
        assert_composition_unchanged(original, expanded);
        assert_eq!(expanded.graph_atoms(), original.atoms_including_hydrogens(), "{}", case.id);
        assert_eq!(expanded.bracket_hydrogens(), 0, "{}", case.id);
        assert_eq!(expanded.inferred_hydrogens(), 0, "{}", case.id);
    }
}

#[test]
fn ordinary_hydrogen_collapse_preserves_composition_and_keeps_isotope_charge_map() {
    for input in ["[H]N([H])[H]", "[2H]N([H])[H]", "[H:7]N", "[H-].[Na+]", "[NH4+]"] {
        let molecule: Smiles = input.parse().unwrap();
        let original = molecule.atom_statistics().unwrap();
        let canonical = molecule.canonicalize();
        let after = canonical.atom_statistics().unwrap();
        assert_composition_unchanged(original, after);
        if input == "[H]N([H])[H]" {
            assert_eq!(original.graph_atoms(), 4);
            assert_eq!(after.graph_atoms(), 1);
            assert_eq!(after.bracket_hydrogens() + after.inferred_hydrogens(), 3);
        }
    }
}

#[test]
fn raw_counts_do_not_validate_valence_or_cap_atom_totals() {
    let overvalent: Smiles = "C(C)(C)(C)(C)C".parse().unwrap();
    let counts = overvalent.atom_statistics().unwrap();
    assert_eq!(counts.carbon_atoms(), 6);
    assert_eq!(counts.carbon_total_connectivity_four(), 5);
    for size in [999, 1000, 1001] {
        let input = std::iter::repeat_n("O", size).collect::<Vec<_>>().join(".");
        let molecule: Smiles = input.parse().unwrap();
        let counts = molecule.atom_statistics().unwrap();
        let expected = u64::try_from(size).unwrap();
        assert_eq!(counts.nitrogen_oxygen_atoms(), expected);
        assert_eq!(counts.hydrogens_on_nitrogen_oxygen(), 2 * expected);
        assert_eq!(counts.atoms_including_hydrogens(), 3 * expected);
        assert_eq!(counts.carbon_total_connectivity_four_fraction().to_bits(), 0);
    }
}
