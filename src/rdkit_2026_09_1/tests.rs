use core::{error::Error, num::NonZeroU16};

use elements_rs::Element;
use molecular_formulas_native::rdkit_2026_09_1::{AtomMassInput, molecular_weights};

use super::{RecordedWeightError, project_records, recorded_weights, validate_shape};
use crate::{
    Smiles,
    atom::{Atom, atom_symbol::AtomSymbol},
};

fn carbon(hydrogens: u8) -> Atom {
    Atom::builder().with_symbol(AtomSymbol::Element(Element::C)).with_hydrogens(hydrogens).build()
}

fn same_bits(actual: super::MolecularWeights, expected: super::MolecularWeights) {
    assert_eq!(actual.average().to_bits(), expected.average().to_bits());
    assert_eq!(actual.exact().to_bits(), expected.exact().to_bits());
}

#[test]
fn projects_each_component_atom_in_recorded_order() {
    let graph = "CC(=O)[O-].[Na+]".parse::<Smiles>().unwrap();
    let records = project_records(graph.nodes(), graph.implicit_hydrogen_counts()).unwrap();
    let expected = [
        AtomMassInput {
            element: Element::C,
            isotope: None,
            formal_charge: 0,
            nonvertex_hydrogens: 3,
        },
        AtomMassInput {
            element: Element::C,
            isotope: None,
            formal_charge: 0,
            nonvertex_hydrogens: 0,
        },
        AtomMassInput {
            element: Element::O,
            isotope: None,
            formal_charge: 0,
            nonvertex_hydrogens: 0,
        },
        AtomMassInput {
            element: Element::O,
            isotope: None,
            formal_charge: -1,
            nonvertex_hydrogens: 0,
        },
        AtomMassInput {
            element: Element::Na,
            isotope: None,
            formal_charge: 1,
            nonvertex_hydrogens: 0,
        },
    ];
    assert_eq!(records, expected);
    same_bits(recorded_weights(&graph).unwrap(), molecular_weights(expected).unwrap());
}

#[test]
fn graph_hydrogens_keep_their_own_isotope_and_charge() {
    let graph = "[2H]N([H])[H]".parse::<Smiles>().unwrap();
    let records = project_records(graph.nodes(), graph.implicit_hydrogen_counts()).unwrap();
    assert_eq!(
        records.iter().map(|atom| atom.element).collect::<alloc::vec::Vec<_>>(),
        [Element::H, Element::N, Element::H, Element::H]
    );
    assert_eq!(records[0].isotope, NonZeroU16::new(2));
    assert!(records.iter().all(|atom| atom.nonvertex_hydrogens == 0));
    let charged = "[H-].[Na+]".parse::<Smiles>().unwrap();
    let records = project_records(charged.nodes(), charged.implicit_hydrogen_counts()).unwrap();
    assert_eq!(records[0].element, Element::H);
    assert_eq!(records[0].formal_charge, -1);
    assert_eq!(records[1].element, Element::Na);
    assert_eq!(records[1].formal_charge, 1);
    assert!(records.iter().all(|atom| atom.nonvertex_hydrogens == 0));
}

#[test]
fn isotope_and_bracket_hydrogens_remain_recorded_values() {
    let graph = "[13CH4]".parse::<Smiles>().unwrap();
    let expected = [AtomMassInput {
        element: Element::C,
        isotope: NonZeroU16::new(13),
        formal_charge: 0,
        nonvertex_hydrogens: 4,
    }];
    assert_eq!(project_records(graph.nodes(), graph.implicit_hydrogen_counts()).unwrap(), expected);
    same_bits(recorded_weights(&graph).unwrap(), molecular_weights(expected).unwrap());
}

#[test]
fn observes_each_transformed_graph_without_changing_it() {
    let graph = "C[C@H](O)C(=O)O".parse::<Smiles>().unwrap();
    for state in [graph.clone(), graph.with_explicit_hydrogens(), graph.canonicalize()] {
        let before = alloc::format!("{state:?}");
        let records = project_records(state.nodes(), state.implicit_hydrogen_counts()).unwrap();
        same_bits(recorded_weights(&state).unwrap(), molecular_weights(&records).unwrap());
        assert_eq!(before, alloc::format!("{state:?}"));
    }
}

#[test]
fn checks_cache_shape_before_projection() {
    let atoms = [carbon(0)];
    assert_eq!(
        project_records(&atoms, &[]),
        Err(RecordedWeightError::HydrogenCacheLength { atoms: 1, hydrogens: 0 })
    );
    assert_eq!(
        project_records(&atoms, &[0, 0]),
        Err(RecordedWeightError::HydrogenCacheLength { atoms: 1, hydrogens: 2 })
    );
    assert_eq!(validate_shape(0, 0), Ok(()));
    assert_eq!(validate_shape(u32::MAX as usize, u32::MAX as usize), Ok(()));
}

#[test]
#[cfg(target_pointer_width = "64")]
fn rejects_native_atom_count_without_large_allocation() {
    let count = u32::MAX as usize + 1;
    assert_eq!(
        validate_shape(count, count),
        Err(RecordedWeightError::Calculation(super::MolecularWeightError::AtomCountOverflow))
    );
    assert_eq!(
        validate_shape(count, 0),
        Err(RecordedWeightError::HydrogenCacheLength { atoms: count, hydrogens: 0 })
    );
}

#[test]
fn preserves_widened_hydrogen_sum_for_calculator_validation() {
    for (bracket, recorded, count) in
        [(128, 128, 256), (255, 127, 382), (255, 128, 383), (255, 255, 510)]
    {
        let records = project_records(&[carbon(bracket)], &[recorded]).unwrap();
        assert_eq!(records[0].nonvertex_hydrogens, count);
        let weights = molecular_weights(&records);
        if count > 382 {
            assert_eq!(
                weights,
                Err(super::MolecularWeightError::NonvertexHydrogensOutOfRange {
                    atom_index: 0,
                    count
                })
            );
        } else {
            assert!(weights.is_ok());
        }
    }
}

#[test]
fn missing_elements_and_zero_isotopes_are_explicit_errors() {
    let wildcard = Atom::builder().build();
    assert_eq!(
        project_records(&[wildcard], &[0]),
        Err(RecordedWeightError::MissingElement { atom_index: 0 })
    );
    let zero_isotope =
        Atom::builder().with_symbol(AtomSymbol::Element(Element::C)).with_isotope(0).build();
    assert_eq!(
        project_records(&[zero_isotope], &[0]),
        Err(RecordedWeightError::ZeroIsotope { atom_index: 0 })
    );
}

#[test]
fn completes_projection_before_checking_calculator_bounds() {
    let wildcard = Atom::builder().build();
    let zero_isotope =
        Atom::builder().with_symbol(AtomSymbol::Element(Element::C)).with_isotope(0).build();
    assert_eq!(
        project_records(&[carbon(255), zero_isotope], &[255, 0]),
        Err(RecordedWeightError::ZeroIsotope { atom_index: 1 })
    );
    assert_eq!(
        project_records(&[carbon(255), wildcard], &[255, 0]),
        Err(RecordedWeightError::MissingElement { atom_index: 1 })
    );
    assert_eq!(
        project_records(&[wildcard, zero_isotope], &[0, 0]),
        Err(RecordedWeightError::MissingElement { atom_index: 0 })
    );
}

#[test]
fn preserves_typed_calculation_source_and_empty_positive_zero() {
    let cause =
        super::MolecularWeightError::NonvertexHydrogensOutOfRange { atom_index: 0, count: 383 };
    let error = RecordedWeightError::Calculation(cause);
    assert_eq!(error.source().unwrap().downcast_ref::<super::MolecularWeightError>(), Some(&cause));
    assert!(RecordedWeightError::ZeroIsotope { atom_index: 0 }.source().is_none());
    let empty = project_records(&[], &[]).unwrap();
    let weights = molecular_weights(&empty).unwrap();
    assert_eq!(weights.average().to_bits(), 0);
    assert_eq!(weights.exact().to_bits(), 0);
}
