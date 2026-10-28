//! Preserve unspecified neighboring geometry when explicit H is retained.
use smiles_rs::Smiles;

#[test]
fn canonicalization_keeps_unspecified_middle_with_isotopic_h() {
    let source = r"F/C([2H])=C(/[2H])C([2H])=C([2H])C(/[2H])=C(\[2H])Cl";
    let count = |molecule: &Smiles| {
        (0..molecule.nodes().len()).flat_map(|a| {
            (a + 1..molecule.nodes().len()).map(move |b| (a, b))
        }).filter(|&(a, b)| molecule.double_bond_stereo_config(a, b).is_some()).count()
    };
    let parsed: Smiles = source.parse().unwrap();
    assert_eq!(count(&parsed), 2);
    let canonical = parsed.canonicalize();
    assert_eq!(count(&canonical), 2);
    let reparsed: Smiles = canonical.to_string().parse().unwrap();
    assert_eq!(count(&reparsed), 2);
    assert_eq!(canonical.to_string(), reparsed.canonicalize().to_string());
}
