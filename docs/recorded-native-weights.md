# Optional recorded native-vintage weights

The `rdkit-2026-09-1-mass` feature adds
`smiles_rs::rdkit_2026_09_1::recorded_weights(&Smiles)`. It calculates the
weights of the current graph's recorded atoms in node order. Graph H vertices
remain separate records; bracket H plus the existing implicit-H cache supply
nonvertex H. No preparation, sanitization, H inference, fragment selection or
formula aggregation is performed.

This local branch pairs with the existing molecular-formulas fork at commit
`81c1b24af252ddce2a9551a952e33d692fbccbc1`, checked out as the sibling directory
`molecular-formulas-native-mass`. That dependency provides the versioned mass
tables and checked arithmetic under its MIT AND BSD-3-Clause license and bundled
notices. Its package version is 0.3.0; this branch does not claim that the registry
release contains the local addition. The relative path is explicit in
Cargo.toml, and a portable checkout must preserve this pairing. Publication or a
registry dependency switch requires a separately qualified source reference.

The existing registry molecular-formulas 0.1.10 dependency, scientific methods and
SMILES-to-ChemicalFormula conversions are unchanged. Enabling this feature can
compile both formula versions; actual package identity and feature resolution
must be recorded for a consumer. The optional bridge reexports only its
versioned result/error types, not a second unqualified ChemicalFormula.

Projection uses a temporary Vec and preserves typed element/isotope/count
errors. Cache shape and native atom-count bounds are checked before allocation;
all records are projected before the calculator checks H limits. A later
projection error consequently precedes an earlier record's calculator error.
Isotope zero is rejected rather than interpreted as natural abundance. Native
H/count limits belong to the calculator and do not establish sanitizability.

Raw, H-expanded and canonicalized graphs are different recorded states. Evaluate
each separately; floating-point rounding need not be invariant when order or H
representation changes. Success is not Cheminee prepared amw/exactmw parity, and
there is no zero-allocation claim. The new helper does not change any current
preparation or default behavior.

This addition is being qualified independently. Source tests and metadata alone
do not establish the final compiled closure, native parity or production
integration.
