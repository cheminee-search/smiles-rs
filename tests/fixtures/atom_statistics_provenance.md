# Raw atom-statistics fixture

`atom_statistics_raw_rdkit_2026_09_1.json` contains 41 declared development
inputs. The expected counts were generated with RDKit 2026.09.1, whose pinned
source is `fece8caa860bdf6c9c82bdb9253f22a0cf87150c`. The JSON records the
declaration and generator SHA256 hashes, version and exact reference phase.

The reference uses `SmilesParserParams` with `sanitize=False`, `removeHs=False`,
`parseName=False` and `allowCXSMILES=False`, followed by
`UpdatePropertyCache(strict=False)`. It sums the recorded atom metadata:

- graph atoms: the number of graph vertices;
- bracket and inferred H: `GetNumExplicitHs()` and `GetNumImplicitHs()`;
- all atoms: graph atoms plus bracket and inferred H;
- heavy atoms: atomic number greater than one;
- carbon-four atoms: carbon with `GetDegree() + bracket H + inferred H == 4`;
- N/O atoms: atomic number seven or eight;
- H on N/O: bracket and inferred H plus neighboring atomic-number-one vertices.

The carbon-four fraction is the ordinary floating-point quotient, or positive
zero without carbon; the fixture stores its binary64 bits. Graph H vertices
include ordinary, isotopic, charged and mapped hydrogen. Their attached-atom
incidences are separate from the all-atoms count.

All reference inputs parsed in that raw native mode. The concrete Rust parser
rejects `*` and `[999C]`; both rows and their reference counts remain in the
fixture with explicit exclusion flags. The test asserts the complete exclusion
list. The three 999/1000/1001-oxygen inputs are direct-counter boundaries, not
SMARTS-match counts. The overvalent carbon input deliberately demonstrates that
counting a supplied graph does not validate its valence.

This is bounded evidence about raw metadata arithmetic. It does not establish
sanitized descriptor parity or a fragment/tautomer/neutralization preparation
contract. Transform tests separately check the actual resulting graph state;
graph-vertex and non-vertex hydrogen subtotals can change when H is materialized
or collapsed. No native code or coefficient table is included in the calculator.

`atom_statistics_graph_h_rdkit_2026_09_1.json` adds four separately declared
development inputs with the same reference settings. `N[H]O` and `N[2H]O` have
one graph H vertex adjacent to both N and O: it contributes once to graph/all
atoms and twice to the N/O attachment count. `N[H-]` and `N[H+]` retain the
attached charged H vertex. These raw graph examples do not claim valid sanitized
valence. Their declaration and generator hashes are recorded separately.
