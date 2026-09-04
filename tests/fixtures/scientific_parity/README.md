# Scientific parity fixtures

These fixtures protect the scientific orchestration extracted into
`baysor::run_segmentation(...)`. Their reference outputs were produced by the
pre-extraction Baysor CLI at commit
`f46a1e1dce1606d0ea644f4f8f1cf682597ba65c`.

The historical output is the independent oracle for two current entry points:

```text
                 recorded pre-extraction result
                      /                 \
                     v                   v
          current Baysor CLI     baysor::run_segmentation(...)
```

This distinction matters because both current entry points link the same
scientific implementation. Comparing them only with each other could miss a
regression introduced in shared code.

## Fixture matrix

| Fixture | Scientific behaviour exercised | Regression detected |
| --- | --- | --- |
| `duplicate_2d` | Two pairs of identical XY coordinates, including the random jitter used before graph construction | A change in legacy random-stream consumption or duplicate-coordinate graph handling |
| `clustering_2d` with `mrf.toml` | Seeded ICA initialization, MRF molecule clustering, molecule types, and their use by BMM | A change in MRF cluster assignment or its downstream segmentation |
| `clustering_2d` with `louvain.toml` | NCV basis construction, Louvain graph clustering, molecule types, and downstream BMM | A change in the Louvain/NCV scientific path |
| `clustering_2d` with `leiden.toml` | NCV basis construction, Leiden graph clustering, molecule types, and downstream BMM | A change in the Leiden/NCV scientific path |
| `three_d` | XYZ loading, 3D BMM dispatch, joined 2D boundaries, and boundaries for each Z layer | A change specific to 3D orchestration or boundary-stack materialization |

The fixtures are intentionally small and locked to seed `1` and one native
thread. They are regression fixtures, not representative biological datasets
or performance benchmarks.

## Tests represented by these fixtures

Every matrix row produces two parameterized GoogleTest cases:

- `CurrentCliMatchesPreExtractionReference` invokes the current CLI in a fresh
  process and checks it against the historical result. This detects changes in
  shared scientific code as well as CLI orchestration.
- `DirectOperationMatchesPreExtractionReference` invokes
  `baysor::run_segmentation(...)` in-process and checks it against the same
  historical result. This detects differences introduced by the extracted
  native operation.

Two focused tests complement that matrix:

- `DefaultAndExplicitSeedUseTheSameScientificStream` verifies both the legacy
  default random sequence and complete scientific results against an explicit
  seed `1`.
- `AreRepeatableAndScientificallyNonInterfering` verifies that requested NCV
  colours are well-formed and repeatable for the locked build and that
  producing them does not change scientific outputs.

## Semantic comparison contract

The tests compare scientific meaning rather than byte layout:

| Product | Comparison |
| --- | --- |
| Molecule table | Stable transcript identity, gene, XY or XYZ coordinates, confidence, assignment confidence, noise status, and cell partition |
| Molecule clustering | Cluster partition for MRF, Louvain, and Leiden cases |
| Cell statistics | Columns and numeric values after matching cells and clusters |
| Count matrix | Gene-by-cell counts after matching cell labels |
| Boundaries | Cell polygons independent of feature and vertex ordering, including named per-Z layers for 3D |
| Resolved settings | Scientific input, prior, segmentation, clustering, and neighborhood-composition settings |

Cell and molecule-cluster numbers are arbitrary identifiers. They are therefore
matched through a bijection: renumbering an unchanged partition is accepted,
whereas splitting or combining a partition fails. Numeric comparisons use
narrow, product-specific tolerances to allow serialization precision without
masking a changed assignment or geometry.

Exact NCV colour strings are deliberately absent from the historical oracle.
Their UMAP orientation and resulting `#RRGGBB` values have no stable biological
identity and can vary across numerical toolchains. The separate NCV test covers
their representation contract and verifies that this visualization path cannot
alter molecule clusters, assignments, confidence, counts, statistics, or
boundaries.

## Reference generation

Build the historical commit in a separate checkout and set `BASELINE_BAYSOR`
to that executable. From the current repository root, generate each case in a
fresh process with one OpenMP thread:

```bash
export OMP_NUM_THREADS=1
export OMP_DYNAMIC=FALSE

"${BASELINE_BAYSOR}" run \
  tests/fixtures/scientific_parity/duplicate_2d/molecules.csv \
  --config tests/fixtures/scientific_parity/duplicate_2d/config.toml \
  --output reference-output/duplicate_2d \
  --output-style legacy --polygon-format FeatureCollection \
  --count-matrix-format tsv --skip-ncv-color

for method in mrf louvain leiden; do
  "${BASELINE_BAYSOR}" run \
    tests/fixtures/scientific_parity/clustering_2d/molecules.csv \
    --config "tests/fixtures/scientific_parity/clustering_2d/${method}.toml" \
    --output "reference-output/${method}" \
    --output-style legacy --polygon-format FeatureCollection \
    --count-matrix-format tsv --skip-ncv-color
done

"${BASELINE_BAYSOR}" run \
  tests/fixtures/scientific_parity/three_d/molecules.csv \
  --config tests/fixtures/scientific_parity/three_d/config.toml \
  --output reference-output/three_d \
  --output-style legacy --polygon-format FeatureCollection \
  --count-matrix-format tsv --skip-ncv-color
```

Run this recipe twice into different empty output roots and verify that
`segmentation.csv`, cell statistics, counts, 2D polygons, and (for the 3D
case) the polygon stack are byte-identical. The first comment in
`segmentation_params.dump.toml` records invocation paths and is expected to
differ; remove only that comment when storing the file as
`resolved_params.toml`.

## Interpreting failures

- If both current entry points fail the same fixture, first investigate shared
  scientific code or random-state plumbing.
- If only the direct operation fails, investigate extracted orchestration,
  ownership, or result serialization.
- If only the CLI fails, investigate CLI option translation or CLI output
  materialization.
- A seed-test failure indicates that the documented legacy random-stream
  compatibility has changed.
- An NCV-colour test failure indicates either unstable colour generation or
  scientific coupling to an optional visualization product.

A platform-specific numerical difference must be characterized rather than
silenced with a broad tolerance. Do not regenerate a reference merely to make a
failing test pass. Reference replacement is appropriate only after the
scientific change has been explained, reviewed, and accepted; the baseline
revision and generation recipe must then be updated together.

Run the focused fixture and contract tests with:

```bash
ctest --preset tests --output-on-failure \
  -R '^(NativeBaseline|SegmentationOperation|LegacySeedCompatibility|NeighborhoodCompositionColors|N2b/)'
```
