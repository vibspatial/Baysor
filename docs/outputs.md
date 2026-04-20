# Outputs

The C++ branch currently supports two output styles:

- `legacy`
- `parquet`

Select them with:

```bash
./build/baysor run --output-style legacy ...
./build/baysor run --output-style parquet ...
```

`legacy` is the default.

For file-by-file definitions, including exact column order, optional fields,
and HDF5 / Parquet layout details, see [Output Files](output_files.md).

## Legacy Output

The legacy bundle mirrors the familiar Baysor output set:

- `segmentation.csv`
- `segmentation_cell_stats.csv`
- `segmentation_polygons_2d.json`
- `segmentation_polygons_3d.json`
- `segmentation_counts.loom` or `segmentation_counts.tsv`
- `segmentation_params.dump.toml`
- `segmentation_log.log`

Optional:

- `diagnostic_report.html`
- `segmentation_plot.html`

These files are defined in [Output Files](output_files.md#legacy-bundle).

### Xenium Compatibility In Legacy Output

For Xenium-origin inputs, `legacy` automatically adds the extra fields needed by
`xeniumranger import-segmentation`:

- `segmentation.csv`
  - includes `transcript_id`
  - writes `is_noise` as `true` / `false`
- `segmentation_polygons_2d.json`
  - uses GeoJSON `FeatureCollection`
  - includes `properties.cell`

This is the recommended output style when the result will be handed off to
Xenium Ranger / Xenium Explorer.

## Parquet Output

The `parquet` style is an interop-oriented bundle:

- `molecules.parquet`
- `cells.parquet`
- `cell_boundaries.parquet`
- `cell_boundaries_3d.parquet`
- `feature_matrix.h5`
- `run_params.toml`
- `run.log`

Optional:

- `diagnostic_report.html`
- `segmentation_plot.html`

These files are defined in [Output Files](output_files.md#parquet-bundle).

Boundary outputs use GeoParquet-compatible metadata.

The count matrix is written as a 10x-style HDF5 feature-barcode matrix.

## Legacy-Only Flags

These remain meaningful only for `legacy` output:

- `--polygon-format`
- `--count-matrix-format`

For `parquet`, they are ignored with a warning.

## Choosing Between Styles

Use `legacy` when:

- you want the familiar Baysor outputs
- you want Xenium Ranger compatibility
- you need GeoJSON polygons

Use `parquet` when:

- you want easier downstream analysis in Python / R / DuckDB
- you want GeoParquet boundaries
- you want a 10x-style HDF5 matrix bundle
