# Output Files

This page defines the files written by `baysor run`.

For bundle-level selection, see [Outputs](outputs.md).

## Conventions

- `N` = number of retained molecules after input filtering and cropping
- `C` = number of final cells in the segmentation
- `G` = number of genes after input filtering
- cell names are written as `cell_<1-based-id>`
- noise / unassigned molecules are written with cell name `0`

## Legacy Bundle

### `segmentation.csv`

Shape:
- `N` rows
- one row per retained molecule

Exact column order:
- `[transcript_id,] cell, gene, x, y [,z] [,confidence] [,cluster] [,ncv_color] [,assignment_confidence], is_noise`

Columns:
- `transcript_id`
  - present only for Xenium-origin inputs that preserve source transcript IDs
  - string
- `cell`
  - assigned cell name
  - string
  - `0` means noise / unassigned
- `gene`
  - gene name written from the current input panel
  - string
- `x`, `y`
  - molecule coordinates
  - float
- `z`
  - present for 3D inputs only
  - float
- `confidence`
  - molecule confidence estimated during noise modeling
  - float
- `cluster`
  - molecule-cluster label from the molecule clustering stage
  - integer
- `ncv_color`
  - per-molecule NCV color as hex RGB, for example `#4F9BD5`
  - string
- `assignment_confidence`
  - per-molecule posterior assignment confidence
  - float
- `is_noise`
  - last column, always present
  - for Xenium-origin inputs: literal `true` / `false`
  - otherwise: `1` / `0`

### `segmentation_cell_stats.csv`

Shape:
- `C` rows
- one row per final cell

Exact column order:
- `cell`
- `x`, `y`
- optional `z`
- optional `cluster`
- `n_transcripts`
- `density`
- `elongation`
- `area`
- `avg_confidence`
- optional `avg_assignment_confidence`
- optional `max_cluster_frac`
- optional `lifespan`

Columns:
- `cell`
  - cell name, for example `cell_17`
  - string
- `x`, `y`, optional `z`
  - centroid of molecules assigned to the cell
  - float
- `cluster`
  - dominant / assigned cell-level cluster label when available
  - float-valued column containing integer labels
- `n_transcripts`
  - number of molecules assigned to the cell
  - float-valued column containing counts
- `density`
  - `n_transcripts / area` when area is available
  - float
- `elongation`
  - ratio of principal covariance eigenvalues for the cell molecule cloud
  - float
- `area`
  - 2D convex-hull area used for cell summary statistics
  - float
- `avg_confidence`
  - mean molecule confidence within the cell
  - float
- `avg_assignment_confidence`
  - mean posterior assignment confidence within the cell
  - float
- `max_cluster_frac`
  - fraction of molecules in the most frequent molecule-cluster label within the cell
  - float
- `lifespan`
  - number of traced iterations for the component GUID when tracing is available
  - float-valued column containing integer values

### `segmentation_polygons_2d.json`

When `--polygon-format FeatureCollection`:
- root object type: `FeatureCollection`
- one feature per cell with a valid 2D polygon
- each feature has:
  - `id`: cell name
  - `geometry.type`: `Polygon`
  - `geometry.coordinates`: one closed outer ring in data coordinates
  - `properties.cell`: cell name

When `--polygon-format GeometryCollection`:
- root object type: `GeometryCollection`
- one geometry object per cell
- each geometry stores:
  - `type: Polygon`
  - `coordinates`
  - `cell`

When `--polygon-format none`:
- file is omitted

### `segmentation_polygons_3d.json`

Shape:
- JSON object keyed by layer name

Structure:
- each key is a layer name such as `"z_003"` or another layer label emitted by polygon estimation
- each value is a 2D polygon collection in the same schema as `segmentation_polygons_2d.json`
- omitted when no per-layer polygon stack is written

### `segmentation_counts.loom`

HDF5 layout written by the current C++ implementation:
- `/matrix`
  - dataset type: `float32`
  - shape: `(G, C)`
- `/attrs/LOOM_SPEC_VERSION`
  - variable-length UTF-8 string
  - value: `"3.0.0"`
- `/row_attrs/Name`
  - UTF-8 string array
  - length: `G`
  - gene names
- `/col_attrs/Name`
  - UTF-8 string array
  - length: `C`
  - cell names
- `/col_attrs/CellID`
  - `float64`
  - length: `C`
  - values `1, 2, ..., C`
- optional extra `/col_attrs/<key>`
  - written when additional column attributes are supplied
  - string arrays or float64 arrays

Notes:
- matrix values are counts
- rows are genes and columns are cells, matching Loom row/column attributes

### `segmentation_counts.tsv`

Shape:
- dense tab-separated matrix

Layout:
- first column header: `gene`
- remaining column headers: one column per cell name
- one data row per gene
- matrix shape on disk: `(G rows) x (1 + C columns)`

This is the text alternative to the Loom count matrix.

### `segmentation_params.dump.toml`

Shape:
- TOML document

Contents:
- resolved run parameters after config-file loading and CLI overrides

### `segmentation_log.log`

Shape:
- line-oriented text log

Contents:
- stage progress
- iteration summaries
- convergence messages
- save-step messages

### `diagnostic_report.html`

Shape:
- self-contained HTML document

Contents:
- diagnostic plots and summary panels for `run --plot`

### `segmentation_plot.html`

Shape:
- self-contained HTML document

Contents:
- interactive molecule / cell visualization for `run --plot`

## Parquet Bundle

### `molecules.parquet`

Shape:
- `N` rows
- one row per retained molecule

Exact column order:
- `cell, gene, x, y [,z] [,confidence] [,cluster] [,ncv_color] [,assignment_confidence], is_noise`

Column types:
- `cell`: UTF-8 string
- `gene`: UTF-8 string
- `x`, `y`, optional `z`: float64
- `confidence`: float64
- `cluster`: int32
- `ncv_color`: UTF-8 string
- `assignment_confidence`: float64
- `is_noise`: boolean

Notes:
- unlike Xenium-friendly legacy CSV output, `transcript_id` is not currently written here

### `cells.parquet`

Shape:
- `C` rows
- one row per final cell

Columns:
- `cell`
- then the same cell-stat columns and order as `segmentation_cell_stats.csv`

Column types:
- `cell`: UTF-8 string
- all statistic columns: float64

### `cell_boundaries.parquet`

Shape:
- one row per cell with a valid 2D polygon

Columns:
- `cell`
  - UTF-8 string
- `n_vertices`
  - int32
  - number of polygon vertices excluding the duplicated closing vertex
- `geometry`
  - binary WKB polygon

GeoParquet metadata:
- primary geometry column: `geometry`
- encoding: `WKB`
- geometry type: `Polygon`
- CRS: unset / `null`

### `cell_boundaries_3d.parquet`

Shape:
- one row per `(cell, layer)` polygon in the polygon stack

Columns:
- `cell`
  - UTF-8 string
- `layer`
  - UTF-8 string
- `n_vertices`
  - int32
- `geometry`
  - binary WKB polygon

GeoParquet metadata:
- same `geometry` metadata as `cell_boundaries.parquet`

### `feature_matrix.h5`

10x-style HDF5 layout:
- `/matrix/barcodes`
  - UTF-8 string array
  - length: `C`
  - cell names
- `/matrix/data`
  - int32
  - length: `nnz`
  - nonzero values of the sparse matrix
- `/matrix/indices`
  - int32
  - length: `nnz`
  - row indices into the feature axis
- `/matrix/indptr`
  - int64
  - length: `C + 1`
  - column pointer array
- `/matrix/shape`
  - int64
  - length: `2`
  - value: `[G, C]`
- `/matrix/features/id`
  - UTF-8 string array
  - length: `G`
- `/matrix/features/name`
  - UTF-8 string array
  - length: `G`
- `/matrix/features/feature_type`
  - UTF-8 string array
  - length: `G`
  - current value for every gene: `Gene Expression`
- `/matrix/features/genome`
  - UTF-8 string array
  - length: `G`
  - current value for every gene: empty string

Notes:
- the matrix is stored as CSC after transposing the internal `C x G` count matrix to `G x C`
- current implementation writes `id = name = gene_names`

### `run_params.toml`

Shape:
- TOML document

Contents:
- resolved run parameters for the Parquet bundle

### `run.log`

Shape:
- line-oriented text log

Contents:
- run progress and stage messages for the Parquet bundle

### `diagnostic_report.html`

Shape:
- self-contained HTML document

Contents:
- diagnostic plots and summary panels for `run --plot`

### `segmentation_plot.html`

Shape:
- self-contained HTML document

Contents:
- interactive molecule / cell visualization for `run --plot`
