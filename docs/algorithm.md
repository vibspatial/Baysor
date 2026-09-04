# How Baysor Segmentation Works

Baysor assigns individual spatial transcripts to cells or to noise. It combines
transcript coordinates, gene identities, local molecule neighborhoods, and an
optional prior segmentation in a spatial mixture model. A cell is represented
by both a spatial distribution and a gene-composition distribution.

This page is an implementation-oriented overview of the native Baysor
segmentation pipeline. It explains how the main stages fit together and points
to their implementations. For the scientific method and its evaluation, see the
[Baysor paper](citation.md). The source code remains authoritative for exact
numerical behavior.

## Pipeline at a Glance

```text
transcript table + optional prior segmentation
                    |
                    v
          load, validate, and filter
                    |
                    v
       resolve scale and initial cell count
                    |
                    v
       estimate molecule signal confidence
                    |
                    v
        build the molecule-neighborhood graph
                    |
                    v
        optional molecule-type clustering
                    |
                    v
       initialize spatial cell components
                    |
                    v
     iteratively update assignments and cells
                    |
                    v
       refine the assignment through history
                    |
                    v
   assignments, statistics, counts, and boundaries
```

The native library entry point is `baysor::run_segmentation(...)`, implemented
in [`src/segmentation/segmentation.cpp`](../src/segmentation/segmentation.cpp).
It returns an owned scientific result. Serialization and presentation happen
after this operation and are not part of the segmentation model.

## 1. Load and Prepare the Molecule Cloud

The input is a CSV or Parquet transcript table, or a Xenium experiment resolved
to its transcript table. Each retained molecule has at least:

- an `x` and `y` coordinate;
- an optional `z` coordinate;
- a gene identity; and
- optional source identity, quality, cluster, or prior-segmentation fields.

Loading applies configured column mappings, spatial crops, gene exclusions,
quality filters, and minimum-frequency filters. These operations define the
molecule cloud seen by every later stage. Stable source transcript identifiers
are retained where the input provides them.

The loaders live under [`src/data_loading`](../src/data_loading). For accepted
formats and columns, see [Input Data](inputs.md).

## 2. Incorporate an Optional Prior Segmentation

A prior can come from a molecule-table column, a labelled image, or cell
boundaries. Baysor converts it into a per-molecule prior label. The prior guides
the model; it is not simply copied into the output.

It can influence several parts of a run:

- estimate cell scale when no explicit scale is supplied;
- adjust the automatically inferred initial cell count;
- place a lower bound on signal confidence for molecules covered by the prior;
  and
- favor assignments that remain coherent with prior segments while still
  permitting corrections.

The strength of the last two effects is controlled by
`prior_segmentation_confidence`. A larger value gives the prior more influence,
but Baysor continues to combine it with spatial and gene-composition evidence.

See [Prior Segmentation Inputs](priors.md) and
[`src/data_loading/prior_segmentation.cpp`](../src/data_loading/prior_segmentation.cpp).

## 3. Resolve Data-Dependent Parameters

Some settings cannot be resolved until the retained molecule cloud and prior
have been loaded.

### Cell scale

`scale` describes the approximate cell radius and regularizes the spatial shape
of inferred cells. It can be supplied explicitly or estimated from a prior
segmentation. `scale_std` controls the permitted variation around it.

### Initial number of cells

When `n_cells_init` is not positive, Baysor first computes:

```text
default_initial_cells =
    2 * floor(number_of_molecules / min_molecules_per_cell)
```

If a prior is present, it also computes a prior-aware estimate from the number
of active prior segments and unassigned molecules. The selected count is capped
by that estimate, but never below the number of active prior segments. This is
an initialization heuristic: components can subsequently disappear, and
disconnected assignment fragments can be returned to noise during fitting.

Defaults and validation are implemented in
[`src/utils/options.cpp`](../src/utils/options.cpp). The prior-aware cell-count
calculation and the complete run-level resolution are in
[`src/segmentation/segmentation.cpp`](../src/segmentation/segmentation.cpp).

## 4. Estimate Molecule Confidence and Noise

Baysor estimates whether each molecule looks like local signal or spatial
background before fitting cells.

For each molecule, it measures a local nearest-neighbor distance. It then fits a
two-component model to those distances: a denser signal population and a more
dispersed noise population. A molecule's signal posterior becomes its
confidence. Neighboring molecules influence this fit through an
MRF-regularized expectation step, and a prior segmentation can provide a
minimum confidence for covered molecules.

Confidence is subsequently used to balance the probability of assigning a
molecule to a cell against assigning it to noise. It is not the same as the
final assignment confidence, which is derived later from assignment history.

This stage is implemented in
[`src/processing/data_processing/noise_estimation.cpp`](../src/processing/data_processing/noise_estimation.cpp).

## 5. Build the Molecule-Neighborhood Graph

Baysor builds a weighted spatial adjacency graph between nearby molecules. Edge
weights decrease with spatial distance.

The graph serves two related purposes:

- it supplies the local Markov-random-field evidence used during assignment;
  and
- it limits the candidate cells for a molecule to cells already represented in
  its local neighborhood, plus the noise state where applicable.

This keeps the fitting process spatially coherent and avoids comparing every
molecule with every cell component. Graph construction is implemented in
[`src/processing/data_processing/initialization.cpp`](../src/processing/data_processing/initialization.cpp),
with graph data stored by
[`include/baysor/processing/models/adj_list.h`](../include/baysor/processing/models/adj_list.h).

## 6. Optionally Estimate Coarse Molecule Types

Molecule clustering supplies a coarse gene-composition compatibility signal to
the cell model. It does not segment cells by itself.

Supported methods are:

- `mrf`, which initializes gene profiles with ICA and fits a categorical
  mixture regularized by the molecule graph;
- `louvain` and `leiden`, which build neighborhood-composition vectors, cluster
  an anchor graph, and transfer the resulting labels to all molecules; and
- `none`, which disables this prior.

During cell fitting, a molecule is penalized when its coarse type conflicts with
the dominant type of a candidate cell. The clustering seed is independent of
the core segmentation stream so that its randomness is explicit and
reproducible under the documented execution contract.

The dispatcher and graph-clustering implementations are in
[`src/processing/bmm_algorithm/molecule_clustering_louvain.cpp`](../src/processing/bmm_algorithm/molecule_clustering_louvain.cpp).
The ICA/MRF implementation is in
[`src/processing/bmm_algorithm/molecule_clustering.cpp`](../src/processing/bmm_algorithm/molecule_clustering.cpp).

## 7. Initialize Cell Components

Each candidate cell is represented by a component containing:

- a multivariate normal spatial distribution in two or three dimensions;
- a smoothed categorical distribution over genes;
- its current molecule count and mixture weight;
- optional prior-segment bookkeeping; and
- an optional dominant coarse molecule type.

Initial centers are distributed across the molecule cloud. Initial spatial
covariances are regularized toward the configured cell scale, while gene
profiles start from an uninformative distribution and are learned from assigned
molecules.

Initialization is implemented by `initialize_bmm_data<N>(...)` in
[`src/processing/data_processing/initialization.cpp`](../src/processing/data_processing/initialization.cpp).
The cell-component contract is in
[`include/baysor/processing/models/component.h`](../include/baysor/processing/models/component.h).

## 8. Fit Assignments and Cell Components

The central algorithm alternates between assignment and component updates.

### Assignment step

For every molecule, Baysor evaluates the candidate cells represented in its
graph neighborhood. A candidate's weight combines:

- the cell's current mixture weight;
- spatial likelihood under the cell's multivariate normal distribution;
- gene likelihood under the cell's smoothed composition distribution;
- agreement with neighboring molecule assignments;
- optional prior-segmentation evidence;
- optional coarse molecule-type compatibility; and
- the molecule's estimated signal confidence.

When the molecule's confidence is below one, noise is also a candidate. The
next assignment is sampled from the resulting weights rather than always taking
the largest one. This stochastic update is one reason segmentation has an
explicit random seed.

### Component-update step

After applying the new assignments, Baysor re-estimates each cell's:

- center and spatial covariance;
- smoothed gene-composition distribution;
- molecule count and mixture weight; and
- dominant coarse molecule type, when clustering is enabled.

The covariance estimate is regularized toward the configured cell-scale prior,
which stabilizes small or irregular components.

### Component maintenance and convergence

During fitting, Baysor periodically removes disconnected fragments from a cell
assignment and returns those molecules to noise. Components containing fewer
than two molecules are dropped in the compatibility implementation. Assignment
history and component-count traces are recorded for refinement and diagnostics.

When a positive tolerance is configured, convergence is reached when the worst
fraction of changed assignments across the latest 20 iterations is below that
tolerance. A final refinement selects assignments using their recent history,
produces assignment-confidence values, and updates the surviving components.

The loop is implemented in
[`src/processing/bmm_algorithm/bmm_algorithm.cpp`](../src/processing/bmm_algorithm/bmm_algorithm.cpp).
Spatial and gene likelihoods are implemented under
[`src/processing/distributions`](../src/processing/distributions).

## 9. Materialize Scientific Results

The fitted state is converted into an owning `SegmentationResult`. Cell label
`0` represents noise; positive labels identify inferred cells.

The result can contain:

- retained molecules and their stable source identities;
- cell/noise assignments and molecule signal confidence;
- assignment confidence derived from recent assignment history;
- optional coarse molecule-cluster labels;
- per-cell statistics;
- a sparse cell-by-gene count matrix;
- two-dimensional or z-layer boundary polygons;
- neighborhood-composition colors; and
- report-neutral diagnostics and run provenance.

The assignments are the central segmentation result. Statistics, count matrices,
boundaries, colors, and diagnostics are derived from the fitted assignments or
the retained molecule cloud.

In particular, boundary estimation happens after model fitting. It reconstructs
polygons from the complete set of molecule coordinates and final assignments;
those polygons do not feed back into the BMM. The implementation is in
[`src/processing/data_processing/boundary_estimation.cpp`](../src/processing/data_processing/boundary_estimation.cpp).

Neighborhood-composition colors are also post-segmentation outputs. Their UMAP
embedding is stochastic, but the dedicated color stream prevents requesting
colors from consuming the core assignment stream. The color implementation is
in [`src/reporting/color_utils.cpp`](../src/reporting/color_utils.cpp).

For the exact native result and ownership contract, see the
[Native Segmentation API Contract](../contract.md). For file formats, see
[Outputs](outputs.md) and [Output Files](output_files.md).

## Two-Dimensional and Three-Dimensional Runs

The same mixture-model loop is instantiated for two or three spatial
dimensions. A 3D run includes `z` in component centers, covariance estimates,
and spatial likelihoods.

Boundary products remain reporting geometries. Baysor produces a joined
two-dimensional boundary collection and, for 3D inputs, can additionally
estimate boundary collections for individual or binned z layers. The boundary
representation should not be interpreted as the full volumetric likelihood
model used during fitting.

## Randomness, Threads, and Reproducibility

One segmentation request supplies a master seed. Baysor derives versioned,
named streams for:

- core scientific operations;
- molecule clustering;
- neighborhood-composition output; and
- diagnostics.

Separating these streams prevents enabling a diagnostic-only product from
changing scientific assignments by consuming random values intended for the
core algorithm.

Deterministic semantic replay is guaranteed only when the prepared input,
options, build, platform, master seed, and explicit one-thread execution are the
same. Parallel runs are not guaranteed to produce identical results, even with
the same seed and thread count, because dynamic scheduling and floating-point
evaluation order can vary.

The authoritative guarantee and the meaning of semantic equality are defined in
the [Native Segmentation API Contract](../contract.md#randomness-and-native-threads).

## Why the Input Context Matters

Baysor does not treat molecules as independent rows. Several stages depend on
the surrounding molecule cloud:

- automatic parameter resolution uses molecule and prior counts;
- confidence estimation uses nearest-neighbor distance distributions;
- the assignment candidate set and MRF evidence come from the spatial graph;
- optional clustering uses local neighborhood composition; and
- boundaries use the final collection of coordinates assigned to each cell.

Consequently, running the algorithm on a spatial subset can differ from running
it on the complete field, especially near subset edges or when the subset is too
small to represent the dataset's density and expression composition. Any tiled
execution strategy must account for that context through overlap, globally
resolved calibration where appropriate, and post-segmentation reconciliation.

## Implementation Map

| Stage | Principal implementation |
| --- | --- |
| Native orchestration and result materialization | [`src/segmentation/segmentation.cpp`](../src/segmentation/segmentation.cpp) |
| Input loading | [`src/data_loading`](../src/data_loading) |
| Option defaults and validation | [`src/utils/options.cpp`](../src/utils/options.cpp) |
| Prior loading and scale estimation | [`src/data_loading/prior_segmentation.cpp`](../src/data_loading/prior_segmentation.cpp) |
| Confidence and noise estimation | [`src/processing/data_processing/noise_estimation.cpp`](../src/processing/data_processing/noise_estimation.cpp) |
| Graph construction and BMM initialization | [`src/processing/data_processing/initialization.cpp`](../src/processing/data_processing/initialization.cpp) |
| ICA/MRF molecule clustering | [`src/processing/bmm_algorithm/molecule_clustering.cpp`](../src/processing/bmm_algorithm/molecule_clustering.cpp) |
| Louvain/Leiden molecule clustering | [`src/processing/bmm_algorithm/molecule_clustering_louvain.cpp`](../src/processing/bmm_algorithm/molecule_clustering_louvain.cpp) |
| Iterative BMM | [`src/processing/bmm_algorithm/bmm_algorithm.cpp`](../src/processing/bmm_algorithm/bmm_algorithm.cpp) |
| Cell component model | [`include/baysor/processing/models/component.h`](../include/baysor/processing/models/component.h) |
| Boundary estimation | [`src/processing/data_processing/boundary_estimation.cpp`](../src/processing/data_processing/boundary_estimation.cpp) |
| Native serializers | [`src/reporting/output.cpp`](../src/reporting/output.cpp) |

## Related Documentation

- [Running Baysor](run.md)
- [Configuration](configuration.md)
- [Prior Segmentation Inputs](priors.md)
- [Outputs](outputs.md)
- [Native Segmentation API Contract](../contract.md)
- [Citation](citation.md)
