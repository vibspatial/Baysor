# Configuration

The main protocol presets live under [configs/](../configs/):

- [example_config.toml](../configs/example_config.toml)
- [xenium.toml](../configs/xenium.toml)
- [iss.toml](../configs/iss.toml)
- [starmap.toml](../configs/starmap.toml)
- [osm_fish.toml](../configs/osm_fish.toml)

Use them with:

```bash
./build/baysor run -c configs/xenium.toml ...
```

## Config Structure

The example config is organized into sections such as:

- `[data]`
- `[segmentation]`
- `[plotting]`

Common fields include:

### `[data]`

- column names (`x`, `y`, `z`, `gene`)
- filtering:
  - `min_molecules_per_gene`
  - `exclude_genes`
  - `min_molecules_per_cell`
- `force_2d`

### `[segmentation]`

- `scale`
- `scale_std`
- `cluster_method`
- `estimate_scale_from_centers`
- `n_clusters`
- `cluster_resolution`
- `cluster_graph_k`
- `cluster_n_dims`
- `cluster_basis_sample_size`
- `prior_segmentation_confidence`
- `iters`
- `n_cells_init`
- `nuclei_genes`
- `cyto_genes`
- `unassigned_prior_label` (accepted for compatibility; prefer
  `[prior].unassigned_label` in new configs)

### `[prior]`

- `type`
- `path`
- `column_name`
- `unassigned_label`
- `min_molecules_per_segment`
- `estimate_scale_from_prior`

### `[plotting]`

- `gene_composition_neigborhood`
- `min_pixels_per_cell`
- `max_plot_size`
- `ncv_method`

## Clustering Options

The segmentation pre-clustering step can be configured through
`[segmentation]`.

- `cluster_method = "mrf" | "louvain" | "leiden" | "none"`
  - `mrf` is the legacy Baysor molecule clustering path and remains the default.
  - `louvain` and `leiden` cluster basis anchors in NCV space, then transfer the
    final labels back to all molecules.
  - `none` disables the clustering prior entirely.
- `n_clusters`
  - for `mrf`, this is the exact number of molecule clusters used by the
    clustering model
  - for `louvain` and `leiden`, this is the target final number of coarse
    clusters after anchor communities are merged
  - if omitted or set to `0`, Baysor uses:
    - `4` for `mrf`
    - `10` for `louvain` and `leiden`
- `cluster_resolution`
  - advanced parameter for `louvain` and `leiden`
  - controls the initial overclustering on the anchor graph before merging back
    to `n_clusters`
- `cluster_graph_k`
  - number of nearest neighbors used in NCV space for:
    - the Louvain/Leiden anchor graph
    - 2D and 3D NCV UMAPs
    - NCV interpolation
  - the default is `15` and usually does not need to be changed
- `cluster_n_dims`
  - number of low-dimensional NCV coordinates used by graph clustering
- `cluster_basis_sample_size`
  - maximum number of basis anchors used by graph clustering and shared NCV
    embedding

For very large runs, especially high-gene-panel Xenium runs such as 5K panels,
the recommended clustering prior is:

```toml
[segmentation]
cluster_method = "louvain"
n_clusters = 10
```

Treat `10` as a practical starting point for the final coarse cluster count.
Increase or decrease it only if the diagnostic report suggests the expression
structure is clearly under- or over-clustered.

## NCV Neighborhoods

Two different neighborhood sizes matter for NCV-based reporting and
graph clustering:

- `gene_composition_neigborhood` / `gene_composition_neighborhood`
  - spatial neighborhood size used to compute the NCV itself
  - if not set, Baysor falls back to
    `max(n_genes / 10, min_molecules_per_cell, 3)`
- `cluster_graph_k`
  - neighbor count in the low-dimensional NCV space, after the NCV has already
    been computed

Both spellings of `gene_composition_neigborhood` are accepted. The misspelled
form is kept for compatibility with existing Julia-era configs.

## Prior Options

Use `[prior]` for prior-input-specific settings in new configs:

```toml
[prior]
unassigned_label = "UNASSIGNED"
estimate_scale_from_prior = true
```

`[segmentation].unassigned_prior_label` and
`[segmentation].estimate_scale_from_centers` are still accepted for
compatibility with older configs. CLI flags such as
`--unassigned-prior-label` override both config forms.

## CLI vs Config

Config values set defaults. CLI flags override them.

Typical pattern:

```bash
./build/baysor run \
  -c configs/xenium.toml \
  --x-min 0 --x-max 2000 \
  --y-min 0 --y-max 2000 \
  data/experiment.xenium :cell_id
```

## Protocol-Specific Notes

### Xenium

[configs/xenium.toml](../configs/xenium.toml) maps Xenium columns:

- `x_location`
- `y_location`
- `z_location`
- `feature_name`
- `qv`

and sets:

- Xenium-style transcript filtering
- Xenium prior defaults

### ISS / STARmap / osm-FISH

The protocol config files mainly adjust:

- minimum molecules per gene / cell
- plotting defaults
- expected coordinate / gene column layout

## Recommendation

For protocol-specific work, start from the closest existing config and override
only the values that are truly dataset-specific. For a fuller list of fields and
comments, use [example_config.toml](../configs/example_config.toml).
