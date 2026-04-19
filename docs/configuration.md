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
- `estimate_scale_from_centers`
- `n_clusters`
- `prior_segmentation_confidence`
- `iters`
- `n_cells_init`
- `nuclei_genes`
- `cyto_genes`
- `unassigned_prior_label`

### `[plotting]`

- `gene_composition_neigborhood`
- `min_pixels_per_cell`
- `max_plot_size`
- `ncv_method`

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
