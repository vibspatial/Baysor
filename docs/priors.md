# Prior Segmentation Inputs

Baysor can use optional prior segmentation to guide the run.

The `prior_segmentation` positional argument may be:

- a transcript-native label column: `:column_name`
- an image mask
- a boundary CSV / Parquet file

## Transcript-Native Priors

If prior labels are already attached to each molecule row, pass them as:

```bash
./build/baysor run molecules.parquet :cell_id
```

This is often the fastest and simplest prior mode, because no extra
point-in-polygon or mask lookup is required.

For Xenium, `:cell_id` is usually the preferred prior mode.

## Image Mask Priors

Image masks are useful when segmentation is available as labeled pixels, such
as DAPI-based or watershed-derived masks.

Example:

```bash
./build/baysor run -c configs/iss.toml molecules.csv dapi_mask.tif
```

## Boundary Priors

Polygon / boundary priors are useful when segmentation is published as cell
outlines rather than masks.

Examples:

```bash
./build/baysor run -c configs/xenium.toml data/experiment.xenium data/cell_boundaries.parquet
./build/baysor run -c configs/xenium.toml data/experiment.xenium data/nucleus_boundaries.parquet
```

## Prior Confidence

The weight of the prior is controlled by:

```text
--prior-segmentation-confidence
```

and by the corresponding config value:

```toml
[segmentation]
prior_segmentation_confidence = ...
```

Higher values make the segmentation adhere more strongly to the prior.

## Unassigned Prior Labels

When using transcript-native priors, the unassigned label can be set with:

```text
--unassigned-prior-label
```

For Xenium, the default unassigned label in [configs/xenium.toml](../configs/xenium.toml)
is:

```toml
[prior]
unassigned_label = "UNASSIGNED"
```

## When To Use Which Prior

Use `:column_name` when:

- prior assignments already exist per transcript
- you want the fastest prior mode
- you want clean Xenium Ranger handoff later

Use boundary priors when:

- you want Baysor to follow a published polygon segmentation more explicitly
- you have cell or nucleus boundaries but no transcript-native labels

Use no prior when:

- you want a fully de novo segmentation
- you have no trustworthy prior
- you are comfortable setting `--scale` explicitly
