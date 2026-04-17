# Example run on STARmap data

Data is published in [Wang X, et al. - Three-dimensional intact-tissue sequencing of single-cell transcriptional states](https://doi.org/10.1126/science.aat5691) and available at [the website](https://www.starmapresources.com/data/).


## Get the data

Create directories:

```bash
mkdir -p data output_dapi output_no_dapi
```

### Preprocessing raw data

In this example we used "visual_1020/20180410-BY3_1kgenes" dataset. To convert it to a proper format you need to download files "goodPoints.mat", "labels.npz" and "cell_barcode_names.csv" (which, indeed, contains gene barcodes) to the "data" folder. From there you need to run the converting script:

```bash
python3 ../convert_to_csv.py goodPoints.mat cell_barcode_names.csv labels.npz
```

*OR*
### Download pre-processed data

```bash
wget -P data http://pklab.med.harvard.edu/viktor/baysor/starmap/molecules.csv
wget -P data http://pklab.med.harvard.edu/viktor/baysor/starmap/segmentation.tiff
```

## CLI run

### With prior segmentation

```bash
baysor run --iters 500 --plot -c ../../configs/starmap.toml -o ./output_dapi ./data/molecules.csv ./data/segmentation.tiff
```

### Without prior segmentation

```bash
baysor run --iters 500 --plot --scale 89.995 --scale-std 16.398 -c ../../configs/starmap.toml -o ./output_no_dapi ./data/molecules.csv
```
