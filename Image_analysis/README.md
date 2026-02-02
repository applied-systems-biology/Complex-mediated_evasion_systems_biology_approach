## Image analysis

The image analysis was performed in **JIPipe v5.3.0**: https://www.jipipe.org/

This repository contains the following analysis files:

### JIPipe workflow
- **`LL37 analysis.jip`**  
  Main analysis file that should be loaded in JIPipe.

### Cellpose model (retrained)
- **`calb_LL37_CP_preprocessed.zip`** and **`calb_LL37_CP_preprocessed.zip.json`**  
  Retrained **CellPose** network for segmenting *C. albicans* yeast cells, based on the pretrained **`cyto3`** network.  
  This is loaded in the `LL37 analysis.jip` workflow.

### Manual segmentations (training / test data)
- **`s1_F50.tif`**  
  Manual segmentation of frame 50 from `20241106_LL-37_s1.ome.tiff`. Used as test data during retraining of the CellPose network in `LL37 analysis.jip`.
- **`s2_F50.tif`**  
  Manual segmentation of frame 50 from `20241106_LL-37_s2.ome.tiff`. Used for retraining of the CellPose network in `LL37 analysis.jip`.

### Downstream tracking analysis (R)
- **`LL37 cell death from tracks.Rmd`**  
  R Markdown file for analysis of the track files to count alive and dead cells from the time-lapse data.
