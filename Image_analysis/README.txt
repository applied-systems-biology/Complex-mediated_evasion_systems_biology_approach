The image analysis is done in JIPipe (version 5.3.0, https://www.jipipe.org/). The image analysis consist of the following files

LL37 analysis.jip: The main analysis file that should be loaded by JIPipe
calb_LL37_CP_preprocessed.zip/calb_LL37_CP_preprocessed.zip.json: The retrained CellPose network for the segmentation of C. albicans yeast cells based on the pretrained 'cyto3' network. This is loaded by the in the LL37 analysis.jip workflow
s1_F50.tif: Manual segmentation of the 50th frame of the data 20241106_LL-37_s1.ome.tiff. Used as test data during  retraining of the CellPose network in LL37 analysis.jip. 
s2_F50.tif: Manual segmentation of the 50th frame of the data 20241106_LL-37_s2.ome.tiff. Used for retraining of the CellPose network in LL37 analysis.jip.
LL37 cell death from tracks.Rmd: R-Markdown file for analysis of the track files to count alive and dead cells from the time lapse data.