Use NeMO Analytics Single Cell Workbench (https://nemoanalytics.org/sc_workbench.html) to the 8th step (clustering) and download the H5AD file. 
This file will contain all of the paramiters entered into NeMO including dimentions and PCA. It also contains all genes and their Ensembl ID.

Run the inspect_h5ad.py file via the command center: 
open the py file and change the file path to the H5AD file path (shift right-click on file and 'copy as path') paste into the "" leaving the r. Do not have multiple "" sets.

For my set up the commad window requires the input of 'python C:\code\NEMO\inspect_h5ad.py'
Once the command center/window displays save sucessful, the exported .mat can be renamed and moved. This python code only needs to be run once per H5AD file.

Run the matlab code using the new .mat file path:
The nemo_umap_gene_browser('filepath') runs a GUI app similar to NeMO's website.
The nemo_umap_multigene_panels('filepath', "gene1,gene2,gene3") runs an export to tiff version of the browser (there is no limit to the number of enterable genes in panels but they must be in the "x,y,z" format).

.mat files must be named with the following convention "region sub-region descriptor.mat" the region is a two letter code corresponding to: thalamus "TH", Cortex "CX", and Ganglion Eminance "GE".
These must be present to ensure appropriate selection of marker genes and cluster annotations.
The following 4 letters (including the space after region) should correspond to the sub-region in order to ensure file and figure titles are named correctly.

Example: 'CX Prefrontal GW18 Matlab.mat' uses CX for annotations and will generate file/figures titled 'CX Pre_' followed by the 
figure title 'AtlasClusterMap' or 'MultiUMAP'.
This same nameing convention will be seen in the browser titles.

If adding or changing clustering annotations use format: ' "annotation name", ["gene1", "gene2"]; '
Ensure changes are made within the appropriate markerDB{ } functions. 
Add this under the Markers secion (beginning line 52).

To add more regions add a new region map where ' regionMap('region abbreviation')="region name"; '
Add this under the Region Map secion (beginning line 38).

Note Not All Genes Will Be Found, If Not Found Try An Alias (ie. COUTF2 = NR2F2, NKX2.1 = NKX2-1).
