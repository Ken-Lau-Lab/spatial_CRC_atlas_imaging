# spatial_CRC_atlas_imaging

Code to reproduce results from [Heiser, *et al.* (2023)](https://doi.org/10.1101/2023.03.09.530832) that employ multiplex immunofluorescence (MxIF) imaging data and their integration with Visium spatial transcriptomics (ST).

## Purpose

These analyses integrate multiplex immunofluorescence (MxIF) and spatial transcriptomics (ST) to generate models and metrics presented in our manuscript, "Molecular cartography uncovers evolutionary and microenvironmental dynamics in sporadic colorectal tumors" (currently in review at *Cell*; bioRxiv preprint available: [doi.org/10.1101/2023.03.09.530832](https://doi.org/10.1101/2023.03.09.530832)). Functions herein include:

* Reading in and pre-processing data
* Single-cell segmentation and quantification of MxIF images
* Spatial registration of MxIF and ST
* Summaries and visualizations of integrated spatial data
* Immune Exclusion analysis in MxIF

## Order of Operations

Code in this repository is intended to be run in order, as some analyses depend on outputs from previous tools. We therefore divide up the analysis notebooks into steps, provided as ordered directories.

1. Download publicly available data using the scripts in the `data/` directory first
2. Proceed to the notebooks in the `step1/` directory, then `step2/`, and `step3/` and so forth. All notebooks residing in the same directory can be run in any order, but they *all* must be successfully run before proceeding to the next step.
