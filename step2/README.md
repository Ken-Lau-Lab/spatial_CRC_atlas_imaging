# Step 2

`mesmer_to_MIRIAM.m` - MATLAB script for running MIRIAM segmentation and quantification on mesmer outputs from [`../step1/apply_mesmer_local.ipynb`](../step1/apply_mesmer_local.ipynb). **NOTE:** In order to run this MATLAB script, you need to install the [Statistics and Machine Learning Toolbox](https://www.mathworks.com/products/statistics.html), and add the source code in [`../resources/MIRIAM_functions/`](../resources/MIRIAM_functions/) to MATLAB's search path.

`tile_mxif.ipynb` - Python notebook for "tiling" or stitching together MxIF regions into whole-slide space using the "canvas" from [`../step1/compile_region_canvas.ipynb`](../step1/compile_region_canvas.ipynb) as a guide
