# Step 4

`PosStats_tiling.ipynb` - Python notebook for "tiling" out segmentation results ("PosStats") from individual image tile coordinates into whole-slide space

`MxIF_visium_registration.ipynb` - Python notebook for spatially registering and cropping whole-slide MxIF MILWRM `.npz`s from `../step3/MxIF_maskgen_imgcreation_#.ipynb` to Visium ST space using Justin Shao's `napari` functions. Results are affine transformations matrices used to transform images downstream. *Note: `napari` plugin for registration code is in development at [github.com/codyheiser/napari-registration](github.com/codyheiser/napari-registration).*

`overlays.ipynb` - Python notebook for plotting whole-slide MxIF image overlays using MILWRM functionality
