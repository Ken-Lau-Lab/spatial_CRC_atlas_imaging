# Step 6

`export_IES_images_MxIF_size.ipynb` - Uses MILWRM ST image registration to create MILWRM `.npz` images of immune exclusion signature (IES) scores from ST in MxIF pixel space following registration to ST in [`../step5/MxIF_visium_transformation.ipynb`](../step5/MxIF_visium_transformation.ipynb)

`ST_wholeslide_transform.ipynb` - Python notebook for combining multiple Visium ST `.h5ad`s into a single, whole-slide `.h5ad` using the inverse affine transformations from [`../step4/MxIF_visium_registration.ipynb`](../step4/MxIF_visium_registration.ipynb) and whole-slide virtual H&E images from MxIF
