from MILWRM.ST import *
from MILWRM.MxIF import *
from skimage import transform, measure, io
from skimage.morphology import skeletonize, binary_dilation
from pylab import *
from napari.utils.transforms import Affine


def flip_coords(coords, vert=True, horiz=True):
    """
    Invert coordinates in cartesian (x,y) layout

    Parameters
    ----------
    coords : np.array
        x, y coordinates
    vert : bool
        Flip vertically
    horiz : bool
        Flip horizontally
    """
    if vert:
        mx = coords[:, 1].max()
        mn = coords[:, 1].min()
        coords[:, 1] = (mx - coords[:, 1]) + mn
    if horiz:
        mx = coords[:, 0].max()
        mn = coords[:, 0].min()
        coords[:, 0] = (mx - coords[:, 0]) + mn
    return coords


def transform_coords(
    coords,
    affine,
    flip=False,
):
    """
    Transform coordinates using affine

    Parameters
    ----------
    coords : np.array
        x, y coordinates
    affine : np.array
        Affine transformation coordinates
    shape : np.array
        Shape of final image
    flip : bool
        Flip image horizontally before transforming
    """
    if flip:
        coords_new = transform(
            src=flip_coords(coords, horiz=True, vert=False), m=affine
        )
    else:
        coords_new = transform(src=coords, m=affine)
    return coords_new


def transform_image(
    image,
    affine,
    shape,
    flip=False,
    plot=True,
    dilation_iterations=5,
):
    """
    Transform image using affine

    Parameters
    ----------
    image : np.array
        The image to transform
    affine : np.array
        Affine transformation coordinates
    shape : np.array
        Shape of final image
    flip : bool
        Flip image horizontally before transforming
    plot : bool
        imshow a dilated version of output for visualization
    """
    if flip:
        if image.ndim == 2:
            image = transform.warp(
                image[::-1, ::-1], transform.AffineTransform(affine), output_shape=shape
            )
        elif image.ndim == 3:
            image = transform.warp(
                image[::-1, ::-1, :],
                transform.AffineTransform(affine),
                output_shape=shape,
            )
    else:
        image = transform.warp(
            image, transform.AffineTransform(affine), output_shape=shape
        )
    if plot:
        kernel = np.ones((5, 5), np.uint8)
        plt.imshow(dilate(image, kernel=kernel, iterations=dilation_iterations))
    return image


def reconcile_image_shape(ref_img, alter_img):
    print(
        "Images differ by {} pixel in x and {} pixels in y.".format(
            ref_img.shape[0] - alter_img.shape[0],
            ref_img.shape[1] - alter_img.shape[1],
        )
    )
    if (ref_img.shape[0] - alter_img.shape[0]) > 0:
        alter_img = np.pad(
            alter_img, ((0, ref_img.shape[0] - alter_img.shape[0]), (0, 0))
        )
        print(
            "Padding x by {} - new shape: {}".format(
                ref_img.shape[0] - alter_img.shape[0], alter_img.shape
            )
        )
    elif (ref_img.shape[0] - alter_img.shape[0]) < 0:
        alter_img = alter_img[
            : ref_img.shape[0],
        ].copy()
        print(
            "Removed {} pixels from x - new shape: {}".format(
                abs(ref_img.shape[0] - alter_img.shape[0]), alter_img.shape
            )
        )
    if (ref_img.shape[1] - alter_img.shape[1]) > 0:
        alter_img = np.pad(
            alter_img, ((0, 0), (0, ref_img.shape[1] - alter_img.shape[1]))
        )
        print(
            "Padding y by {} - new shape: {}".format(
                ref_img.shape[1] - alter_img.shape[1], alter_img.shape
            )
        )
    elif (ref_img.shape[1] - alter_img.shape[1]) < 0:
        alter_img = alter_img[:, : ref_img.shape[1]].copy()
        print(
            "Removed {} pixels from y - new shape: {}".format(
                abs(ref_img.shape[1] - alter_img.shape[1]), alter_img.shape
            )
        )
    return alter_img


def tile_coords(
    coords,
    canvas,
    name,
    regions,
    scale=8,
):
    """
    Transform x-y cartesian coordinates to whole-slide (ws) space

    Parameters
    ----------
    coords : pd.DataFrame
        Observations dataframe from Anndata containing x-y cartesian coordinates in pixel space
    canvas : str
        Path to "region_canvas.csv" file
    name : str
        Slide identifier (e.g. "WD868055")
    regions : list of int
        List of regions that make up the whole image for slide with identifier `name`.
        e.g. `["region_001","region_002","region_003"]` for a slide with 3 regions.
    scale : float (optional, default=8)
        Factor by which to downsample image regions before stitching

    Returns
    -------
    coords : np.array
        Adjusted x-y cartesian coordinates in whole-slide pixel space
    """
    c = pd.read_csv(canvas)  # read in region_canvas file
    # calculate extent of x & y directions
    x_adj = (
        abs(c.loc[c.SlideID == name, "du"]) - abs(c.loc[c.SlideID == name, "du"]).min()
    )
    y_adj = (
        abs(c.loc[c.SlideID == name, "dv"]) - abs(c.loc[c.SlideID == name, "dv"]).min()
    )
    x = (x_adj + c.loc[c.SlideID == name, "totalWidth"]).max()
    y = (y_adj + c.loc[c.SlideID == name, "totalHeight"]).max()

    outshape = (ceil(x), ceil(y))

    # loop through region tiles, get their new x and y origins, paste into `output`
    coords["x_ws"] = 0
    coords["y_ws"] = 0
    for i, j, tile in zip(
        abs(c.loc[c.SlideID == name, "du"]) - abs(c.loc[c.SlideID == name, "du"]).min(),
        abs(c.loc[c.SlideID == name, "dv"]) - abs(c.loc[c.SlideID == name, "dv"]).min(),
        regions,
    ):
        coords.loc[coords.region == tile, "x_ws"] = (
            coords.loc[coords.region == tile, "x"] + i
        )
        coords.loc[coords.region == tile, "y_ws"] = (
            coords.loc[coords.region == tile, "y"] + j
        )

    coords["x_ws"] = coords["x_ws"] / scale
    coords["y_ws"] = coords["y_ws"] / scale
    return coords


def identify_cells(adata, celldict, thresh=0.4):
    """
    Make cell type ID calls from positive marker expression

    Parameters
    ----------
    adata : anndata.AnnData
        Object with thresholded marker values in .obs
    celldict : dict
        Dictionary with cell types as keys and lists of marker names as values.
        Order of keys matters, as cells will be called sequentially. Therefore,
        if one cell type shares markers with another, place the cell type with
        more markers later in the dict.

    Returns
    -------
    adata is edited in place, adding "Cell Type" to .obs
    """
    adata.obs["Cell Type"] = ""
    for celltype in celldict.keys():
        print(celltype)
        for marker in celldict[celltype]:
            if "{}_pos".format(marker) not in adata.obs.columns:
                print("\tThresholding {}".format(marker))
                adata.obs["{}_pos".format(marker)] = 0
                adata.obs.loc[
                    adata.obs["Area_Thresh_{}".format(marker)] / adata.obs.Cell_Area
                    > thresh,
                    "{}_pos".format(marker),
                ] = 1
        if len(celldict[celltype]) == 1:
            adata.obs.loc[
                adata.obs["{}_pos".format(celldict[celltype][0])] == 1, "Cell Type"
            ] = celltype
        if len(celldict[celltype]) == 2:
            adata.obs.loc[
                (adata.obs["{}_pos".format(celldict[celltype][0])] == 1)
                & (adata.obs["{}_pos".format(celldict[celltype][1])] == 1),
                "Cell Type",
            ] = celltype
        if len(celldict[celltype]) == 3:
            adata.obs.loc[
                (adata.obs["{}_pos".format(celldict[celltype][0])] == 1)
                & (adata.obs["{}_pos".format(celldict[celltype][1])] == 1)
                & (adata.obs["{}_pos".format(celldict[celltype][2])] == 1),
                "Cell Type",
            ] = celltype


def process_segmentation(
    input_csv, output_h5ad, tile_kwargs, id_cells_kwargs, plot_out=True, spot_size=50
):
    # read segmentation quant from .csv
    print("Reading from {}".format(input_csv))
    a = pd.read_csv(input_csv, index_col=0)
    a.index = "cell_" + a.index.astype(str)  # avoid numerical indices for AnnData
    # make AnnData object
    numeric = a.columns[a.columns.str.startswith("Area_Thresh_")]
    obs = a
    a = sc.AnnData(a.loc[:, numeric])
    a.obs = obs
    # get coords to ws space
    a.obs = tile_coords(a.obs, **tile_kwargs)
    # add coords to .obsm
    a.obsm["spatial_mm"] = a.obs[["x_mm", "y_mm"]].values
    a.obsm["spatial_ws"] = a.obs[["x_ws", "y_ws"]].values
    # ID cells in adata
    print("Making cell type calls with identify_cells function")
    identify_cells(a, **id_cells_kwargs)
    # subset to only cells ID'ed above
    a_sub = a[a.obs["Cell Type"] != "", :].copy()
    if plot_out:
        sc.pl.spatial(
            a_sub,
            basis="spatial_ws",
            color=["Cell Type"],
            ncols=4,
            img=None,
            spot_size=spot_size,
            scale_factor=1,
            frameon=False,
            cmap="viridis",
        )
    # write h5ad file
    print("Writing AnnData output to {}".format(output_h5ad))
    a_sub.write(output_h5ad, compression="gzip")


def quantify_cells_per_clone(
    visium_adata,
    seg_adata,
    affine,
    shape,
    img_obj,
    visium_pita_save=None,
    seg_pita_save=None,
    flip=False,
):
    master_dict = {}
    shape = shape.astype(int)
    # assemble pita
    celltype_pita = None
    p = assemble_pita(
        visium_adata,
        features="CNV Clone",
        use_rep="obs",
        plot_out=True,
        label="CNV Clone",
        save_to=visium_pita_save,
    )
    for celltype in seg_adata.obs["Cell Type"].cat.categories:
        print(celltype)
        # translate cell centroids back into an image
        a_celltype = np.zeros(img_obj.img.shape[:2])
        errcounter = 0
        for x, y in ceil(
            seg_adata[seg_adata.obs["Cell Type"] == celltype, :].obsm["spatial_ws"]
        ):
            try:
                a_celltype[int(y), int(x)] = 1
            except IndexError:
                errcounter += 1
        if errcounter > 0:
            print(
                "{} {} cells skipped because they were outside the Visium area".format(
                    errcounter, celltype
                )
            )
        # transform whole-slide image to Visium space
        a_celltype = transform_image(
            a_celltype,
            affine,
            shape,
            flip=flip,
            plot=False,
        )
        # downsample
        a_celltype = measure.block_reduce(a_celltype, block_size=2, func=np.max)
        # get images same shape
        a_celltype = reconcile_image_shape(p[0], a_celltype)
        # bring fractional values up to 1.0
        a_celltype[np.where(a_celltype != 0)] = 1
        # mask cell image with visium pita
        a_celltype[np.where(np.isnan(p[0]))] = 0
        # get single pixels per cell back
        a_celltype = skeletonize(a_celltype)
        print("Number of centroids", end=" - ")
        print(np.unique(a_celltype, return_counts=True))

        # mask and count values
        p2 = p[0].copy()  # deep copy np array image
        p2[np.where(a_celltype == 0)] = np.nan
        counts = dict(
            zip(
                np.unique(p2, return_counts=True)[0],
                np.unique(p2, return_counts=True)[1],
            )
        )
        # reconfigure dict with interpretable keys
        new_counts = {}
        for key in counts.keys():
            if not np.isnan(key):
                print(key)
                new_counts[p[1][0][1][key]] = counts[key]
        # add to master dict
        master_dict[celltype] = new_counts

        # add to celltype_pita for plotting later
        kernel = np.ones((5, 5), np.uint8)
        if celltype_pita is None:
            celltype_pita = binary_dilation(a_celltype, footprint=kernel)
        else:
            celltype_pita = np.dstack(
                [celltype_pita, binary_dilation(a_celltype, footprint=kernel)]
            )

    # plot celltype_pita
    _ = show_pita(
        celltype_pita,
        label=list(seg_adata.obs["Cell Type"].cat.categories),
        ncols=3,
        save_to=seg_pita_save,
    )
    return master_dict
