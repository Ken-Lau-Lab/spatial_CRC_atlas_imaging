import os
import pandas as pd
from math import ceil
from skimage.io import imread
from PIL import Image

Image.MAX_IMAGE_PIXELS = 933120000

mychannels = [
    "VIMENTIN",
    "SOX9",
    "SMA",
    "PSTAT3",
    "PEGFR",
    "PCNA",
    "PANCK",
    "OLFM4",
    "NAKATPASE",
    "MUC5AC",
    "MUC2",
    "LYSOZYME",
    "HLAA",
    "GAMMAACTIN",
    "FOXP3",
    "ERBB2",
    "COLLAGEN",
    "CGA",
    "CDX2",
    "CD68",
    "CD45",
    "CD20",
    "CD11B",
    "CD8",
    "CD4_",
    "CD3D",
    "BCATENIN",
]


def tile_tiffs(
    tiffdir,
    canvas,
    name,
    regions,
    channels=mychannels,
    common_strings=None,
    RGB=False,
    scale=8,
):
    """
    Parameters
    ----------
    tiffdir : str
        Path to directory containing `.tif` files for a multiplexed image
    canvas : str
        Path to "region_canvas.csv" file
    name : str
        Slide identifier (e.g. "WD868055")
    regions : list of int
        List of regions that make up the whole image for slide with identifier `name`.
        e.g. `[1,2,3,4,5]` for a slide with 5 regions.
    channels : tuple of str, optional (default=`mychannels` defined above)
        List of channels present in `.tif` file names (case-sensitive)
        corresponding to img.shape[2] e.g. `("ACTG1","BCATENIN","DAPI",...)`
    common_strings : str, list of str, or `None`, optional (default=None)
        Strings to look for in all `.tif` files in `tiffdir` corresponding to
        `channels` e.g. `("WD86055_", "_region_001.tif")` for files named
        "WD86055_[MARKERNAME]_region_001.tif". If `None`, assume that only 1 image
        for each marker in `channels` is present in `tiffdir`.
    RGB : bool, optional (default=`False`)
        Are the images 3-channel RGB files, or single-channel grayscale? Only set to
        `True` if you're stitching virtual H&E files.
    scale : float (optional, default=8)
        Factor by which to downsample image regions before stitching

    Returns
    -------
    Writes downsampled and stitched images for each channel to
    "{tiffdir}/{channel}_stitched_downsample{scale}.tif"
    """
    c = pd.read_csv(canvas)  # read in region_canvas file
    if common_strings is not None:
        # coerce single string to list
        if isinstance(common_strings, str):
            common_strings = [common_strings]
    A = []  # list for dumping numpy arrays
    # loop through channels, stitch images and save to `tiffdir``
    for channel in channels:
        try:
            print(channel)
            A = []  # list for dumping numpy arrays
            # loop through regions and read images into list
            for region in regions:
                if common_strings is None:
                    # find file matching all common_strings and channel name
                    f = [
                        f
                        for f in os.listdir(tiffdir)
                        if all(x in f for x in [channel, "region_00{}.tif".format(region)])
                    ]
                    print(f)
                else:
                    # find file matching all common_strings and channel name
                    f = [
                        f
                        for f in os.listdir(tiffdir)
                        if all(
                            x in f
                            for x in common_strings
                            + [channel, "region_00{}".format(region)]
                        )
                    ]
                # assertions so we only get one file per channel
                assert len(f) != 0, "No file found with channel {}".format(channel)
                assert (
                    len(f) == 1
                ), "More than one match found for file with channel {}".format(channel)
                f = os.path.join(tiffdir, f[0])  # get full path to file for reading
                print("Reading marker {} from {}".format(channel, f))
                if RGB:
                    # read in 3-channel image (histology or VHE)
                    #tmp = Image.open(f)
                    tmp = imread(f)
                    tmp = Image.fromarray(tmp, mode="RGB")  # read in .tif file
                else:
                    # read in single channel image
                    tmp = imread(f)
                    tmp = Image.fromarray(tmp, mode="I;16")  # read in .tif file
                A.append(tmp)  # append numpy array to list
            # calculate extent of x & y directions
            x_adj = (
                abs(c.loc[c.SlideID == name, "du"])
                - abs(c.loc[c.SlideID == name, "du"]).min()
            )
            y_adj = (
                abs(c.loc[c.SlideID == name, "dv"])
                - abs(c.loc[c.SlideID == name, "dv"]).min()
            )
            x = (x_adj + c.loc[c.SlideID == name, "totalWidth"]).max()
            y = (y_adj + c.loc[c.SlideID == name, "totalHeight"]).max()
            # make new image
            outshape = (ceil(x), ceil(y))
            if RGB:
                # create new image size of stitched whole-slide image (histology or VHE)
                output = Image.new("RGB", outshape)
            else:
                # create new image size of stitched whole-slide image
                output = Image.new("I;16", outshape)
            # loop through region tiles, get their new x and y origins, paste into `output`
            for i, j, tile in zip(
                abs(c.loc[c.SlideID == name, "du"])
                - abs(c.loc[c.SlideID == name, "du"]).min(),
                abs(c.loc[c.SlideID == name, "dv"])
                - abs(c.loc[c.SlideID == name, "dv"]).min(),
                A,
            ):
                # paste region tile into `output`
                # this part is sketchy because of rounding i and j... 
                # not sure of another way to do it to avoid fractional pixels.
                #print("Tile shape: {}\t{}, {}".format(tile.size, ceil(i), ceil(j)))
                output.paste(tile, (ceil(i), ceil(j)))
            # downsample to manageable size and save image
            output_small = output.resize(
                (outshape[0] // scale, outshape[1] // scale),
                resample=Image.NEAREST,
            )
            output_small.save(
                "{}/{}_stitched_downsample{}.tif".format(tiffdir, channel, scale)
            )
        except AssertionError as err:
            print("{} - moving to next marker".format(err))
