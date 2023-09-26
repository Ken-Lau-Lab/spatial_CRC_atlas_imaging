import napari
import numpy as np
from napari.utils.transforms import Affine


def registration_viewer(bg_img, fg_img):
    """
    Run napari registration plugin

    Parameters
    ----------
    bg_img : np.array
        Background image to register `fg_img` to
    fg_img : np.array
        Foreground image to register to `bg_img`
    """
    viewer = napari.view_image(bg_img)  # whole-slide image to register fg_img onto
    viewer.add_image(fg_img, opacity=0.5)  # fg_img on top with transparency
    # set up interaction box for moving/scaling fg_img
    viewer.layers.selection.active.interactive = False
    viewer.overlays.interaction_box.points = viewer.layers.selection.active.extent.world
    viewer.overlays.interaction_box.show = True
    viewer.overlays.interaction_box.show_vertices = True
    viewer.overlays.interaction_box.show_handle = True
    viewer.overlays.interaction_box.allow_new_selection = False
    # define inner function for capturing affine event
    def on_transform_changed_drag(event):
        """record affine transform as ST is manipulated"""
        viewer.layers.selection.active.affine = event.value
    viewer.overlays.interaction_box.events.transform_drag.connect(on_transform_changed_drag)
    # open viewer
    napari.run()
    # return viewer object for capturing live information
    return viewer


def extract_affine(viewer, fg_img):
    """
    Extract the inverse of the applied affine transformation from napari `registration_viewer`

    Parameters
    ----------
    viewer : napari.viewer
        Napari viewer object initialized by `registration_viewer`
    fg_img : np.array
        Foreground image to register to `bg_img`

    Returns
    -------
    affine : np.array
        Affine matrix describing transformation
    scale : tuple
        Scale from affine transformation
    shape : tuple
        New shape of transformed image
    """
    # extract inverse of affine
    cut = Affine(
        rotate=viewer.layers.selection.active.affine.inverse.rotate,
        translate=viewer.layers.selection.active.affine.translate[::-1],
    )
    affine = cut.affine_matrix
    # extract scale factor from affine object
    scale = (
        viewer.layers.selection.active.affine.scale[0],
        viewer.layers.selection.active.affine.scale[1],
        1,  # for 3D image
    )
    # calculate new shape from original and scale factor
    shape = np.array(
        np.array(fg_img.shape, dtype = np.float64) * scale,
        dtype = int,
    )  # keep original aspect ratio
    return affine, scale, shape
