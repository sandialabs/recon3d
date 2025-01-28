"""
This module contains functions related to the `downscale` command line argument 
with a provided input file.

Functions
---------
padded_size(img_stack_dim_size, target_res, original_res, tolerance, limit_factor)
    Determine the expected pad size along a dimension.

pad_amount(img_stack_dim_size, target_res, original_res, tolerance, limit_factor)
    Determine the amount of padding to add along a dimension.

apply_bbox(image_stack, threshold)
    Crop the image stack to the smallest bounding box over the threshold.

bbox_range(image_stack, threshold)
    Calculate the bounding box range in all dimensions.

save_downscale_stack(image_stack, path, folder_suffix)
    Save the new stack as a tiff image stack.

downscale(path_file_input)
    Downscale the image stack based on the provided input file.

main()
    Runs the module from the command line, invoked from pyproject.toml with 'downscale' command.
"""

import argparse
import itertools
import math
from pathlib import Path
from typing import Tuple

import numpy as np
from scipy import ndimage
from pyevtk.hl import gridToVTK

from recon3d.types import CartesianAxis3D
import recon3d.utility as ut


def apply_bbox(image_stack: np.ndarray, min_threshold: float) -> np.ndarray:
    """
    Crop the image stack to the smallest bounding box over the threshold.

    Parameters
    ----------
    image_stack : np.ndarray
        The stack of images as a data cube.
    min_threshold : float
        The minimum threshold value for cropping (values greater than this threshold will be retained).

    Returns
    -------
    np.ndarray
        The cropped image stack.

    Examples
    --------
    >>> image_stack = np.random.rand(10, 10, 10)
    >>> cropped_stack = apply_bbox(image_stack, 0.5)
    Cropping image stack to bounding box...
    """
    print("Cropping image stack to bounding box...")

    (z_start, z_end, y_start, y_end, x_start, x_end) = bbox_range(
        image_stack, min_threshold
    )

    # Slice the image stack to the bbox array
    # Add one to include the full range of the bounding box because of the slicing convention
    image_stack = image_stack[
        z_start : z_end + 1,
        y_start : y_end + 1,
        x_start : x_end + 1,
    ]
    return image_stack


def bbox_range(
    image_stack: np.ndarray, min_threshold: float
) -> Tuple[int, int, int, int, int, int]:
    """
    Calculate the bounding box range in all dimensions.

    Parameters
    ----------
    image_stack : np.ndarray
        The stack of images as a data cube.
    min_threshold : float
        The minimum threshold value for calculating the bounding box (values greater than this threshold will be retained).


    Returns
    -------
    tuple[int, int, int, int, int, int]
        The bounding box range in all dimensions.

    Examples
    --------
    >>> image_stack = np.ones((100, 100, 100))
    >>> bbox_range(image_stack, 0.5)
    (0, 99, 0, 99, 0, 99)
    """

    image_stack = image_stack > min_threshold
    # Determine the bbox range, see link
    # https://stackoverflow.com/questions/31400769/bounding-box-of-numpy-array

    # TODO: check this function with RGB image stack as input
    # N = image_stack.ndim
    N = 3
    bbox = []
    for ax in itertools.combinations(reversed(range(N)), N - 1):
        nonzero = np.any(image_stack, axis=ax)
        bbox.extend(np.where(nonzero)[0][[0, -1]])
    (z_start, z_end, y_start, y_end, x_start, x_end) = tuple(bbox)

    return (z_start, z_end, y_start, y_end, x_start, x_end)


def downscale(path_file_input: str) -> bool:
    """
    Downscale the image stack based on the provided input file.

    This function reads the input file, processes the image stack, and saves the downscaled stack.

    Parameters
    ----------
    path_file_input : str
        The path to the input file.

    Returns
    -------
    bool
        True if the processing is successful, False otherwise.
    """

    processed = False

    print(f"Processing file: {path_file_input}")

    db = ut.yaml_to_dict(Path(path_file_input))

    output_stack_type = db["output_stack_type"]
    if not (
        (output_stack_type == "downscaled")
        or (output_stack_type == "bounding_box")
        or (output_stack_type == "padded")
    ):
        raise ValueError(
            f'Invalid "output_stack_type" of "{output_stack_type}" input. Valid options inclue "downscaled", "bounding_box", or "padded".'
        )

    segmented_stack = ut.read_images(
        Path(db["image_dir"]).expanduser(), db["image_type"]
    )
    # zyx dimension ordering
    print(f"Original array size: {segmented_stack.shape}")

    dim_list = ["dz", "dy", "dx"]

    z_pad, y_pad, x_pad = [
        pad_amount(
            segmented_stack.shape[dim],
            db["resolution_output"][dim_list[dim]],
            db["resolution_input"][dim_list[dim]],
            db["downscale_tolerance"],
            db["image_limit_factor"],
        )
        for dim in range(0, 3)  # TODO should we restrict the dimensions to 3?
    ]

    padded_stack = np.pad(segmented_stack, (z_pad, y_pad, x_pad))
    print(f"New array size: {padded_stack.shape}")

    # Use list comprehension
    z_ratio, y_ratio, x_ratio = [
        db["resolution_input"][dim_list[dim]] / db["resolution_output"][dim_list[dim]]
        for dim in range(0, 3)  # TODO should we restrict the dimensions to 3
    ]

    downscaled_stack = ndimage.zoom(
        padded_stack,
        (z_ratio, y_ratio, x_ratio),
        order=0,
        mode="grid-constant",
        grid_mode=True,
    )

    print(f"Downscaled array size: {downscaled_stack.shape}")

    output_stack_type = db["output_stack_type"]
    if output_stack_type == "downscaled":
        output_stack = downscaled_stack

    else:
        cropped_stack = apply_bbox(downscaled_stack, 0)

        if output_stack_type == "bounding_box":
            output_stack = cropped_stack

        else:  # output_stack_type == "padded"
            # Pad after cropping
            (box_pad_z, box_pad_y, box_pad_x) = (
                db["padding"]["nz"],
                db["padding"]["ny"],
                db["padding"]["nx"],
            )
            padded_stack = np.pad(
                cropped_stack, ((box_pad_z,), (box_pad_y,), (box_pad_x,))
            )
            output_stack = padded_stack

    # output new downscaled stack
    out_dir = db["out_dir"]
    folder_suffix = f"{int(db['resolution_output']['dx'])}_dx"
    save_downscale_stack(output_stack, Path(out_dir), folder_suffix)

    # TODO: functionalize and test
    if db["save_npy"]:
        # npy_path = Path(out_dir).joinpath(f"{folder_suffix}.npy")
        npy_path = Path(out_dir).expanduser().joinpath(f"{folder_suffix}.npy")
        np.save(str(npy_path), output_stack)
        print(".npy file saved.")

    # TODO: functionalize and test
    if db["writeVTR"]:
        data = output_stack
        nx, ny, nz = data.shape[0], data.shape[1], data.shape[2]

        x = np.arange(0, nx + 1)
        y = np.arange(0, ny + 1)
        z = np.arange(0, nz + 1)

        # vtk_path = Path(out_dir).joinpath(f"{folder_suffix}")
        vtk_path = Path(out_dir).expanduser().joinpath(f"{folder_suffix}")
        gridToVTK(str(vtk_path), x, y, z, cellData={"imagedata": data})
        print(".vtr file saved.")

    processed = True  # overwrite, if we reach this point, all code is successful

    print(f"Finished processing file: {path_file_input}")

    return processed


def pad_amount(
    img_stack_dim_size: int,
    target_res: float,
    original_res: float,
    tolerance: float,
    limit_factor: float,
) -> tuple[int, int]:
    """
    Determine the amount of padding to add along a dimension.

    If the amount of padding is odd, the larger value will be padded to the front of the dimension.

    Parameters
    ----------
    img_stack_dim_size : int
        The size of the image stack along a dimension.
    target_res : float
        The target resolution.
    original_res : float
        The original resolution.
    tolerance : float
        The tolerance for the padding calculation.
    limit_factor : float
        The limit factor for the padding calculation.

    Returns
    -------
    tuple[int, int]
        The amount of padding to add to the front and back of the dimension.

    Examples
    --------
    >>> ds.pad_amount(100,0.62,1.0,0.01,2.0)
    (12, 12)
    """
    the_padded_size = padded_size(
        img_stack_dim_size,
        target_res,
        original_res,
        tolerance,
        limit_factor,
    )
    total_pad = the_padded_size - img_stack_dim_size
    pad = (math.ceil(total_pad / 2), math.floor(total_pad / 2))

    return pad


def padded_size(
    img_stack_dim_size: int,
    target_res: float,
    original_res: float,
    tolerance: float,
    limit_factor: float,
) -> int:
    """
    Determine the expected pad size along a dimension.

    Parameters
    ----------
    img_stack_dim_size : int
        The size of the image stack along a dimension.
    target_res : float
        The target resolution.
    original_res : float
        The original resolution.
    tolerance : float
        The tolerance for the padding calculation.
    limit_factor : float
        The limit factor for the padding calculation.

    Returns
    -------
    int
        The new dimension size after padding.

    Examples
    --------
    >>> padded_size(100, 0.5, 1.0, 0.01, 2.0)
    100
    """

    downscale_factor = float(target_res) / float(original_res)

    # start at the prior dimension size to enter the while loop
    new_dim = img_stack_dim_size - 1

    _too_high, _too_low = True, True

    _b = 0.0  # what we are comparing the result to

    while _too_high and _too_low:
        if new_dim > (img_stack_dim_size * int(limit_factor)):
            raise ValueError(
                f'Could not find even padding within "{tolerance}" for a {limit_factor}x image size. Increase tolerance or choose a different "target_res"'
            )
        new_dim += 1

        # to avoid the modulo (new_dim % downscale_factor) float error, so we must do manually
        _a_high = new_dim - (math.ceil(new_dim / downscale_factor) * downscale_factor)
        _a_low = new_dim - (math.floor(new_dim / downscale_factor) * downscale_factor)

        _too_high = not math.isclose(
            _a_high,
            _b,
            abs_tol=tolerance,
        )
        _too_low = not math.isclose(
            _a_low,
            _b,
            abs_tol=tolerance,
        )

    return new_dim


def save_downscale_stack(
    image_stack: np.ndarray, path: Path, folder_suffix: str
) -> bool:
    """
    Save the new stack as a tiff image stack.

    Parameters
    ----------
    image_stack : np.ndarray
        The stack of images as a data cube.
    path : Path
        The save path.
    folder_suffix : str
        The resolution in the suffix of the save folder "images_at_resolution_{RES_dx}".

    Returns
    -------
    bool
        True for success, False otherwise.

    Examples
    --------
    >>> image_stack = np.random.rand(100, 100, 100)
    >>> save_downscale_stack(image_stack, Path("/path/to/save"), "0.5_dx")
    Saving cropped image stack in /path/to/save > images_at_resolution_0.5_dx
    True
    """

    image_folder_name = f"images_at_resolution_{folder_suffix}"
    print(f"Saving cropped image stack in {path} > {image_folder_name}")
    ut.ndarray_to_img(
        data=image_stack,
        slice_axis=CartesianAxis3D.Z,
        parent_dir=path,
        folder_name=image_folder_name,
    )

    return True


def main():
    """
    Runs the module from the command line, invoked from pyproject.toml with 'downscale' command.

    This function sets up the command line argument parser, parses the input arguments,
    and calls the `downscale` function with the provided input file.

    Parameters
    ----------
    None

    Returns
    -------
    None

    Examples
    --------
    To run this module from the command line with the 'downscale' command:

        $ downscale input_file.yml
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="the .yml user input file")
    args = parser.parse_args()
    input_file = args.input_file

    downscale(path_file_input=input_file)


if __name__ == "__main__":
    main()
