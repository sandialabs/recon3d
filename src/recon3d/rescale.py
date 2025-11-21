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

# TODO fix docstring

import argparse
import itertools
import math
from pathlib import Path
from typing import Tuple, Union, Dict, List, Any

import numpy as np
from scipy import ndimage
from pyevtk.hl import gridToVTK
from skimage.transform import resize as sk_resize

import recon3d.types as rtt
import recon3d.utility as ut


def parse_config(d: dict) -> rtt.RescaleConfig:
    """
    Turn the raw YAML dict into a Config NamedTuple.
    Both `padding:` and `final_size:` are optional in the YAML.
    """
    # padding is only used by PAD_AFTER_CROP
    raw_pad = d.get("padding", {})  # e.g. {"nz":4, "ny":[5,4], ...}
    padding: Dict[str, Tuple[int, int]] = {}
    for axis in ("nz", "ny", "nx"):
        if axis not in raw_pad:
            continue
        v = raw_pad[axis]
        if isinstance(v, (list, tuple)):
            if len(v) != 2:
                raise ValueError(f"padding.{axis} must be an int or 2‐tuple, got {v}")
            padding[axis] = (int(v[0]), int(v[1]))
        else:
            # single int --> uniform padding both sides
            n = int(v)
            padding[axis] = (n, n)

    # final_size is only used by PAD_TO_SIZE
    raw_fs = d.get("final_size", None)
    final_size = None
    if raw_fs is not None:
        final_size = rtt.FinalSize(
            nz=int(raw_fs["nz"]),
            ny=int(raw_fs["ny"]),
            nx=int(raw_fs["nx"]),
        )

    raw_pp = d.get("post_process", [])
    if isinstance(raw_pp, str):
        post_process = [raw_pp]
    else:
        post_process = list(raw_pp)

    # 2) histogram stretch params (only used if in post_process)
    clip = d.get("clip_percentiles", [1.0, 99.0])
    clip_percentiles = (float(clip[0]), float(clip[1]))

    raw_ignore = d.get("ignore_values", [0])
    if isinstance(raw_ignore, (list, tuple)):
        ignore_values = tuple(int(x) for x in raw_ignore)
    else:
        ignore_values = (int(raw_ignore),)

    return rtt.RescaleConfig(
        image_dir=Path(d["image_dir"]).expanduser(),
        image_type=d["image_type"],
        out_dir=Path(d["out_dir"]).expanduser(),
        resolution_input=d["resolution_input"],
        resolution_output=d["resolution_output"],
        rescale_tolerance=d["rescale_tolerance"],
        image_limit_factor=d["image_limit_factor"],
        interpolation_mode=rtt.InterpolationMode.from_string(d["interpolation_mode"]),
        output_stack_type=rtt.OutputStackType.from_string(d["output_stack_type"]),
        padding=padding,
        final_size=final_size,
        save_npy=d["save_npy"],
        writeVTR=d["writeVTR"],
        bbox_threshold=d.get("bbox_threshold", 0.0),
        post_process=post_process,
        clip_percentiles=clip_percentiles,
        ignore_values=ignore_values,
    )


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


def save_rescale_stack(image_stack: np.ndarray, path: Path, folder_suffix: str) -> bool:
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
        slice_axis=rtt.CartesianAxis3D.Z,
        parent_dir=path,
        folder_name=image_folder_name,
    )

    return True


def pad_to_final_size(arr: np.ndarray, final_size: rtt.FinalSize) -> np.ndarray:
    """
    Pad (only spatial dims) so that arr.shape[:3] --> final_size exactly,
    splitting extra voxels front/back as evenly as possible.
    """
    # arr.shape[:3] --> (z,y,x)
    z, y, x = arr.shape[:3]
    # unpack
    fz, fy, fx = final_size.nz, final_size.ny, final_size.nx

    def split(delta: int) -> Tuple[int, int]:
        return (math.ceil(delta / 2), math.floor(delta / 2))

    pad_z = split(fz - z)
    pad_y = split(fy - y)
    pad_x = split(fx - x)

    padded = np.pad(
        arr, (pad_z, pad_y, pad_x, (0, 0)), mode="constant", constant_values=0
    )
    print(f"Padded to final_size {final_size} --> {padded.shape}")
    return padded


def write_vtr(arr: np.ndarray, out_dir: Path, suffix: str):
    """Permute to x,y,z,c order and write via gridToVTK."""
    data = arr.transpose(2, 1, 0, 3)
    nx, ny, nz, nc = data.shape
    x = np.arange(nx + 1, dtype=np.float32)
    y = np.arange(ny + 1, dtype=np.float32)
    z = np.arange(nz + 1, dtype=np.float32)

    cell_fields = {}
    if nc == 1:
        cell_fields["imagedata"] = data[..., 0]
    elif nc == 3:
        cell_fields["imagedata"] = (data[..., 0], data[..., 1], data[..., 2])
    else:
        for c in range(nc):
            cell_fields[f"imagedata_c{c}"] = data[..., c]

    vtk_path = out_dir.joinpath(suffix)
    gridToVTK(str(vtk_path), x, y, z, cellData=cell_fields)
    print(f".vtr written to {vtk_path!s}")


def rescale_stack(stack: np.ndarray, cfg: rtt.RescaleConfig) -> np.ndarray:
    """
    1) Compute & apply input‐padding so that zoom yields integer dims
    2) Zoom
    3) Depending on cfg.output_stack_type:
        - crop to bbox
        - pad after crop (cfg.padding)
        - pad to cfg.final_size
    """
    # 1) input‐pad
    z_in, y_in, x_in = stack.shape[:3]
    z_pad = pad_amount(
        z_in,
        cfg.resolution_output["dz"],
        cfg.resolution_input["dz"],
        cfg.rescale_tolerance,
        cfg.image_limit_factor,
    )
    y_pad = pad_amount(
        y_in,
        cfg.resolution_output["dy"],
        cfg.resolution_input["dy"],
        cfg.rescale_tolerance,
        cfg.image_limit_factor,
    )
    x_pad = pad_amount(
        x_in,
        cfg.resolution_output["dx"],
        cfg.resolution_input["dx"],
        cfg.rescale_tolerance,
        cfg.image_limit_factor,
    )
    c_pad = (0, 0)  # no channel padding
    padded = np.pad(
        stack, (z_pad, y_pad, x_pad, c_pad), mode="constant", constant_values=0
    )
    print(f"Padded input to {padded.shape}")

    # 2) zoom
    zf = cfg.resolution_input["dz"] / cfg.resolution_output["dz"]
    yf = cfg.resolution_input["dy"] / cfg.resolution_output["dy"]
    xf = cfg.resolution_input["dx"] / cfg.resolution_output["dx"]
    factors = (zf, yf, xf, 1.0)
    # zoomed = ndimage.zoom(
    #     padded,
    #     factors,
    #     order=cfg.interpolation_mode.value,
    #     mode="grid-constant",
    #     grid_mode=True,
    # )
    zoomed = sk_resize(
        padded,
        output_shape=(
            int(padded.shape[0] * factors[0]),
            int(padded.shape[1] * factors[1]),
            int(padded.shape[2] * factors[2]),
            padded.shape[3],
        ),
        order=cfg.interpolation_mode.value,  # 0=nearest,1=bilinear,3=bicubic, etc.
        mode="constant",
        cval=0,
        clip=True,
        preserve_range=True,
        anti_aliasing=False,
    ).astype(padded.dtype)
    print(f"Zoomed to {zoomed.shape}")

    # 3) post‐processing
    mode = cfg.output_stack_type
    if mode is rtt.OutputStackType.RESCALED:
        return zoomed

    # crop to bounding box first
    cropped = apply_bbox(zoomed, cfg.bbox_threshold)
    print(f"Cropped--> {cropped.shape}")

    if mode is rtt.OutputStackType.BOUNDING_BOX:
        return cropped

    if mode is rtt.OutputStackType.PAD_AFTER_CROP:
        # padding specified in cfg.padding as {"nz":(before,after), ...}
        pad_spec = (
            cfg.padding["nz"],
            cfg.padding["ny"],
            cfg.padding["nx"],
            (0, 0),
        )
        padded2 = np.pad(cropped, pad_spec, mode="constant", constant_values=0)
        print(f"Padded after crop --> {padded2.shape}")
        return padded2

    # PAD_TO_SIZE
    if mode is rtt.OutputStackType.PAD_TO_SIZE:
        if cfg.final_size is None:
            raise ValueError("final_size must be set for PAD_TO_SIZE mode")
        return pad_to_final_size(cropped, cfg.final_size)

    # unreachable
    raise RuntimeError(f"Unhandled output_stack_type: {mode}")


def plan_uniform_downscale(
    dims_list: List[Tuple[int, int, int, int]],
    res_list: List[Tuple[float, float, float]],
    target_res: Tuple[float, float, float],
    tolerance: float,
    limit_factor: float,
) -> Tuple[Tuple[int, int, int], List[Dict[str, Any]]]:
    """
    Given a list of input volumes (dims,res) and a desired target_res,
    figure out for each volume
      1) how much to pad the input so the zoom factors are integral
      2) what the exact zoomed shape will be
    then take the elementwise max over all zoomed shapes to get a
    single uniform final (Z,Y,X).

    Returns:
      final_shape:  (Fz, Fy, Fx)
      plans:        list of per‐volume dicts:
        {
          "orig_dims":   (nz,ny,nx,c),
          "in_res":      (dz,dy,dx),
          "pad_input":   {"z":(fz,bz), "y":(fy,by), "x":(fx,bx)},
          "zoomed":      (Z, Y, X),
        }
    """
    plans: List[Dict[str, Any]] = []

    for (nz, ny, nx, _c), (dz, dy, dx) in zip(dims_list, res_list):
        # compute how much to pad the *input* so that padded_n / downscale_factor
        # is within tolerance of an integer
        pad_z = pad_amount(nz, target_res[0], dz, tolerance, limit_factor)
        pad_y = pad_amount(ny, target_res[1], dy, tolerance, limit_factor)
        pad_x = pad_amount(nx, target_res[2], dx, tolerance, limit_factor)

        padded_nz = nz + pad_z[0] + pad_z[1]
        padded_ny = ny + pad_y[0] + pad_y[1]
        padded_nx = nx + pad_x[0] + pad_x[1]

        # compute the zoomed shape:
        # zoom_factor = original_res / target_res
        Z = int(round(padded_nz * dz / target_res[0]))
        Y = int(round(padded_ny * dy / target_res[1]))
        X = int(round(padded_nx * dx / target_res[2]))

        plans.append(
            {
                "orig_dims": (nz, ny, nx, _c),
                "in_res": (dz, dy, dx),
                "pad_input": {"z": pad_z, "y": pad_y, "x": pad_x},
                "zoomed": (Z, Y, X),
            }
        )

    # elementwise maximum of all zoomed shapes --> uniform final shape
    final_Z = max(p["zoomed"][0] for p in plans)
    final_Y = max(p["zoomed"][1] for p in plans)
    final_X = max(p["zoomed"][2] for p in plans)

    return (final_Z, final_Y, final_X), plans


def filter_by_resolution_range(
    metas: List[Dict[str, Any]], axis: str, min_val: float, max_val: float
) -> List[Dict[str, Any]]:
    """
    Return only those metadata dicts whose resolution[axis] is in [min_val, max_val].

    metas is the list of dicts you built in standardize_geometry, each like:
      { "folder": Path(...),
        "resolution": (dz,dy,dx),
        "dimensions": (nz,ny,nx,c), }
    axis must be one of "dz","dy","dx".
    """
    idx = {"dz": 0, "dy": 1, "dx": 2}[axis]
    out = []
    for m in metas:
        r = m["resolution"][idx]
        if min_val <= r <= max_val:
            out.append(m)
    return out


def rescale_from_yaml(yaml_path: Union[str, Path]) -> bool:
    """
    Read the YAML, load images, rescale, save, optionally write VTR/npy.
    """
    yaml_path = Path(yaml_path)
    cfg = parse_config(ut.yaml_to_dict(yaml_path))
    stack = ut.read_images(cfg.image_dir, cfg.image_type)
    print(f"Original array size: {stack.shape}")

    out_stack = rescale_stack(stack, cfg)

    # postprocessing
    # now apply any post‐processing steps in order
    for step in cfg.post_process:
        if step.lower() == "histogram_stretch":
            # import your updated function
            from recon3d.process import histogram_stretch

            out_stack = histogram_stretch(
                out_stack,
                clip_percentiles=cfg.clip_percentiles,
                ignore_values=cfg.ignore_values,
            )
            print(
                f"Applied histogram_stretch; clip={cfg.clip_percentiles}, ignore={cfg.ignore_values}"
            )

        else:
            raise ValueError(f"Unknown post_process step '{step}'")

    # save to TIFFs
    suffix = f"{int(cfg.resolution_output['dx'])}_dx"
    save_rescale_stack(out_stack, cfg.out_dir, suffix)

    if cfg.save_npy:
        npy_path = cfg.out_dir.joinpath(f"{suffix}.npy")
        np.save(str(npy_path), out_stack)
        print(f".npy saved to {npy_path!s}")

    if cfg.writeVTR:
        write_vtr(out_stack, cfg.out_dir, suffix)

    print(f"Finished processing {yaml_path!s}")
    return True


def main():
    """
    Runs the module from the command line, invoked from pyproject.toml with 'rescale' command.

    This function sets up the command line argument parser, parses the input arguments,
    and calls the `rescale` function with the provided input file.

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

    rescale_from_yaml(yaml_path=input_file)


if __name__ == "__main__":
    main()
