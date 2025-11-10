"""
This module holds utilties that can be reused across other modules
"""

# Standard library imports
from datetime import datetime
from itertools import cycle, repeat, tee
from pathlib import Path
from typing import Iterable, Tuple
import glob
import argparse
import shutil
from typing import Optional

# Third-party library imports
import h5py
import numpy as np
import yaml
from PIL import Image
from scipy import ndimage

# Local imports
import recon3d.types as rtt


def binary_with_pores_to_semantic(input_path: Path, output_path: Path) -> dict:
    """
    Convert a folder of segmented/binarized TIFF images to a Semantic Image Stack.

    This function is designed for a specific use case (NOMAD, pores in AM tensile bars)
    and converts binary images with metal as 1 and all else as 0 into a semantic image stack.

    Parameters
    ----------
    input_path : Path
        Folder containing binary images with metal as 1 and all else as 0.
        Figure 1:
        *---*---*---*---*
        | 1 | 1 | 1 | 1 |
        *---*---*---*---*
        | 1 | 0 | 1 | 0 |
        *---*---*---*---*
        | 1 | 1 | 1 | 0 |
        *---*---*---*---*
    output_path : Path
        Directory where the semantic images will be saved.
        Figure 2:
        *---*---*---*---*
        | 1 | 1 | 1 | 1 |
        *---*---*---*---*
        | 1 | 2 | 1 | 0 |
        *---*---*---*---*
        | 1 | 1 | 1 | 0 |
        *---*---*---*---*

    Returns
    -------
    dict
        Dictionary of class labels:
            class_labels:
                air:
                    value: 0
                metal:
                    value: 1
                pore:
                    value: 2

    Raises
    ------
    ValueError
        If the input images contain more than two phases.

    Examples
    --------
    >>> input_path = Path("path/to/binary_images")
    >>> output_path = Path("path/to/save/semantic_images")
    >>> binary_with_pores_to_semantic(input_path, output_path)
    {'class_labels': {'air': {'value': 0}, 'metal': {'value': 1}, 'pore': {'value': 2}}}

    """

    # read in binary image stack
    bw_data = read_images(input_path)

    # check data is binary, only 0 or 1
    phase_ids = np.unique(bw_data)
    n_phases = len(phase_ids)

    if n_phases != 2:
        raise ValueError(
            f"Only two phases expected in segmented image, {n_phases} phase found"
        )

    # preallocate the output data
    output_data = np.zeros(shape=bw_data.shape, dtype=np.uint8)

    # fill holes and assign as 'metal', 1
    print("\tIsolating Sample...")
    sample = np.zeros_like(bw_data, dtype=bool)
    for c in range(bw_data.shape[-1]):
        sample[..., c] = ndimage.binary_fill_holes(bw_data[..., c])
    np.place(output_data, sample, 1)

    # isolate holes within 'metal', assign as 'pore', 2
    print("\tIsolating Voids...")
    voids = np.zeros_like(bw_data, dtype=bool)
    for c in range(bw_data.shape[-1]):
        voids[..., c] = np.logical_xor(sample[..., c], bw_data[..., c])
    np.place(output_data, voids, 2)

    # thus, everything else is 'air', 0

    class_labels = {
        "class_labels": {
            "air": {"value": 0},
            "metal": {"value": 1},
            "pore": {"value": 2},
        }
    }

    ndarray_to_img(
        data=output_data,
        slice_axis=rtt.CartesianAxis3D.Z,
        parent_dir=output_path,
        folder_name="",
    )

    return class_labels


def validate_yml(yml_input_file: Path, cli_entry_point: str) -> tuple[Path, Path]:
    """
    Verify that the YAML file contains the required arguments.

    This function checks if the specified YAML file contains the necessary arguments
    and validates the CLI entry point. It also verifies the existence of input and output directories.

    Parameters
    ----------
    yml_input_file : Path
        The path to the YAML input file.
    cli_entry_point : str
        The CLI entry point to validate against the YAML file.

    Returns
    -------
    tuple[Path, Path, dict]
        A tuple containing:
        - Path to the input directory.
        - Path to the output directory.
        - Dictionary of YAML values.

    Raises
    ------
    ValueError
        If the CLI entry point is not found in the YAML file or if the input directory does not exist.

    Examples
    --------
    >>> yml_input_file = Path("path/to/input.yml")
    >>> cli_entry_point = "process_images"
    >>> validate_yml(yml_input_file, cli_entry_point)
    (PosixPath('path/to/images'), PosixPath('path/to/output'), {'cli_entry_points': ['process_images'], 'image_dir': 'path/to/images', 'out_dir': 'path/to/output'})
    """

    print(f"Processing specification file: {yml_input_file}")
    yml_vals = yaml_to_dict(yml_input_file)

    # check cli_entry_points is valid
    if cli_entry_point not in yml_vals["cli_entry_points"]:
        raise ValueError(
            f"""Error. Incorrect yml format.
                         This function requires the "cli_entry_points" key to contain {cli_entry_point},
                         but currently contains the following options: {yml_vals["cli_entry_points"]} """
        )

    path_input = Path(yml_vals["image_dir"]).expanduser()
    if not path_input.is_dir():
        raise ValueError(f"Error, 'image_dir', {path_input} not found.")
    print(f"Input path: {path_input}")

    path_output = Path(yml_vals["out_dir"]).expanduser()
    path_output.mkdir(parents=True, exist_ok=True)
    print(f"Output path: {path_output}")

    return (path_input, path_output, yml_vals)


def binarize(data: np.ndarray, val: int) -> np.ndarray:
    """
    Binarize the data based on a specified value.

    This function converts the values within the data matching the specified value to 1,
    and all other values to 0.

    Parameters
    ----------
    data : np.ndarray
        The input array to be binarized.
    val : int
        The value in the data to be binarized (converted to 1).

    Returns
    -------
    np.ndarray
        The binarized array with the same shape as the input data.

    Examples
    --------
    >>> data = np.array([[1, 2, 3], [4, 1, 6], [7, 8, 1]])
    >>> binarize(data, 1)
    array([[1, 0, 0],
           [0, 1, 0],
           [0, 0, 1]], dtype=int8)
    """

    bw_data = np.zeros(shape=data.shape, dtype=np.bool_)
    np.place(bw_data, data == val, 1)

    return bw_data


def semantic_to_binary(yml_input_file: Path) -> bool:
    """
    Convert a semantic image stack to a binary image stack using a YAML configuration file.

    This function reads a YAML file to prepare a binary image stack from a semantic image stack.
    The selected class from the semantic stack is binarized based on the specified value.

    Parameters
    ----------
    yml_input_file : Path
        The path to the YAML input file containing configuration settings.

    Returns
    -------
    bool
        True if the binary image stack was successfully created, False otherwise.

    Examples
    --------
    >>> yml_input_file = Path("path/to/config.yml")
    >>> semantic_to_binary(yml_input_file)
    True
    """

    [input_path, output_path, params] = validate_yml(
        yml_input_file, "semantic_to_binary"
    )

    # params = ut.yaml_to_dict(yml_input_file)

    sel_class = params["selected_class"]
    sel_class_value = params["class_labels"][sel_class]["value"]

    input_data = read_images(input_path)

    bw_data = binarize(input_data, sel_class_value)

    ndarray_to_img(
        data=bw_data,
        slice_axis=rtt.CartesianAxis3D.Z,
        parent_dir=output_path,
        folder_name="",
    )


def binary_to_semantic(yml_input_file: Path) -> bool:
    """
    Convert a binary image stack to a semantic image stack using a YAML configuration file.

    This function reads a YAML file to prepare a semantic image stack from a binary image stack.
    It is designed for a specific use case (NOMAD, AM tensile bars) as of 31 May 2024.

    Parameters
    ----------
    yml_input_file : Path
        The path to the YAML input file containing configuration settings.

    Returns
    -------
    bool
        True if the semantic image stack was successfully created, False otherwise.

    Examples
    --------
    >>> yml_input_file = Path("path/to/config.yml")
    >>> binary_to_semantic(yml_input_file)
    True
    """

    [input_path, output_path, params] = validate_yml(
        yml_input_file, "binary_to_semantic"
    )

    class_labels = binary_with_pores_to_semantic(
        input_path=input_path, output_path=output_path
    )

    print(f"class labels for semantic stack:\n{class_labels}")


def main_binary_to_semantic():
    """
    Run the binary to semantic conversion module from the command line.

    This function serves as the entry point for terminal-based access to the binary to semantic conversion module.
    It is invoked from the command line using the 'binary_to_semantic' command specified in pyproject.toml.
    The function processes a YAML input file to convert a binary image stack to a semantic image stack.

    Parameters
    ----------
    None

    Returns
    -------
    None

    Examples
    --------
    To run the binary to semantic conversion, use the following command in the terminal:
    $ binary_to_semantic path/to/input.yml
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="the .yml user input file")
    args = parser.parse_args()
    yml_input_file = args.input_file

    binary_to_semantic(yml_input_file=yml_input_file)

    print(f"{yml_input_file} processed!")


def main_semantic_to_binary():
    """
    Run the semantic to binary conversion module from the command line.

    This function serves as the entry point for terminal-based access to the semantic to binary conversion module.
    It is invoked from the command line using the 'semantic_to_binary' command specified in pyproject.toml.
    The function processes a YAML input file to convert a semantic image stack to a binary image stack.

    Parameters
    ----------
    None

    Returns
    -------
    None

    Examples
    --------
    To run the semantic to binary conversion, use the following command in the terminal:
    $ semantic_to_binary path/to/input.yml
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="the .yml user input file")
    args = parser.parse_args()
    yml_input_file = args.input_file

    semantic_to_binary(yml_input_file=yml_input_file)

    print(f"{yml_input_file} processed!")


def hdf_to_instance_properties(
    hdf_path: Path, group_path: str
) -> rtt.InstanceProperties:
    """
    Read instance analysis data from an HDF5 file and create an InstanceProperties object.

    This function reads data from an HDF5 file that has had instance analysis performed
    with a pre-defined internal structure to create an InstanceProperties object.
    The group name should be the internal 'folder' in the HDF5 file containing the instance
    analysis property data.

    Parameters
    ----------
    hdf_path : Path
        The path to the HDF5 file containing the instance analysis data.
    group_path : str
        The internal 'folder' in the HDF5 file containing the instance analysis property data.

    Returns
    -------
    InstanceProperties
        An InstanceProperties object containing the instance analysis data.

    Examples
    --------
    >>> hdf_path = Path("path/to/instance_analysis.h5")
    >>> group_path = "instance_properties"
    >>> instance_properties = hdf_to_instance_properties(hdf_path, group_path)
    >>> print(instance_properties)
    InstanceProperties(source_name='instance_properties', labels=InstanceLabels(data=array([0, 1, 2, ...])), n_voxels=[NVoxel(value=...), ...], equivalent_sphere_diameters=[Length(value=..., unit=<Units.MICRON: 'micron'>), ...], centroids=Centroids(data=[Centroid(cx=Length(value=..., unit=<Units.MICRON: 'micron'>), cy=Length(value=..., unit=<Units.MICRON: 'micron'>), cz=Length(value=..., unit=<Units.MICRON: 'micron'>)), ...]), ellipsoids=BestFitEllipsoids(data=[BestFitEllipsoid(a=EllipsoidAxis(length=Length(value=..., unit=<Units.MICRON: 'micron'>), orientation=UnitVector(u=..., v=..., w=...)), b=EllipsoidAxis(length=Length(value=..., unit=<Units.MICRON: 'micron'>), orientation=UnitVector(u=..., v=..., w=...)), c=EllipsoidAxis(length=Length(value=..., unit=<Units.MICRON: 'micron'>), orientation=UnitVector(u=..., v=..., w=...))), ...]), surface_areas=EllipsoidSurfaceAreas(data=[Area(value=..., unit_squared=<Units.MICRON: 'micron'>), ...]), volumes=EllipsoidVolumes(data=[Volume(value=..., unit_cubed=<Units.MICRON: 'micron'>), ...]))
    """

    with h5py.File(hdf_path, "r") as f:
        instance_data = f[group_path]
        hdf_n_voxels = np.squeeze(instance_data["num_voxels"][:])
        hdf_equiv_diam = np.squeeze(instance_data["equivalent_sphere_diameters"][:])
        equiv_diam_unit = rtt.Units(
            instance_data["equivalent_sphere_diameters"].attrs["units"]
        )
        hdf_centroids = np.squeeze(instance_data["centroids"][:])
        centroid_unit = rtt.Units(instance_data["centroids"].attrs["units"])
        hdf_semi_axes = np.squeeze(instance_data["semi-axis_lengths"][:])
        semi_axes_unit = rtt.Units(instance_data["semi-axis_lengths"].attrs["units"])
        hdf_vectors = np.squeeze(instance_data["axis_vectors"][:])
        hdf_surface_areas = np.squeeze(instance_data["ellipsoid_surface_areas"][:])
        surface_area_unit = rtt.Units(
            instance_data["ellipsoid_surface_areas"].attrs["units_squared"]
        )
        hdf_volumes = np.squeeze(instance_data["ellipsoid_volumes"][:])
        volume_unit = rtt.Units(instance_data["ellipsoid_volumes"].attrs["units_cubed"])

    labels = np.arange(0, len(hdf_equiv_diam), dtype=int)
    n_voxels = [rtt.NVoxel(value=i) for i in hdf_n_voxels]
    equivalent_sphere_diameters = [
        rtt.Length(value=i, unit=equiv_diam_unit) for i in hdf_equiv_diam
    ]
    centroids = rtt.Centroids(
        data=[
            rtt.Centroid(
                cx=rtt.Length(value=i[0], unit=centroid_unit),
                cy=rtt.Length(value=i[1], unit=centroid_unit),
                cz=rtt.Length(value=i[2], unit=centroid_unit),
            )
            for i in hdf_centroids
        ]
    )

    ellipsoids = rtt.BestFitEllipsoids(
        [
            rtt.BestFitEllipsoid(
                a=rtt.EllipsoidAxis(
                    length=rtt.Length(value=i[0], unit=semi_axes_unit),
                    orientation=rtt.UnitVector(u=j[0], v=j[1], w=j[2]),
                ),
                b=rtt.EllipsoidAxis(
                    length=rtt.Length(value=i[1], unit=semi_axes_unit),
                    orientation=rtt.UnitVector(u=j[3], v=j[4], w=j[5]),
                ),
                c=rtt.EllipsoidAxis(
                    length=rtt.Length(value=i[2], unit=semi_axes_unit),
                    orientation=rtt.UnitVector(u=j[6], v=j[7], w=j[8]),
                ),
            )
            for i, j in zip(hdf_semi_axes, hdf_vectors)
        ]
    )

    surface_areas = rtt.EllipsoidSurfaceAreas(
        [rtt.Area(value=i, unit_squared=surface_area_unit) for i in hdf_surface_areas]
    )

    volumes = rtt.EllipsoidVolumes(
        [rtt.Volume(value=i, unit_cubed=volume_unit) for i in hdf_volumes]
    )

    instance_props = rtt.InstanceProperties(
        source_name=group_path,
        labels=rtt.InstanceLabels(data=labels),
        n_voxels=n_voxels,
        equivalent_sphere_diameters=equivalent_sphere_diameters,
        centroids=centroids,
        ellipsoids=ellipsoids,
        surface_areas=surface_areas,
        volumes=volumes,
    )

    # n_voxels: list[NVoxel]  # could have a InstanceImageStack
    # equivalent_sphere_diameters: list[Length]
    # centroids: Centroids  # could have a InstanceImageStack
    # ellipsoids: BestFitEllipsoids
    # surface_areas: SurfaceAreas)

    return instance_props


def hdf_to_metadata(hdf_path: Path, dataset_path: str) -> rtt.MetaData:
    """
    Extract metadata from an HDF5 dataset containing string metadata attributes.

    This function reads metadata attributes from a specified dataset within an HDF5 file
    and returns a MetaData object.

    Parameters
    ----------
    hdf_path : Path
        The path to the HDF5 file containing the dataset.
    dataset_path : str
        The internal path to the dataset within the HDF5 file.

    Returns
    -------
    MetaData
        A MetaData object containing the extracted metadata.

    Examples
    --------
    >>> hdf_path = Path("path/to/data.h5")
    >>> dataset_path = "dataset"
    >>> metadata = hdf_to_metadata(hdf_path, dataset_path)
    >>> print(metadata)
    MetaData(data_volume=DataVolume(x_width=..., y_height=..., z_image_count=...), resolution=Resolution(dx=Length(value=..., unit=<Units.MICRON: 'micron'>), dy=Length(value=..., unit=<Units.MICRON: 'micron'>), dz=Length(value=..., unit=<Units.MICRON: 'micron'>)), pixel_units=<Units.MICRON: 'micron'>, origin=Origin(x0=Length(value=..., unit=<Units.MICRON: 'micron'>), y0=Length(value=..., unit=<Units.MICRON: 'micron'>), z0=Length(value=..., unit=<Units.MICRON: 'micron'>)))
    """

    with h5py.File(hdf_path, "r") as f:
        dataset = f[dataset_path]

        x_width = int(dataset.attrs["x_width (pixels)"])
        y_height = int(dataset.attrs["y_height (pixels)"])
        z_image_count = int(dataset.attrs["z_image_count (pixels)"])
        data_volume = rtt.DataVolume(
            x_width=x_width, y_height=y_height, z_image_count=z_image_count
        )

        pixel_units = rtt.Units(dataset.attrs["Resolution, units"])
        dx = rtt.Length(value=float(dataset.attrs["Resolution, dx"]), unit=pixel_units)
        dy = rtt.Length(value=float(dataset.attrs["Resolution, dy"]), unit=pixel_units)
        dz = rtt.Length(value=float(dataset.attrs["Resolution, dz"]), unit=pixel_units)
        resolution = rtt.Resolution(
            dx=dx,
            dy=dy,
            dz=dz,
        )

        origin_units = rtt.Units(dataset.attrs["Origin, units"])
        x0 = rtt.Length(value=float(dataset.attrs["Origin, x0"]), unit=origin_units)
        y0 = rtt.Length(value=float(dataset.attrs["Origin, y0"]), unit=origin_units)
        z0 = rtt.Length(value=float(dataset.attrs["Origin, z0"]), unit=origin_units)
        origin = rtt.Origin(
            x0=x0,
            y0=y0,
            z0=z0,
        )

    metadata = rtt.MetaData(
        data_volume=data_volume,
        resolution=resolution,
        pixel_units=pixel_units,
        origin=origin,
    )

    return metadata


def centroids_to_ndarray(centroids: rtt.Centroids) -> np.ndarray:
    """
    Convert centroid data type into a NumPy array.

    This function converts a Centroids object into a NumPy array with each row representing
    the (cx, cy, cz) coordinates of a centroid.

    Parameters
    ----------
    centroids : Centroids
        The Centroids object containing the centroid data.

    Returns
    -------
    np.ndarray
        A NumPy array with shape (n, 3), where n is the number of centroids, and each row
        contains the (cx, cy, cz) coordinates of a centroid.

    Examples
    --------
    >>> centroids = Centroids(data=[
    ...     Centroid(cx=Length(0.0, Units.MICRON), cy=Length(0.0, Units.MICRON), cz=Length(0.0, Units.MICRON)),
    ...     Centroid(cx=Length(1.0, Units.MICRON), cy=Length(1.0, Units.MICRON), cz=Length(1.0, Units.MICRON)),
    ...     Centroid(cx=Length(2.0, Units.MICRON), cy=Length(2.0, Units.MICRON), cz=Length(2.0, Units.MICRON))
    ... ])
    >>> centroids_to_ndarray(centroids)
    array([[0., 0., 0.],
           [1., 1., 1.],
           [2., 2., 2.]])
    """

    data = centroids.data
    cx = [i.cx.value for i in data]
    cy = [i.cy.value for i in data]
    cz = [i.cz.value for i in data]

    ndarray = np.array((cx, cy, cz), dtype=float).T

    return ndarray


def ellipsoids_to_ndarray(
    ellipsoids: rtt.BestFitEllipsoids,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convert ellipsoid data type into NumPy arrays.

    This function converts a BestFitEllipsoids object into two NumPy arrays:
    one for the axis lengths and one for the axis vectors.

    Parameters
    ----------
    ellipsoids : BestFitEllipsoids
        The BestFitEllipsoids object containing the ellipsoid data.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        A tuple containing:
        - A NumPy array with shape (n, 3) for the axis lengths, where n is the number of ellipsoids.
        - A NumPy array with shape (n, 9) for the axis vectors, where n is the number of ellipsoids.

    Examples
    --------
    >>> ellipsoids = BestFitEllipsoids(data=[
    ...     BestFitEllipsoid(
    ...         a=EllipsoidAxis(length=Length(5.0, Units.MICRON), orientation=UnitVector(u=1.0, v=0.0, w=0.0)),
    ...         b=EllipsoidAxis(length=Length(3.0, Units.MICRON), orientation=UnitVector(u=0.0, v=1.0, w=0.0)),
    ...         c=EllipsoidAxis(length=Length(2.0, Units.MICRON), orientation=UnitVector(u=0.0, v=0.0, w=1.0))
    ...     ),
    ...     BestFitEllipsoid(
    ...         a=EllipsoidAxis(length=Length(6.0, Units.MICRON), orientation=UnitVector(u=0.707, v=0.707, w=0.0)),
    ...         b=EllipsoidAxis(length=Length(4.0, Units.MICRON), orientation=UnitVector(u=-0.707, v=0.707, w=0.0)),
    ...         c=EllipsoidAxis(length=Length(3.0, Units.MICRON), orientation=UnitVector(u=0.0, v=0.0, w=1.0))
    ...     )
    ... ])
    >>> axis_lengths, axis_vectors = ellipsoids_to_ndarray(ellipsoids)
    >>> print(axis_lengths)
    array([[5., 3., 2.],
           [6., 4., 3.]])
    >>> print(axis_vectors)
    array([[ 1.   ,  0.   ,  0.   ,  0.   ,  1.   ,  0.   ,  0.   ,  0.   ,  1.   ],
           [ 0.707,  0.707,  0.   , -0.707,  0.707,  0.   ,  0.   ,  0.   ,  1.   ]])
    """

    data = ellipsoids.data
    # axis lengths:
    a = [i.a.length.value for i in data]
    b = [i.b.length.value for i in data]
    c = [i.c.length.value for i in data]
    axis_lengths = np.array((a, b, c), dtype=float).T

    # axes vectors
    a_u = [i.a.orientation.u for i in data]
    a_v = [i.a.orientation.v for i in data]
    a_w = [i.a.orientation.w for i in data]

    b_u = [i.b.orientation.u for i in data]
    b_v = [i.b.orientation.v for i in data]
    b_w = [i.b.orientation.w for i in data]

    c_u = [i.c.orientation.u for i in data]
    c_v = [i.c.orientation.v for i in data]
    c_w = [i.c.orientation.w for i in data]

    axis_vectors = np.array(
        (a_u, a_v, a_w, b_u, b_v, b_w, c_u, c_v, c_w), dtype=float
    ).T

    return axis_lengths, axis_vectors


def surface_areas_to_ndarray(surface_areas: rtt.EllipsoidSurfaceAreas) -> np.ndarray:
    """
    Convert surface area data type into a NumPy array.

    This function converts an EllipsoidSurfaceAreas object into a NumPy array
    containing the surface area values.

    Parameters
    ----------
    surface_areas : EllipsoidSurfaceAreas
        The EllipsoidSurfaceAreas object containing the surface area data.

    Returns
    -------
    np.ndarray
        A NumPy array containing the surface area values.

    Examples
    --------
    >>> surface_areas = EllipsoidSurfaceAreas(data=[
    ...     Area(value=100.0, unit_squared=Units.MICRON),
    ...     Area(value=200.0, unit_squared=Units.MICRON),
    ...     Area(value=300.0, unit_squared=Units.MICRON)
    ... ])
    >>> surface_areas_to_ndarray(surface_areas)
    array([100., 200., 300.])
    """

    data = surface_areas.data
    areas = [i.value for i in data]

    ndarray = np.array((areas), dtype=float).T

    return ndarray


def volumes_to_ndarray(ellipsoid_volumes: rtt.EllipsoidVolumes) -> np.ndarray:
    """
    Convert ellipsoid volume data type into a NumPy array.

    This function converts an EllipsoidVolumes object into a NumPy array
    containing the volume values.

    Parameters
    ----------
    ellipsoid_volumes : EllipsoidVolumes
        The EllipsoidVolumes object containing the volume data.

    Returns
    -------
    np.ndarray
        A NumPy array containing the volume values.

    Examples
    --------
    >>> ellipsoid_volumes = EllipsoidVolumes(data=[
    ...     Volume(value=500.0, unit_cubed=Units.MICRON),
    ...     Volume(value=1000.0, unit_cubed=Units.MICRON),
    ...     Volume(value=1500.0, unit_cubed=Units.MICRON)
    ... ])
    >>> volumes_to_ndarray(ellipsoid_volumes)
    array([ 500., 1000., 1500.])
    """

    data = ellipsoid_volumes.data
    volumes = [i.value for i in data]

    ndarray = np.array((volumes), dtype=float).T

    return ndarray


def rmdir(directory: Path) -> None:
    """
    Recursively delete a directory and all its contents.

    This function deletes the specified directory and all its contents, including
    subdirectories and files. If the directory does not exist, the function does nothing.

    Parameters
    ----------
    directory : Path
        The path to the directory to be deleted.

    Returns
    -------
    None

    Examples
    --------
    >>> from pathlib import Path
    >>> directory = Path("path/to/directory")
    >>> rmdir(directory)
    """

    shutil.rmtree(str(directory), ignore_errors=True)


def compare_files(file1: Path, file2: Path, ignore_words: list[str]) -> bool:
    """
    Compare two files line-by-line, ignoring lines that contain specified words.

    This function compares two files line-by-line and ignores lines that contain
    any of the words specified in the `ignore_words` list. If the files are identical
    up to the allowed differences, the function returns True. Otherwise, it returns False.

    Parameters
    ----------
    file1 : Path
        Path to the first file for comparison.
    file2 : Path
        Path to the second file for comparison.
    ignore_words : list of str
        List of words that cause a particular line to be ignored during comparison.
        Use an empty list for strict word-by-word comparison. Use a non-empty list,
        e.g., ["date", "time"], to ignore lines with these words, as they may differ
        from file to file.

    Returns
    -------
    bool
        True if the two files are the same (up to any allowed differences in the
        ignore_words list), False if the two files are different.

    Examples
    --------
    >>> file1 = Path("path/to/file1.txt")
    >>> file2 = Path("path/to/file2.txt")
    >>> ignore_words = ["date", "time"]
    >>> compare_files(file1, file2, ignore_words)
    True
    """

    with open(file=file1, mode="r", encoding="utf-8") as fa:
        with open(file=file2, mode="r", encoding="utf-8") as fb:
            dba = fa.readlines()
            dbb = fb.readlines()

            if len(dba) == len(dbb):
                for k, line in enumerate(dba):
                    if line != dbb[k]:
                        # if "Autogenerated" in line:
                        iw_in_fa_line = [iw in line for iw in ignore_words]
                        iw_in_fb_line = [iw in dbb[k] for iw in ignore_words]
                        in_either_line = iw_in_fa_line + iw_in_fb_line
                        # if any(iw in [line + dbb[k]] for iw in ignore_words):
                        if any(in_either_line):
                            # skip comparison of the ignore_words case(s)
                            print("Skip comparison of word in ignore_word.")
                        else:
                            print(f"Files differ at line {k}:")
                            print(f"\t file1: {line}")
                            print(f"\t file2: {dbb[k]}")
                            return False

            else:
                print("Files have non-equal line length.")
                print("Line-by-line comparison not performed.")
                return False

    return True  # files are the same


def date_time_utc() -> str:
    """
    Returns the current date and time in ISO format in the UTC time zone.

    This function returns the current date and time in ISO 8601 format, with the time zone
    specified as UTC. The colons (:) and periods (.) in the time string are replaced with
    underscores (_) for compatibility with file naming conventions.

    Returns
    -------
    str
        The current date and time in ISO format in the UTC time zone.

    Examples
    --------
    >>> date_time_utc()
    '2024-08-13_UTC_12_34_56_789012'
    """

    # ts = datetime.utcnow().isoformat()  # The "naive object in UTC", utcnow is deprecated
    ts = datetime.now(datetime.UTC).isocalendar()  # The "naive object in UTC"
    ts = ts.replace(":", "_").replace(".", "_")  # Overwrite ":", "." with "_"
    ts = ts.replace("T", "_UTC_")  # Overwrite T with UTC time zone indication
    return ts


def underline() -> str:
    """
    Return a commented underline composed of 78 dashes ('-') and a newline character ('\\n').

    This function generates a string that represents a commented underline
    consisting of 78 dashes prefixed with a comment character ('* ') and
    followed by a newline character.

    Returns
    -------
    str
        A string representing a commented underline.

    Examples
    --------
    >>> underline()
    '* ------------------------------------------------------------------------------\\n'
    """

    return "* " + "".join(repeat("-", 78)) + "\n"


def in_a_but_not_in_b(*, a: Iterable, b: Iterable) -> Iterable:
    """
    Return items in Iterable `a` that are not in Iterable `b`.

    This function takes two iterables `a` and `b`, and returns a tuple containing
    all the items that are present in `a` but not in `b`.

    Parameters
    ----------
    a : Iterable
        The first iterable to compare.
    b : Iterable
        The second iterable to compare.

    Returns
    -------
    Tuple
        A tuple containing items that are in `a` but not in `b`.

    Examples
    --------
    >>> a = [1, 2, 3, 4]
    >>> b = [3, 4, 5, 6]
    >>> in_a_but_not_in_b(a=a, b=b)
    (1, 2)
    """

    result = ()  # empty tuple
    for item in a:
        if item not in b:
            result = result + (item,)

    return result


def pairwise(x: Iterable) -> Iterable:
    """
    Return successive overlapping pairs taken from the input iterable.
    This appears to be implemented in Python 3.10
    https://docs.python.org/3/library/itertools.html#itertools.pairwise
    but we currently use 3.9, so we implement `pairwise` here.

    The number of 2-tuples in the output iterator will be one fewer than the
    number of inputs. It will be empty if the input iterable has fewer than
    two values.

    Parameters
    ----------
    x : Iterable
        The input iterable from which to generate pairs.

    Returns
    -------
    Iterable[Tuple]
        An iterable of 2-tuples containing successive overlapping pairs from the input iterable.

    Examples
    --------
    >>> list(pairwise('ABCDEFG'))
    [('A', 'B'), ('B', 'C'), ('C', 'D'), ('D', 'E'), ('E', 'F'), ('F', 'G')]
    >>> list(pairwise([1, 2, 3, 4]))
    [(1, 2), (2, 3), (3, 4)]
    """

    a, b = tee(x)
    next(b, None)
    return zip(a, b)


def pairwise_circular(x: Iterable) -> Iterable:
    """
    Return successive overlapping pairs taken from the input iterable.

    The number of 2-tuples in the output iterator will be one fewer than the
    number of inputs. It will be empty if the input iterable has fewer than
    two values.

    Parameters
    ----------
    x : Iterable
        The input iterable from which to generate pairs.

    Returns
    -------
    Iterable[Tuple]
        An iterable of 2-tuples containing successive overlapping pairs from the input iterable.

    Examples
    --------
    >>> list(pairwise('ABCDEFG'))
    [('A', 'B'), ('B', 'C'), ('C', 'D'), ('D', 'E'), ('E', 'F'), ('F', 'G')]
    >>> list(pairwise([1, 2, 3, 4]))
    [(1, 2), (2, 3), (3, 4)]
    """

    a = cycle(x)
    next(a)
    return zip(x, a)


def metadata_to_dict(metadata: rtt.MetaData) -> dict:
    """
    Convert MetaData to a dictionary.

    This function converts a MetaData object to a dictionary representation.

    Parameters
    ----------
    metadata : MetaData
        The metadata of the Data Volume.

    Returns
    -------
    dict
        A dictionary containing the metadata.

    Examples
    --------
    >>> metadata = MetaData(
    ...     data_volume=DataVolume(x_width=256, y_height=256, z_image_count=10),
    ...     resolution=Resolution(
    ...         dx=Length(1.0, Units.MICRON),
    ...         dy=Length(1.0, Units.MICRON),
    ...         dz=Length(1.0, Units.MICRON)
    ...     ),
    ...     pixel_units=Units.MICRON,
    ...     origin=Origin(
    ...         x0=Length(0.0, Units.MICRON),
    ...         y0=Length(0.0, Units.MICRON),
    ...         z0=Length(0.0, Units.MICRON)
    ...     )
    ... )
    >>> metadata_to_dict(metadata)
    {
        'x_width (pixels)': '256',
        'y_height (pixels)': '256',
        'z_image_count (pixels)': '10',
        'Resolution, dx': '1.0',
        'Resolution, dy': '1.0',
        'Resolution, dz': '1.0',
        'Resolution, units': 'micron',
        'Origin, x0': '0.0',
        'Origin, y0': '0.0',
        'Origin, z0': '0.0',
        'Origin, units': 'micron'
    }
    """

    data_volume = metadata.data_volume
    resolution = metadata.resolution
    # pixel_units = metadata.pixel_units
    origin = metadata.origin

    x_width = data_volume.x_width
    y_height = data_volume.y_height
    z_image_count = data_volume.z_image_count

    dx = resolution.dx
    dy = resolution.dz
    dz = resolution.dz

    meta_dict = {
        "x_width (pixels)": str(x_width),
        "y_height (pixels)": str(y_height),
        "z_image_count (pixels)": str(z_image_count),
        "Resolution, dx": str(dx.value),
        "Resolution, dy": str(dy.value),
        "Resolution, dz": str(dz.value),
        "Resolution, units": str(dx.unit.value),
        "Origin, x0": str(origin.x0.value),
        "Origin, y0": str(origin.y0.value),
        "Origin, z0": str(origin.z0.value),
        "Origin, units": str(origin.x0.unit.value),
    }
    # print(meta_dict)

    return meta_dict


def yaml_to_dict(path_file_input: Path) -> dict:
    """
    Convert a YAML file to a dictionary.

    This function reads a YAML file and converts its contents to a dictionary.

    Parameters
    ----------
    path_file_input : Path
        The full path to the input YAML file.

    Returns
    -------
    dict
        A dictionary containing the contents of the YAML file.

    Raises
    ------
    TypeError
        If the file type is not supported.
    OSError
        If there is an error opening or decoding the YAML file.

    Examples
    --------
    >>> path_file_input = Path("path/to/config.yml")
    >>> yaml_to_dict(path_file_input)
    {'key1': 'value1', 'key2': 'value2', ...}
    """

    file_type = Path(path_file_input).suffix.casefold()
    supported_types = (".yaml", ".yml")

    if file_type not in supported_types:
        raise TypeError("Only file types .yaml, and .yml are supported.")

    try:
        with open(file=path_file_input, mode="r", encoding="utf-8") as stream:
            db = yaml.load(stream, Loader=yaml.SafeLoader)  # Load YAML file
    except yaml.YAMLError as error:
        print(f"Error with YAML file: {error}")
        print(f"Could not open or decode: {path_file_input}")
        raise OSError from error

    print(f"Success: database created from file: {path_file_input}")
    print("key, value, type")
    print("---, -----, ----")
    for key, value in db.items():
        print(f"{key}, {value}, {type(value)}")

    return db


def dict_to_yaml(db: dict, file: str) -> Path:
    """
    Convert a dictionary to a YAML file.

    This function writes the contents of a dictionary to a YAML file.

    Parameters
    ----------
    db : dict
        The dictionary to be converted to YAML.
    file : str
        The path to the output YAML file.

    Returns
    -------
    Path
        The path to the created YAML file.

    Examples
    --------
    >>> db = {'key1': 'value1', 'key2': 'value2'}
    >>> file = "path/to/output.yml"
    >>> dict_to_yaml(db, file)
    PosixPath('path/to/output.yml')
    """

    with open(file, "w", encoding="utf-8") as out_file:
        yaml.dump(
            db,
            out_file,
            default_flow_style=False,
            sort_keys=False,
        )  # Write dictionary to YAML

    return Path(file)


def ndarray_to_img(
    *,
    data: np.ndarray,
    slice_axis: rtt.CartesianAxis3D,
    parent_dir: Path,
    folder_name: str,
    pad_length: int = 4,
    file_type: str = ".tif",
) -> bool:
    """
    Convert an ndarray to an image stack and save it to a specified directory.

    This function takes a NumPy ndarray and creates an image stack, saving the images
    into a user-specified directory.

    Parameters
    ----------
    data : np.ndarray
        The semantic labels.
    slice_axis : CartesianAxis3D
        The axis along which to slice the ndarray:
        - 0 for Z axis
        - 1 for Y axis
        - 2 for X axis
    parent_dir : Path
        The parent directory to save the image folder.
    folder_name : str
        The folder name to save images.
    pad_length : int, optional
        The number of digits to pad the file names with (default is 4).
    file_type : str, optional
        The image file type (default is ".tif").

    Returns
    -------
    bool
        True if the images were successfully created.

    Examples
    --------
    >>> data = np.random.randint(0, 255, (10, 256, 256), dtype=np.uint8)
    >>> slice_axis = CartesianAxis3D.Z
    >>> parent_dir = Path("path/to/save")
    >>> folder_name = "image_stack"
    >>> ndarray_to_img(data=data, slice_axis=slice_axis, parent_dir=parent_dir, folder_name=folder_name)
    True
    """

    img_dir = parent_dir.joinpath(folder_name).expanduser()
    img_dir.mkdir(parents=True, exist_ok=True)
    n_slices = data.shape[slice_axis.value]

    for i in range(n_slices):
        fname = f"{img_dir}/{i:0{pad_length}}{file_type}"

        # slice out a [H,W] or [H,W,C] array
        match slice_axis:
            case rtt.CartesianAxis3D.Z:
                arr = data[i, ...]
            case rtt.CartesianAxis3D.Y:
                arr = data[:, i, ...]
            case rtt.CartesianAxis3D.X:
                arr = data[:, :, i]
            case _:
                raise ValueError(f"Bad slice_axis {slice_axis}")

        # if it's [H,W] we leave it
        # if it's [H,W,1] we drop the trailing channel
        # if it's [H,W,C] with C=3 or 4 we leave it
        if arr.ndim == 3 and arr.shape[2] == 1:
            arr = arr[:, :, 0]

        # optional: if you want to force every slice to RGB or L
        # uncomment the next two lines:
        # mode = "RGB" if (arr.ndim == 3 and arr.shape[2] == 3) else "L"
        # img = Image.fromarray(arr.astype("uint8"), mode=mode)

        # otherwise just let Pillow pick a mode for you
        img = Image.fromarray(arr)

        img.save(fname)

    return True


def read_images(
    file_dir: Path,
    file_type: str = ".tif",
    convert_mode: Optional[str] = None,
) -> np.ndarray:
    """
    Read all images of a given suffix from a directory, inspect their PIL mode,
    optionally convert them to a uniform mode, and return a single
    NumPy array of shape (N, H, W, C).

    Parameters
    ----------
    file_dir : Path
        Directory containing your images.
    file_type : str, optional
        File extension (default=".tif").
    convert_mode : str or None, optional
        If not None, will do `img = img.convert(convert_mode)` on every image.
        E.g. "L" for 8-bit grayscale, "RGB" for 8-bit 3-channel.

    Returns
    -------
    np.ndarray
        An array of shape (N, H, W, C), where C=1 for grayscale or C=3 for RGB.
        dtype comes straight from PIL (uint8, uint16 or float32).

    Raises
    ------
    FileNotFoundError
        If no files matching `*{file_type}` are found in `file_dir`.
    """

    # 1) find & sort
    files = sorted(glob.glob(f"{file_dir.as_posix()}/*{file_type}"))
    if not files:
        raise FileNotFoundError(f"No '*{file_type}' files found in '{file_dir}'")

    img_arrays = []
    first_file = True
    for p in files:
        with Image.open(p) as img:
            if first_file:
                print(
                    f"Reading {p!r}  mode={img.mode!r}  size={img.size}  info={img.info}"
                )

            # 2) optional uniform conversion
            if convert_mode is not None:
                img = img.convert(convert_mode)
                if first_file:
                    print(f"  -> converted to mode={img.mode!r}")

            first_file = False

            # 3) pull out numpy array (PIL will give you uint8 for "L"/"RGB",
            #    uint16 for "I;16" TIFF, float32 for "F", etc.)
            arr = np.array(img)

            # 4) if we got (H, W) grayscale, make it (H, W, 1)
            if arr.ndim == 2:
                arr = arr[..., None]

            # 5) if we got RGBA, drop alpha to force shape (H, W, 3)
            if arr.ndim == 3 and arr.shape[2] == 4:
                arr = arr[..., :3]

            img_arrays.append(arr)

    # 6) stack them into one array
    stack = np.stack(img_arrays, axis=0)

    print(
        f"Stacked {len(img_arrays)} images -> "
        f"shape={stack.shape}, dtype={stack.dtype}"
    )
    return stack
