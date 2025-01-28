"""
HDF5 Processing Module
=======================

This module provides a set of functions for creating, modifying, and interacting with HDF5 files,
as well as converting between image stacks and voxel data. It includes command line interfaces
for some of the key functionalities.

Functions
---------
create_hdf(hdf_path, base_container)
    Prescribe foundational structure of the HDF5 file.

modify_hdf_dataset(hdf_path, dataset_loc, data, dtype, operation, mutable)
    Modify an existing HDF5 file to alter data within a container.

write_attr_dict(hdf_path, d, dataset_loc)
    Write attributes to a specified dataset in the HDF5 file.

write_h5(hdf_path, data)
    Write data to an HDF5 file.

add_to_h5(data, hdf_path, hdf_group)
    Add data to an HDF5 file based on its type.

image_to_voxel(yml_path)
    Populate the HDF5 file with the semantic segmentation image stack specified in the YAML file, including metadata.

image_to_voxel_command_line()
    The command line wrapper for the `image_to_voxel` function.

voxel_to_image(yml_path)
    Save the image data within the HDF5 file as TIFFs in a new directory.

voxel_to_image_command_line()
    The command line wrapper for the `voxel_to_image` function.

Examples
--------
To create an HDF5 file with a base container:

    >>> create_hdf(Path("output.h5"), "base_group")

To modify an HDF5 dataset:

    >>> data = np.random.rand(10, 5)
    >>> modify_hdf_dataset(Path("output.h5"), "group/dataset", data, np.float64, "create", True)

To add data to an HDF5 file based on its type:

    >>> centroids = Centroids(...)
    >>> add_to_h5(centroids, Path("output.h5"), "group_name")

To populate an HDF5 file with image stack data from a YAML file:

    >>> yml_path = Path("config.yml")
    >>> hdf5_path = image_to_voxel(yml_path)

To convert voxel data in an HDF5 file to image files:

    >>> yml_path = Path("config.yml")
    >>> image_dir = voxel_to_image(yml_path)

To run the `image_to_voxel` function from the command line:

    $ python -m your_module_name input_file.yml

To run the `voxel_to_image` function from the command line:

    $ python -m your_module_name input_file.yml
"""

import argparse
from pathlib import Path
import numpy as np
import h5py
from functools import singledispatch

# from recon3d.feature_analysis import SemanticImageStack
from recon3d.types import *
import recon3d.utility as ut
import recon3d.instance_analysis as ia
import recon3d.types as cs


### BASE/HELPER FUNCTIONS ###
def create_hdf(hdf_path: Path, base_container: str) -> bool:
    """
    Prescribe foundational structure of the HDF5 file.

    This function creates an HDF5 file at the specified path and initializes it with a base container group.

    Parameters
    ----------
    hdf_path : Path
        The path to the location of the HDF5 file.
    base_container : str
        The name of the base container group in the HDF5 file.

    Returns
    -------
    bool
        True if the HDF5 file is created successfully, False otherwise.

    Examples
    --------
    >>> create_hdf(Path("output.h5"), "base_group")
    True
    """

    with h5py.File(hdf_path, "w") as file:
        file.create_group(base_container)


def modify_hdf_dataset(
    hdf_path: Path,
    dataset_loc: str,
    data: np.ndarray,
    dtype: type,
    operation: str,
    mutable: bool,
) -> bool:
    """
    Modify an existing HDF5 file to alter data within a container.

    This function modifies an HDF5 file by creating, appending, or overwriting a dataset at the specified location.

    Parameters
    ----------
    hdf_path : Path
        The path to the location of the HDF5 file.
    dataset_loc : str
        The internal path to the dataset in the HDF5 file (e.g., "file/container/container2/dataset").
    data : np.ndarray
        The array to be written to the HDF5 file at the specified dataset location.
    dtype : type
        The data type of the dataset (e.g., np.float64, np.int, np.uint16, for string: h5py.special_dtype(vlen=str)).
    operation : str
        The operation to perform on the dataset:
        - "create": create a new dataset
        - "append": append to an existing dataset along dimension 0
        - "overwrite": overwrite an existing dataset along dimension 0 (e.g., shrinking dataset)
    mutable : bool
        If True, the initial shape is zero in the first dimension, with no size limit.

    Returns
    -------
    bool
        True if the operation is successful, False otherwise.

    Raises
    ------
    ValueError
        If the dataset already exists when trying to create a new one.
    KeyError
        If an unsupported operation is requested.

    Examples
    --------
    >>> data = np.random.rand(10, 5)
    >>> modify_hdf_dataset(Path("output.h5"), "group/dataset", data, np.float64, "create", True)
    True

    >>> new_data = np.random.rand(5, 5)
    >>> modify_hdf_dataset(Path("output.h5"), "group/dataset", new_data, np.float64, "append", True)
    True

    >>> overwrite_data = np.random.rand(8, 5)
    >>> modify_hdf_dataset(Path("output.h5"), "group/dataset", overwrite_data, np.float64, "overwrite", True)
    True
    """

    # initial shape is zero in first dimension, with no size limit
    shape, maxshape = list(data.shape), list(data.shape)
    if mutable:
        # shape[0], maxshape[0] = 0, None
        maxshape[0] = None
    shape, maxshape = tuple(shape), tuple(maxshape)

    with h5py.File(hdf_path, "r+") as hdf_file:
        # with h5py.File(hdf_path, "r+", locking=False) as hdf_file:
        if operation == "create":
            if not dataset_loc in hdf_file:

                dataset = hdf_file.create_dataset(
                    name=dataset_loc,
                    data=data,
                    shape=shape,
                    dtype=dtype,
                    # compression="gzip", #TODO 4/2/24 compression not allowing read in HDF_View
                    maxshape=maxshape,
                    chunks=True,
                )
                return True
            else:
                raise ValueError(f'Dataset already exists, cannot "{operation}"')

        dataset = hdf_file[dataset_loc]

        if operation == "append":
            # appending to dataset
            dataset.resize(dataset.shape[0] + data.shape[0], axis=0)
            dataset[-data.shape[0] :] = data
            return True

        elif operation == "overwrite":
            # shrinking dataset
            dataset.resize(data.shape[0], axis=0)
            dataset[:] = data
            return True

        else:
            raise KeyError(
                f'operation of "{operation}" requested, only, {"create"}, {"expand"}, and {"overwrite"} currently supported.'
            )


def write_attr_dict(hdf_path: Path, d: dict, dataset_loc: str) -> bool:
    """Write attributes to group or dataset in the h5 file
    from the provided dictionary

    hdf_path : path object to location of hdf file
    d: dict of attributes to add
    dataset_loc: internal path to the dataset in the hdf file (e.g. "file/container/container2/dataset")

    """

    with h5py.File(hdf_path, "r+") as hdf_file:
        loc = hdf_file[dataset_loc]
        # print(d, type(d))
        for key, value in d.items():
            # print(key, value)
            if isinstance(value, Path):
                value = str(value)
            dt = h5py.special_dtype(vlen=str)

            loc.attrs.create(name=key, data=value, dtype=dt)

    return


### INITIAL DATA CREATION ###
def write_h5(hdf_path: Path, semantic_stack: SemanticImageStack) -> bool:
    """Write new h5 file
    the data will be written into a dataset
    the metadata will be written as attributes of the dataset

    hdf_path : path object to location of hdf file
    semantic_stack: the semantic image stack

    """

    base_container = "VoxelData"
    create_hdf(hdf_path=hdf_path, base_container=base_container)

    img_data = semantic_stack.data
    dtype = type(
        img_data.flat[0]
    )  # data.item(0) will return python dtype, which may be better
    dataset_loc = f"{base_container}/{semantic_stack.name}"
    modify_hdf_dataset(
        hdf_path=hdf_path,
        dataset_loc=dataset_loc,
        data=img_data,
        dtype=dtype,
        operation="create",
        mutable=False,
    )

    # TODO 4/2/2024 write metadata as attributes
    md = semantic_stack.metadata
    # TODO 5/7/2024 make metadata strings add to h5
    metadata_dict = ut.metadata_to_dict(md)
    # d = {"height": }
    write_attr_dict(hdf_path=hdf_path, d=metadata_dict, dataset_loc=dataset_loc)


### ALL OTHER FUNCTIONALITIES ###
@singledispatch
def add_to_h5(
    data,
    h5_path: Path,
    h5_group: str,  # internal "folder" in the hdf to save data
):
    """
    Add data to an HDF5 file.

    This is a generic function that adds data to an HDF5 file. Specific implementations
    are provided for different types of data.

    Parameters
    ----------
    data : object
        The data to be added to the HDF5 file.
    h5_path : Path
        The path to the HDF5 file.
    h5_group : str
        The internal "folder" in the HDF5 file to save data.

    Returns
    -------
    None

    Raises
    ------
    NotImplementedError
        If no handler is registered for the type of `data`.

    Examples
    --------
    >>> add_to_h5(some_data, Path("output.h5"), "group_name")
    NotImplementedError: No handler for type <class 'type'>
    """
    _ = data
    __ = h5_path
    ___ = h5_group
    raise NotImplementedError(f"No handler for type {type(data)}")


@add_to_h5.register(Centroids)
def _(
    data,
    h5_path: Path,
    h5_group: str,
):
    """
    Add centroid data to an HDF5 file.

    This function processes and adds centroid data to an HDF5 file.

    Parameters
    ----------
    data : Centroids
        The centroid data to be added.
    h5_path : Path
        The path to the HDF5 file.
    h5_group : str
        The internal "folder" in the HDF5 file to save data.

    Returns
    -------
    None

    Examples
    --------
    >>> centroids = Centroids(...)
    >>> add_to_h5(centroids, Path("output.h5"), "group_name")
    """

    units = data.data[0].cx.unit.value
    centroid_data = ut.centroids_to_ndarray(centroids=data)

    # add data
    dataset_loc = f"{h5_group}/centroids"
    modify_hdf_dataset(
        hdf_path=h5_path,
        dataset_loc=dataset_loc,
        data=centroid_data,
        dtype=float,
        operation="create",
        mutable=False,
    )
    attrs = {
        "nlabels": str(centroid_data.shape[0]),
        "units": str(units),
        "ordering": str("X, Y, Z"),
    }
    write_attr_dict(
        hdf_path=h5_path,
        d=attrs,
        dataset_loc=dataset_loc,
    )


@add_to_h5.register(BestFitEllipsoids)
def _(
    data,
    h5_path: Path,
    h5_group: str,
):
    """
    Add best fit ellipsoid data to an HDF5 file.

    This function processes and adds best fit ellipsoid data to an HDF5 file.

    Parameters
    ----------
    data : BestFitEllipsoids
        The best fit ellipsoid data to be added.
    h5_path : Path
        The path to the HDF5 file.
    h5_group : str
        The internal "folder" in the HDF5 file to save data.

    Returns
    -------
    None

    Examples
    --------
    >>> ellipsoids = BestFitEllipsoids(...)
    >>> add_to_h5(ellipsoids, Path("output.h5"), "group_name")
    """

    axis_units = data.data[0].a.length.unit.value
    axis_lengths, axis_vectors = ut.ellipsoids_to_ndarray(data)

    # add data
    dataset_loc = f"{h5_group}/semi-axis_lengths"
    modify_hdf_dataset(
        hdf_path=h5_path,
        dataset_loc=dataset_loc,
        data=axis_lengths,
        dtype=float,
        operation="create",
        mutable=False,
    )
    attrs = {
        "notes": str("property of best fit ellipsoid"),
        "nlabels": str(axis_lengths.shape[0]),
        "units": str(axis_units),
        "ordering": str("a, b, c, with a > b > c"),
    }
    write_attr_dict(
        hdf_path=h5_path,
        d=attrs,
        dataset_loc=dataset_loc,
    )

    # add orientation data
    dataset_loc = f"{h5_group}/axis_vectors"
    modify_hdf_dataset(
        hdf_path=h5_path,
        dataset_loc=dataset_loc,
        data=axis_vectors,
        dtype=float,
        operation="create",
        mutable=False,
    )
    attrs = {
        "notes": str("property of best fit ellipsoid"),
        "nlabels": str(axis_lengths.shape[0]),
        "units": "unit vector",
        "ordering": str(
            "u, v, w for each axis a, b, c \n\t(a_u, a_v, a_w, b_u, b_v, b_w, c_u, c_v, c_w)"
        ),
    }
    write_attr_dict(
        hdf_path=h5_path,
        d=attrs,
        dataset_loc=dataset_loc,
    )


@add_to_h5.register(EllipsoidSurfaceAreas)
def _(
    data,
    h5_path: Path,
    h5_group: str,
):
    """
    Add ellipsoid surface area data to an HDF5 file.

    This function processes and adds ellipsoid surface area data to an HDF5 file.

    Parameters
    ----------
    data : EllipsoidSurfaceAreas
        The ellipsoid surface area data to be added.
    h5_path : Path
        The path to the HDF5 file.
    h5_group : str
        The internal "folder" in the HDF5 file to save data.

    Returns
    -------
    None

    Examples
    --------
    >>> surface_areas = EllipsoidSurfaceAreas(...)
    >>> add_to_h5(surface_areas, Path("output.h5"), "group_name")
    """

    area_units = data.data[0].unit_squared.value
    surface_areas = ut.surface_areas_to_ndarray(data)

    # add data
    dataset_loc = f"{h5_group}/ellipsoid_surface_areas"
    modify_hdf_dataset(
        hdf_path=h5_path,
        dataset_loc=dataset_loc,
        data=surface_areas,
        dtype=float,
        operation="create",
        mutable=False,
    )
    attrs = {
        "notes": str("property of best fit ellipsoid"),
        "nlabels": str(surface_areas.shape[0]),
        "units_squared": f"{area_units}",
        "method": "Knud Thomsen approximation for scalene ellipsoids (2004)",
    }
    write_attr_dict(
        hdf_path=h5_path,
        d=attrs,
        dataset_loc=dataset_loc,
    )


@add_to_h5.register(EllipsoidVolumes)
def _(
    data,
    h5_path: Path,
    h5_group: str,
):
    """
    Add ellipsoid volume data to an HDF5 file.

    This function processes and adds ellipsoid volume data to an HDF5 file.

    Parameters
    ----------
    data : EllipsoidVolumes
        The ellipsoid volume data to be added.
    h5_path : Path
        The path to the HDF5 file.
    h5_group : str
        The internal "folder" in the HDF5 file to save data.

    Returns
    -------
    None

    Examples
    --------
    >>> volumes = EllipsoidVolumes(...)
    >>> add_to_h5(volumes, Path("output.h5"), "group_name")
    """
    area_units = data.data[0].unit_cubed.value
    ellipsoid_volumes = ut.volumes_to_ndarray(data)

    # add data
    dataset_loc = f"{h5_group}/ellipsoid_volumes"
    modify_hdf_dataset(
        hdf_path=h5_path,
        dataset_loc=dataset_loc,
        data=ellipsoid_volumes,
        dtype=float,
        operation="create",
        mutable=False,
    )
    attrs = {
        "notes": str("property of best fit ellipsoid"),
        "nlabels": str(ellipsoid_volumes.shape[0]),
        "units_cubed": f"{area_units}",
        "method": "4/3 * pi * a * b * c",
    }
    write_attr_dict(
        hdf_path=h5_path,
        d=attrs,
        dataset_loc=dataset_loc,
    )


@add_to_h5.register(InstanceImageStack)
def _(
    data,
    h5_path: Path,
    h5_group: str,
):
    """
    Add instance image stack data to an HDF5 file.

    This function processes and adds instance image stack data to an HDF5 file.

    Parameters
    ----------
    data : InstanceImageStack
        The instance image stack data to be added.
    h5_path : Path
        The path to the HDF5 file.
    h5_group : str
        The internal "folder" in the HDF5 file to save data.

    Returns
    -------
    None

    Examples
    --------
    >>> image_stack = InstanceImageStack(...)
    >>> add_to_h5(image_stack, Path("output.h5"), "group_name")
    """

    # print(f"data type: {type(data)}")

    instance_image_stack = data

    # add data
    dataset_loc = f"{h5_group}/{instance_image_stack.name}"
    modify_hdf_dataset(
        hdf_path=h5_path,
        dataset_loc=dataset_loc,
        data=instance_image_stack.data,
        dtype=type(instance_image_stack.data.flat[0]),
        operation="create",
        mutable=False,
    )

    # add metadata
    md = instance_image_stack.metadata
    metadata_dict = ut.metadata_to_dict(md)
    write_attr_dict(
        hdf_path=h5_path,
        d=metadata_dict,
        dataset_loc=dataset_loc,
    )

    extra_attrs = {
        "nlabels": str(instance_image_stack.nlabels),
        "min_feature_size": str(instance_image_stack.min_feature_size),
    }
    write_attr_dict(
        hdf_path=h5_path,
        d=extra_attrs,
        dataset_loc=dataset_loc,
    )


@add_to_h5.register(InstanceIndices)
def _(
    data,
    h5_path: Path,
    h5_group: str,
):
    """
    Add instance indices data to an HDF5 file.

    This function processes and adds instance indices data to an HDF5 file.

    Parameters
    ----------
    data : InstanceIndices
        The instance indices data to be added.
    h5_path : Path
        The path to the HDF5 file.
    h5_group : str
        The internal "folder" in the HDF5 file to save data.

    Returns
    -------
    None

    Examples
    --------
    >>> indices = InstanceIndices(...)
    >>> add_to_h5(indices, Path("output.h5"), "group_name")
    """
    raise NotImplementedError(
        "Ability to write out instance indices not yet implemented."
    )
    # TODO: AP, write out as single variable length data instead of individual datasets (currently)
    # https://docs.h5py.org/en/stable/special.html
    # dt = h5py.vlen_dtype(np.dtype('int32'))
    instance_indices = data

    base_container = f"{instance_indices.source_name}_indices"
    # NumPy doesn’t support ragged arrays, and the ‘arrays of arrays’
    # h5py uses as a workaround are not as convenient or efficient as
    # regular NumPy arrays. If you’re deciding how to store data,
    # consider whether there’s a sensible way to do it without a
    # variable-length type.
    ragged_label_indices_array = np.array((len(instance_indices.indices)))
    for each_label in instance_indices.labels.data:
        ragged_label_indices_array[each_label] = instance_indices.indices[each_label]

    dataset_loc = f"{base_container}/label_indices"
    modify_hdf_dataset(
        hdf_path=h5_path,
        dataset_loc=dataset_loc,
        data=ragged_label_indices_array,
        dtype=h5py.vlen_dtype(np.dtype("int32")),
        operation="create",
        mutable=False,
    )
    # # Tree the labels
    # for (
    #     each_label
    # ) in (
    #     instance_indices.labels.data
    # ):  # TODO: why does the instance_indices.labels not return the correct type, instead it provides the ndarray...

    #     dataset_loc = f"{base_container}/label_{each_label:06d}"
    #     modify_hdf_dataset(
    #         hdf_path=h5_path,
    #         dataset_loc=dataset_loc,
    #         data=instance_indices.indices[each_label],
    #         dtype=type(instance_indices.indices[each_label].flat[0]),
    #         operation="create",
    #         mutable=False,
    #     )


@add_to_h5.register(InstanceProperties)
def _(
    data,
    h5_path: Path,
    h5_group: str,
):
    """
    Add instance properties data to an HDF5 file.

    This function processes and adds instance properties data to an HDF5 file.

    Parameters
    ----------
    data : InstanceProperties
        The instance properties data to be added.
    h5_path : Path
        The path to the HDF5 file.
    h5_group : str
        The internal "folder" in the HDF5 file to save data.

    Returns
    -------
    None

    Examples
    --------
    >>> properties = InstanceProperties(...)
    >>> add_to_h5(properties, Path("output.h5"), "group_name")
    """

    # centroids
    add_to_h5(data.centroids, h5_path=h5_path, h5_group=h5_group)

    # ellipsoids
    add_to_h5(data.ellipsoids, h5_path=h5_path, h5_group=h5_group)

    # surface area
    add_to_h5(data.surface_areas, h5_path=h5_path, h5_group=h5_group)

    add_to_h5(data.volumes, h5_path=h5_path, h5_group=h5_group)

    # equiv spherical diameter
    eq_diam = data.equivalent_sphere_diameters
    diam = [i.value for i in eq_diam]
    diam_data = np.array(diam, dtype=float).T
    diam_units = eq_diam[0].unit.value

    dataset_loc = f"{h5_group}/equivalent_sphere_diameters"
    modify_hdf_dataset(
        hdf_path=h5_path,
        dataset_loc=dataset_loc,
        data=diam_data,
        dtype=float,
        operation="create",
        mutable=False,
    )
    attrs = {
        "notes": str("from volume detemined by voxel count and resolution"),
        "nlabels": str(diam_data.shape[0]),
        "units": str(diam_units),
    }
    write_attr_dict(
        hdf_path=h5_path,
        d=attrs,
        dataset_loc=dataset_loc,
    )

    # num voxels
    num_voxel = data.n_voxels
    n_voxels = [i.value for i in num_voxel]
    n_vox_data = np.array(n_voxels, dtype=int).T
    n_vox_units = num_voxel[0].unit.value

    dataset_loc = f"{h5_group}/num_voxels"
    modify_hdf_dataset(
        hdf_path=h5_path,
        dataset_loc=dataset_loc,
        data=n_vox_data,
        dtype=int,
        operation="create",
        mutable=False,
    )
    attrs = {
        "nlabels": str(diam_data.shape[0]),
        "units": str(n_vox_units),
    }
    write_attr_dict(
        hdf_path=h5_path,
        d=attrs,
        dataset_loc=dataset_loc,
    )


@add_to_h5.register(NthNearestNeighbors)
def _(
    data,
    h5_path: Path,
    h5_group: str,
):
    """
    Add nth nearest neighbor data to an HDF5 file.

    This function processes and adds nth nearest neighbor data to an HDF5 file.

    Parameters
    ----------
    data : NthNearestNeighbors
        The nth nearest neighbor data to be added.
    h5_path : Path
        The path to the HDF5 file.
    h5_group : str
        The internal "folder" in the HDF5 file to save data.

    Returns
    -------
    None

    Examples
    --------
    >>> neighbors = NthNearestNeighbors(...)
    >>> add_to_h5(neighbors, Path("output.h5"), "group_name")
    """

    aa = 2
    units = data.distances[0].unit.value

    distance_list = data.distances
    distances = [i.value for i in distance_list]
    distance_array = np.array(distances, dtype=float).T
    neighbor_id = data.instance_id

    dataset_loc = f"{h5_group}/nearest_neighbor_distances"
    modify_hdf_dataset(
        hdf_path=h5_path,
        dataset_loc=dataset_loc,
        data=distance_array,
        dtype=float,
        operation="create",
        mutable=False,
    )
    attrs = {
        "nlabels": str(distance_array.shape[0]),
        "units": str(units),
        "nth_nearest": f"{data.nth_nearest}: 1st nearest is the feature itself",
    }
    write_attr_dict(
        hdf_path=h5_path,
        d=attrs,
        dataset_loc=dataset_loc,
    )

    dataset_loc = f"{h5_group}/nearest_neighbor_IDs"
    modify_hdf_dataset(
        hdf_path=h5_path,
        dataset_loc=dataset_loc,
        data=neighbor_id,
        dtype=int,
        operation="create",
        mutable=False,
    )
    attrs = {
        "nlabels": str(neighbor_id.shape[0]),
        "nth_nearest": f"{data.nth_nearest}: 1st nearest is the feature itself",
    }
    write_attr_dict(
        hdf_path=h5_path,
        d=attrs,
        dataset_loc=dataset_loc,
    )


### CONVERT BETWEEN H5 <-> FOLDER OF IMAGES   ###
def image_to_voxel(yml_path: Path) -> Path:
    """
    Populate the HDF5 file with the semantic segmentation image stack specified in the YAML file, including metadata.

    This function reads the YAML file to obtain the necessary parameters, processes the image stack, and writes the resulting data to an HDF5 file.

    Parameters
    ----------
    yml_path : Path
        The path to the YAML file containing the configuration and parameters.

    Returns
    -------
    Path
        The path to the created HDF5 file.

    Raises
    ------
    ValueError
        If the "cli_entry_points" key in the YAML file does not contain "image_to_voxel".

    Examples
    --------
    >>> yml_path = Path("config.yml")
    >>> hdf5_path = image_to_voxel(yml_path)
    Wrote output file: /path/to/output.h5
    >>> print(hdf5_path)
    /path/to/output.h5
    """

    # Import images and metadata into sematic_stack
    yml_vals = ut.yaml_to_dict(yml_path)

    # check cli_entry_points is valid
    if "image_to_voxel" not in yml_vals["cli_entry_points"]:
        raise ValueError(
            f"""Error. Incorrect yml format. 
                This function requires the "cli_entry_points" key to contain "image_to_voxel", 
                but currently contains the following options: {yml_vals["cli_entry_points"]} """
        )
    semantic_stack_save = ia.process_image_stack(yml_path)

    # Write semantic_stack to hdf
    h5_name = yml_vals["h5_filename"] + ".h5"

    # path_file_output = Path(yml_vals["out_dir"]).joinpath(h5_name)
    path_file_output = Path(yml_vals["out_dir"]).expanduser().joinpath(h5_name)
    write_h5(path_file_output, semantic_stack_save)

    print(f"Wrote output file: {path_file_output}")

    return path_file_output


def image_to_voxel_command_line():
    """
    The command line wrapper for the `image_to_voxel` function.

    This function sets up the command line argument parser, parses the input arguments,
    and calls the `image_to_voxel` function with the provided YAML input file.

    Parameters
    ----------
    None

    Returns
    -------
    None

    Examples
    --------
    To run this function from the command line:

        $ python -m your_module_name input_file.yml
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="the .yml user input file")
    args = parser.parse_args()
    input_file = args.input_file

    image_to_voxel(yml_path=input_file)


def voxel_to_image(yml_path: Path) -> Path:
    """
    Save the image data within the HDF5 file as TIFFs in a new directory.

    This function reads the YAML file to obtain the necessary parameters, extracts the image data from the HDF5 file,
    and saves the images as TIFF files in a specified directory.

    Parameters
    ----------
    yml_path : Path
        The path to the YAML file containing the configuration and parameters.

    Returns
    -------
    Path
        The path to the directory containing the saved images.

    Raises
    ------
    ValueError
        If the "cli_entry_points" key in the YAML file does not contain "voxel_to_image".
        If the specified slicing direction is not valid.

    Examples
    --------
    >>> yml_path = Path("config.yml")
    >>> image_dir = voxel_to_image(yml_path)
    >>> print(image_dir)
    /path/to/output/images
    """

    yml_vals = ut.yaml_to_dict(yml_path)

    # check cli_entry_points is valid
    if "voxel_to_image" not in yml_vals["cli_entry_points"]:
        raise ValueError(
            f"""Error. Incorrect yml format. 
            This function requires the "cli_entry_points" key to contain "voxel_to_image", 
            but currently contains the following options: {yml_vals["cli_entry_points"]} """
        )

    hdf_path = Path(yml_vals["voxel_data_path"]).expanduser()
    # TODO add alternative to ingest npy data dir

    hdf_dataset_location = yml_vals["voxel_data_location"]
    output_image_dir = Path(yml_vals["image_parent_dir"]).expanduser()
    output_image_type = yml_vals["image_output_type"]

    slice_normal = yml_vals["image_slice_normal"]
    valid_slice_normal = set(item.name for item in cs.CartesianAxis3D)
    if slice_normal not in valid_slice_normal:
        raise ValueError(
            f"Error, '{slice_normal}' is not a valid slicing direction, accepted units are: {valid_slice_normal}"
        )
    slice_axis = cs.CartesianAxis3D[slice_normal]

    with h5py.File(hdf_path, "r") as f:
        data = np.squeeze(f[hdf_dataset_location][:])
    ut.ndarray_to_img(
        data=data,
        parent_dir=output_image_dir,
        folder_name=Path(hdf_dataset_location).stem,
        file_type=output_image_type,
        slice_axis=slice_axis,
    )

    return output_image_dir.joinpath(hdf_path.stem)


def voxel_to_image_command_line():
    """
    The command line wrapper for the `voxel_to_image` function.

    This function sets up the command line argument parser, parses the input arguments,
    and calls the `voxel_to_image` function with the provided YAML input file.

    Parameters
    ----------
    None

    Returns
    -------
    None

    Examples
    --------
    To run this function from the command line:

        $ python -m your_module_name input_file.yml
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="the .yml user input file")
    args = parser.parse_args()
    input_file = args.input_file

    new_path = voxel_to_image(yml_path=input_file)

    print(f"\nVoxel data extracted to the following relative directory:'{new_path}'")
