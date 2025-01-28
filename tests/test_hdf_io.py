"""Tests functionality provided by hdf_io.py."""

from pathlib import Path

# Python standard libarires
import glob
from pathlib import Path

# Third party libraries
import numpy as np
import pytest
import h5py

# local libraries
# import recon3d.feature_analysis as fa
import recon3d.constants as cs
import recon3d.instance_analysis as ia
import recon3d.hdf_io as hio
import recon3d.utility as ut
from recon3d.types import *

from tests.test_instance_analysis import create_instance_stack
from recon3d.static_test_paths import *


# @pytest.fixture
# def test_letter_f() -> Path:
#     """Returns the relative path to the letter_f example .tif files."""
#     return Path(__file__).parent.joinpath("data", "letter_f")


# @pytest.fixture
# def test_instance_analysis_yml() -> Path:
#     """The relative path and file string locating the default yml test file."""

#     return Path(__file__).parent.joinpath("examples", "instance_analysis.yml")


### BASE/HELPER FUNCTIONS ###
def test_create_hdf():
    """Test hdf file is created from create_hdf function"""

    hdf_path = Path(__file__).parent.joinpath("data", "test_hdf_io.h5")

    base_container = "test_container"

    hio.create_hdf(hdf_path, base_container)

    assert hdf_path.is_file()

    # remove at end of test
    hdf_path.unlink()


def test_modify_hdf():
    """Test hdf file is created from modify_hdf function"""

    hdf_path = Path(__file__).parent.joinpath("data", "test_hdf_io.h5")

    try:
        hdf_path.unlink()  # delete the test_hdf_io.h5 if it exists aready
    except FileNotFoundError:
        pass

    base_container = "test_container"

    hio.create_hdf(hdf_path, base_container)

    dataset_loc = base_container + "/rand_test_data"
    test_data = np.random.rand(10, 30, 5)
    test_data_dtype = np.float64
    operation = "create"
    mutable = True
    hio.modify_hdf_dataset(
        hdf_path, dataset_loc, test_data, test_data_dtype, operation, mutable
    )

    with h5py.File(hdf_path, "r") as hdf_file:

        found_dataset = hdf_file[dataset_loc]
        assert found_dataset.shape == (10, 30, 5)
        assert found_dataset.dtype == test_data_dtype

    test_data = np.random.rand(12, 30, 5)
    test_data_dtype = np.float64
    operation = "append"
    mutable = True
    hio.modify_hdf_dataset(
        hdf_path, dataset_loc, test_data, test_data_dtype, operation, mutable
    )

    with h5py.File(hdf_path, "r") as hdf_file:

        found_dataset = hdf_file[dataset_loc]
        assert found_dataset.shape == (22, 30, 5)
        assert found_dataset.dtype == test_data_dtype

    operation = "overwrite"
    mutable = True
    hio.modify_hdf_dataset(
        hdf_path, dataset_loc, test_data, test_data_dtype, operation, mutable
    )

    with h5py.File(hdf_path, "r") as hdf_file:

        found_dataset = hdf_file[dataset_loc]
        assert found_dataset.shape == (12, 30, 5)
        assert found_dataset.dtype == test_data_dtype

    # remove at end of test
    hdf_path.unlink()


def test_write_attr_dict():
    """Test attrubutes of dataset in hdf file are added
    from write_attr_dict function"""

    hdf_path = Path(__file__).parent.joinpath("data", "test_hdf_io.h5")

    try:
        hdf_path.unlink()  # delete the test_hdf_io.h5 if it exists aready
    except FileNotFoundError:
        pass

    base_container = "test_container"

    hio.create_hdf(hdf_path, base_container)

    dataset_loc = base_container + "/rand_test_data"
    test_data = np.random.rand(10, 30, 5)
    test_data_dtype = np.float64
    operation = "create"
    mutable = True
    hio.modify_hdf_dataset(
        hdf_path, dataset_loc, test_data, test_data_dtype, operation, mutable
    )

    attr_dict = {"attr1": "dog", "attr2": "cat"}

    hio.write_attr_dict(hdf_path, attr_dict, dataset_loc)

    with h5py.File(hdf_path, "r") as hdf_file:

        found_attrs_keys = dict(hdf_file[dataset_loc].attrs.items())

        assert attr_dict == found_attrs_keys

    # remove at end of test
    hdf_path.unlink()


### ALL OTHER FUNCTIONALITIES ###
def test_add_to_h5():
    """Check the overloaded function can write the different types to the h5 correctly"""

    test_instance_analysis_yml = INSTANCE_ANALYSIS_YML

    output_voxeldata_path = hio.image_to_voxel(test_instance_analysis_yml)

    # InstanceImageStack
    instance_name = "air"
    instance_stack = create_instance_stack(instance_name)
    hio.add_to_h5(
        instance_stack,
        h5_path=output_voxeldata_path,
        h5_group="VoxelData",
    )

    instance_name = "metal"
    instance_stack = create_instance_stack(instance_name)
    hio.add_to_h5(
        instance_stack,
        h5_path=output_voxeldata_path,
        h5_group="VoxelData",
    )

    instance_name = "pore"
    instance_stack = create_instance_stack(instance_name)

    # get instance_indices
    instance_indices = ia.instance_indices(
        instance_stack=instance_stack,
    )
    if instance_stack.min_feature_size > 0:
        # overwrites instance_stack and inst_indices
        instance_stack, instance_indices = ia.minimum_size_filter(
            initial_stack=instance_stack,
            initial_indices=instance_indices,
        )
    hio.add_to_h5(
        instance_stack,
        h5_path=output_voxeldata_path,
        h5_group="VoxelData",
    )

    # Instance dtypes
    # instance_indices = ia.instance_indices(instance_stack)
    # # Centroids
    # centroid_1 = ia.center_of_mass(
    #     indices=instance_indices.indices[2],
    #     resolution=instance_stack.metadata.resolution,
    #     origin=instance_stack.metadata.origin,
    # )
    # centroid_2 = ia.center_of_mass(
    #     indices=instance_indices.indices[3],
    #     resolution=instance_stack.metadata.resolution,
    #     origin=instance_stack.metadata.origin,
    # )
    # centroids = Centroids(data=[centroid_1, centroid_2])
    # hio.add_to_h5(
    #     centroids,
    #     h5_path=output_voxeldata_path,
    #     h5_group="PoreData",
    # )
    # # Ellipsoids
    # ellipsoid_1 = ia.fit_ellipsoid(
    #     indices=instance_indices.indices[2],
    #     centroid=centroid_1,
    #     resolution=instance_stack.metadata.resolution,
    #     origin=instance_stack.metadata.origin,
    # )
    # ellipsoid_2 = ia.fit_ellipsoid(
    #     indices=instance_indices.indices[3],
    #     centroid=centroid_2,
    #     resolution=instance_stack.metadata.resolution,
    #     origin=instance_stack.metadata.origin,
    # )
    # ellipsoids = BestFitEllipsoids(data=[ellipsoid_1, ellipsoid_2])
    # hio.add_to_h5(
    #     ellipsoids,
    #     h5_path=output_voxeldata_path,
    #     h5_group="PoreData",
    # )

    instance_props = ia.instance_properties(
        instance_stack=instance_stack,
        inst_indices=instance_indices,
    )
    hio.add_to_h5(
        instance_props,
        h5_path=output_voxeldata_path,
        h5_group=instance_props.source_name,
    )

    nearest_neighbors = ia.nearest_neighbor_distance(centroids=instance_props.centroids)
    hio.add_to_h5(
        nearest_neighbors,
        h5_path=output_voxeldata_path,
        h5_group=instance_props.source_name,
    )

    # # InstanceIndices
    # instance_indices = ia.instance_indices(instance_stack)
    # hio.add_to_h5(instance_indices, h5_path=output_voxeldata_path)

    output_voxeldata_path.unlink()


### CONVERT BETWEEN H5 <-> FOLDER OF IMAGES   ###
def test_image_to_voxel():
    """Check the images can be properly added to the h5"""

    input_yml = VOXEL_TO_IMAGE_YML  # Path(__file__).parent.joinpath("examples", "voxel_to_image.yml")

    with pytest.raises(ValueError) as err:
        output_voxeldata_path = hio.image_to_voxel(input_yml)
    err_msg = 'Error. Incorrect yml format. \n                This function requires the "cli_entry_points" key to contain "image_to_voxel", \n                but currently contains the following options: [\'voxel_to_image\'] '
    assert err.type == ValueError
    assert err.value.args[0] == err_msg

    # overwrite
    input_yml = IMAGE_TO_VOXEL_YML  # Path(__file__).parent.joinpath("examples", "image_to_voxel.yml")

    output_voxeldata_path = hio.image_to_voxel(input_yml)

    with h5py.File(output_voxeldata_path, "r") as hdf_file:

        found_dataset = hdf_file["/VoxelData/letter_f_test"]
        assert found_dataset.shape == (4, 5, 3)

    # remove at end of test
    output_voxeldata_path.unlink()

    """Check the images from the h5 file can be properly extracted"""

    input_yml = IMAGE_TO_VOXEL_YML  # Path(__file__).parent.joinpath("examples", "image_to_voxel.yml")

    with pytest.raises(ValueError) as err:
        hio.voxel_to_image(input_yml)
    err_msg = 'Error. Incorrect yml format. \n            This function requires the "cli_entry_points" key to contain "voxel_to_image", \n            but currently contains the following options: [\'image_to_voxel\'] '
    assert err.type == ValueError
    assert err.value.args[0] == err_msg

    input_yml = VOXEL_TO_IMAGE_YML  # Path(__file__).parent.joinpath("examples", "voxel_to_image.yml")

    output_images_path = hio.voxel_to_image(input_yml)

    # There should be 4 images written out
    assert len(list(output_images_path.rglob("*.tif"))) == 21

    # remove at end of test
    ut.rmdir(output_images_path)
    ut.rmdir(output_images_path.parent)
