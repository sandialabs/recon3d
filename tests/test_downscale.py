"""This module tests individual functions for the downscale module."""

# python standard libraries
from pathlib import Path

# import shutil # unused import

# 3rd party libraries
import pytest

# import yaml

# local libraries
from recon3d import downscale as ds
import recon3d.utility as ut
from recon3d.static_test_paths import *

# @pytest.fixture
# def test_file() -> Path:
#     """The relative path and file string locating the default yml test file."""

#     return Path(__file__).parent.joinpath("examples", "downscale_example.yml")


def test_downscale():
    """Tests that the test images are downscaled."""
    aa = DOWNSCALE_YML
    db = ut.yaml_to_dict(aa)

    # overwrite so that this test runs on any machine
    db["image_dir"] = str(
        Path(__file__).parent.joinpath("data", "cylinder_machined_semantic")
    )
    db["out_dir"] = str(
        Path(__file__).parent.joinpath("test_output", "cylinder_machined")
    )

    temp_path = Path(__file__).parent.joinpath("test_output", "temp.yml")
    temp_path.parent.mkdir(
        parents=True, exist_ok=True
    )  # makes directory if it doesn't exist
    # make sure the recon3d/examples/temp.yml doesn't already
    # exist from a previous test run or from a test paused
    # midway through
    try:
        temp_path.unlink()  # delete the temp.yml if it exists aready
    except FileNotFoundError:
        pass

    bb = ut.dict_to_yaml(db, str(temp_path))
    success = ds.downscale(path_file_input=str(bb))

    ut.rmdir(temp_path.parent)
    # temp_path.parent.rmdir()  # delete the temp.yml after checking

    assert success


def test_apply_bbox():
    """
    Test the bounding box is working"""
    dataset_path = Path(__file__).parent.joinpath("data", "cylinder_machined_semantic")
    image_type = ".tif"

    image_stack = ut.read_images(dataset_path, image_type)

    threshold = 0
    bb_image_stack = ds.apply_bbox(image_stack, threshold)

    correct_shape = (21, 144, 145)

    assert bb_image_stack.shape == correct_shape


def test_bbox_range():
    """test it calculates the bounding box range correctly"""
    dataset_path = Path(__file__).parent.joinpath("data", "cylinder_machined_semantic")
    image_type = ".tif"

    image_stack = ut.read_images(dataset_path, image_type)
    threshold = 0
    found_bbox_range = ds.bbox_range(image_stack, threshold)

    known_bbox_range = (0, 20, 3, 146, 4, 148)

    assert found_bbox_range == known_bbox_range


def test_downscale_save_stack():
    """check the cropped images are saved"""
    dataset_path = Path(__file__).parent.joinpath("data", "cylinder_machined_semantic")
    image_type = ".tif"

    image_stack = ut.read_images(dataset_path, image_type)
    save_path = Path(__file__).parent.joinpath("test_output", "cylinder_machined")
    folder_suffix = f"{int(100.0)}_dx"

    ds.save_downscale_stack(image_stack, save_path, folder_suffix)

    known_output_folder_path = save_path.joinpath(
        f"images_at_resolution_{folder_suffix}"
    )
    found_image_stack = ut.read_images(known_output_folder_path, image_type)

    # # delete the output files after checking
    # try:
    #     shutil.rmtree(save_path)
    # except OSError as e:
    #     print(f"Failed to delete directory: {e}")

    ut.rmdir(Path(__file__).parent.joinpath("test_output"))

    assert found_image_stack.shape == image_stack.shape


# move to test_utility.py
#
# def test_dict_to_yaml():
#     """Tests that the dict is converted to a yaml file."""
#     temp_path = Path(__file__).parent.joinpath("files", "temp.yml")
#     db = {
#         "image_dir": None,
#         "image_type": ".tif",
#         "resolution_input": {"dx": 10.0, "dy": 10.0, "dz": 10.0},
#         "resolution_output": {"dx": 100.0, "dy": 100.0, "dz": 10.0},
#         "downscale_tolerance": 0.0001,
#         "image_limit_factor": 2.0,
#         "out_dir": None,
#         "padding": {"nx": 1, "ny": 1, "nz": 1},
#         "output_stack_type": "padded",
#         "save_npy": True,
#         "writeVTR": True,
#     }
#     temp_yaml = ds.dict_to_yaml(db, str(temp_path))
#     found_yaml = Path(__file__).parent.joinpath("files", "downscale_example.yml")

#     assert ut.compare_files(temp_yaml, found_yaml, [])


# def test_yaml_to_dict(test_file):
#     """Tests that the yaml file is converted to the dict."""
#     aa = test_file
#     known_db = {
#         "image_dir": None,
#         "image_type": ".tif",
#         "resolution_input": {"dx": 10.0, "dy": 10.0, "dz": 10.0},
#         "resolution_output": {"dx": 100.0, "dy": 100.0, "dz": 10.0},
#         "downscale_tolerance": 0.0001,
#         "image_limit_factor": 2.0,
#         "out_dir": None,
#         "padding": {"nx": 1, "ny": 1, "nz": 1},
#         "output_stack_type": "padded",
#         "save_npy": True,
#         "writeVTR": True,
#     }
#     found_db = ds.yaml_to_dict(aa)

#     assert known_db == found_db


def test_padded_dim_small():
    """test padding with a small image size"""
    img_stack_dim_size = 10
    target_res = 7.0
    original_res = 3.0
    tolerance = 0.01
    limit_factor = 2
    found_padded_size = ds.padded_size(
        img_stack_dim_size, target_res, original_res, tolerance, limit_factor
    )

    known_padded_size = 14

    assert found_padded_size == known_padded_size


def test_padded_dim_large():
    """test padding with a small image size"""
    img_stack_dim_size = 15353  # dim2 14192
    target_res = 7.56
    original_res = 1.08
    tolerance = 0.001
    limit_factor = 2
    found_padded_size = ds.padded_size(
        img_stack_dim_size, target_res, original_res, tolerance, limit_factor
    )

    known_padded_size = 15358

    assert found_padded_size == known_padded_size


def test_pad_amount():
    """test padding with a large image size"""
    img_stack_dim_size = 15353  # dim2 14192
    target_res = 7.56
    original_res = 1.08
    tolerance = 0.001
    limit_factor = 2
    found_pad_dim = ds.pad_amount(
        img_stack_dim_size, target_res, original_res, tolerance, limit_factor
    )

    known_pad_dim = (3, 2)  # known dimension, 15358 is divisible

    assert found_pad_dim == known_pad_dim
