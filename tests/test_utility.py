"""This module tests the utility module.
"""

# python standard libraries
from pathlib import Path
import numpy as np

# third-party libraries
import pytest

# local libraries
import recon3d.types as cs

# from recon3d import downscale as ds
import recon3d.utility as ut

from recon3d.static_test_paths import *

# @pytest.fixture
# def test_file() -> Path:
#     """The relative path and file string locating the default yml test file."""

#     return Path(__file__).parent.joinpath("examples", "downscale_example.yml")


# @pytest.fixture
# def test_dir() -> Path:
#     """Returns the relative path to the temp output dir for tests"""
#     return Path(__file__).parent.joinpath("files", "outputults")


# @pytest.fixture
# def test_binary_to_semantic_yml() -> Path:
#     """The relative path and file string locating the default yml test file."""

#     return Path(__file__).parent.joinpath("examples", "binary_to_semantic.yml")


# @pytest.fixture
# def test_semantic_to_binary_yml() -> Path:
#     """The relative path and file string locating the default yml test file."""

#     return Path(__file__).parent.joinpath("examples", "semantic_to_binary.yml")


@pytest.fixture
def test_data() -> Path:
    """The relative path and file string locating the default yml test file."""

    return Path(__file__).parent.joinpath("data", "cylinder_machined_semantic")


def create_binary_data(test_data) -> Path:
    """create binary data required for other tests

    Returns:
        path to folder of tif images

    """
    input_data = ut.read_images(test_data)

    outputdir_basepath = Path(__file__).parent.joinpath("data")
    outputdir_name = "binarized_cylinder_machined"

    bw_data = ut.binarize(input_data, 1)

    ut.ndarray_to_img(
        data=bw_data,
        slice_axis=cs.CartesianAxis3D.Z,
        parent_dir=outputdir_basepath,
        folder_name=outputdir_name,
    )

    return outputdir_basepath.joinpath(outputdir_name)


def test_binarize(test_data):
    """Test that semantic image data can be converted to binary image stack"""

    input_data = ut.read_images(test_data)

    # air = 0
    found_binary_data = ut.binarize(input_data, 0)
    found_unique, found_counts = np.unique(found_binary_data, return_counts=True)

    assert found_binary_data.shape == (21, 150, 150)
    assert np.allclose(found_unique, [0, 1])
    assert np.allclose(found_counts, [343817, 128683])

    # metal = 1
    found_binary_data = ut.binarize(input_data, 1)
    found_unique, found_counts = np.unique(found_binary_data, return_counts=True)

    assert found_binary_data.shape == (21, 150, 150)
    assert np.allclose(found_unique, [0, 1])
    assert np.allclose(found_counts, [129114, 343386])

    # pore = 2
    found_binary_data = ut.binarize(input_data, 2)
    found_unique, found_counts = np.unique(found_binary_data, return_counts=True)

    assert found_binary_data.shape == (21, 150, 150)
    assert np.allclose(found_unique, [0, 1])
    assert np.allclose(found_counts, [472069, 431])


def test_binary_with_pores_to_semantic(test_data):
    """Test that binary data can be converted to semantic image stack"""

    # Start with semantic stack, generate a binary stack of the metal
    binary_data_path = create_binary_data(test_data)

    # temporary directory for the semantic stack that will be generated
    semantic_data_path = binary_data_path.parent.joinpath("semantic_image_stack")

    # run function
    found_dict = ut.binary_with_pores_to_semantic(
        input_path=binary_data_path,
        output_path=semantic_data_path,
    )

    print(found_dict)

    found_data = ut.read_images(semantic_data_path)
    assert found_data.shape == (21, 150, 150)

    found_unique, found_counts = np.unique(found_data.data, return_counts=True)

    # check the existing semantic segmentation
    known_data = ut.read_images(test_data)
    known_unique, known_counts = np.unique(known_data, return_counts=True)

    assert np.allclose(found_unique, known_unique)
    assert np.allclose(found_counts, known_counts)

    # check with known value
    assert np.allclose(found_unique, [0, 1, 2])
    assert np.allclose(found_counts, [128683, 343386, 431])

    # clean up
    ut.rmdir(binary_data_path)
    ut.rmdir(semantic_data_path)


def test_binary_to_semantic():
    """Test that a yml can generate semantic image stack from binary image stack"""

    ut.binary_to_semantic(BINARY_TO_SEMANTIC_YML)

    yml_vals = ut.yaml_to_dict(BINARY_TO_SEMANTIC_YML)

    out_dir = Path(yml_vals["out_dir"]).expanduser()

    found_files = list(out_dir.iterdir())

    assert len(found_files) == 21

    ut.rmdir(out_dir)


def test_semantic_to_binary():
    """Test that a yml can generate semantic image stack from binary image stack"""

    ut.semantic_to_binary(SEMANTIC_TO_BINARY_YML)

    yml_vals = ut.yaml_to_dict(SEMANTIC_TO_BINARY_YML)

    out_dir = Path(yml_vals["out_dir"]).expanduser()

    found_files = list(out_dir.iterdir())

    assert len(found_files) == 21

    ut.rmdir(out_dir)


def test_dict_to_yaml():
    """Tests that the dict is converted to a yaml file."""
    temp_path = Path(__file__).parent.joinpath("data", "temp.yml")
    db = {
        "cli_entry_points": ["downscale"],
        "downscale_tolerance": 0.0001,
        "image_dir": None,
        "image_limit_factor": 2.0,
        "image_type": ".tif",
        "out_dir": None,
        "output_stack_type": "padded",
        "padding": {"nx": 1, "ny": 1, "nz": 1},
        "resolution_input": {"dx": 10.0, "dy": 10.0, "dz": 10.0},
        "resolution_output": {"dx": 100.0, "dy": 100.0, "dz": 10.0},
        "save_npy": True,
        "writeVTR": True,
    }
    temp_yaml = ut.dict_to_yaml(db, str(temp_path))
    found_yaml = DOWNSCALE_YML  # Path(__file__).parent.joinpath("examples", "downscale_example.yml")

    assert ut.compare_files(temp_yaml, found_yaml, [])

    try:
        temp_path.unlink()  # delete the temp.yml if it exists aready
    except FileNotFoundError:
        pass


def test_hdf_to_instance_properties():
    """Tests that hdf file can be read to construct InstanceProperties"""

    dataset_path = dataset_path = (
        Path(__file__).parent.joinpath("data", "machined_cylinder.h5").expanduser()
    )

    found_props = ut.hdf_to_instance_properties(
        hdf_path=dataset_path, group_path="pores"
    )

    # TODO more testing with assertions
    assert found_props.ellipsoids.data[4].a.orientation.u == pytest.approx(-0.159112797)
    assert found_props.centroids.data[2].cy.value == pytest.approx(108.8947)


def test_nd_to_img():
    """Tests that the images are saved from a ndarray"""

    # TODO check other packages will work with small letter_f example images
    dataset_path = (
        Path(__file__)
        .parent.joinpath("data", "cylinder_machined_semantic")
        .expanduser()
    )
    img_save_path = (
        Path(__file__).parent.joinpath("output_results", "temp_slices").expanduser()
    )

    # remove at beginning of test
    ut.rmdir(img_save_path)

    image_stack = ut.read_images(dataset_path)

    ut.ndarray_to_img(
        data=image_stack,
        slice_axis=cs.CartesianAxis3D.Z,
        parent_dir=img_save_path,
        folder_name="XY_slices",
    )

    # There should be 4 images written out
    assert len(list(img_save_path.rglob("*.tif"))) == 21

    # remove at end of test
    ut.rmdir(img_save_path)

    ut.ndarray_to_img(
        data=image_stack,
        slice_axis=cs.CartesianAxis3D.X,
        parent_dir=img_save_path,
        folder_name="YZ_slices",
    )

    # There should be 3 more images written out
    assert len(list(img_save_path.rglob("*.tif"))) == 150

    # remove at end of test
    ut.rmdir(img_save_path)

    ut.ndarray_to_img(
        data=image_stack,
        slice_axis=cs.CartesianAxis3D.Y,
        parent_dir=img_save_path,
        folder_name="XZ_slices",
    )

    # There should be 5 more images written out
    assert len(list(img_save_path.rglob("*.tif"))) == 150

    # remove at end of test
    ut.rmdir(img_save_path)


def test_read_images():
    """Tests the image stack is read in."""

    # Test with two images
    dataset_path = (
        Path(__file__)
        .parent.joinpath("data", "cylinder_machined_semantic")
        .expanduser()
    )
    image_type = ".tif"

    image_stack = ut.read_images(dataset_path, image_type)

    correct_shape = (21, 150, 150)

    assert image_stack.shape == correct_shape

    # Test with single image
    dataset_path = Path(__file__).parent.joinpath("data", "thunder_gray")

    image_stack = ut.read_images(dataset_path, image_type)

    correct_shape = (1, 348, 734)

    assert image_stack.shape == correct_shape


def test_difference():
    """Given two tuples with overlapping items, tests than the
    returned tuple is the difference of the first tuple less
    the second tuple.
    """
    x1 = (
        1,
        2,
        3,
        4,
    )
    x2 = (
        4,
        2,
        5,
        6,
    )

    y = (
        1,
        3,
    )

    fx = ut.in_a_but_not_in_b(a=x1, b=x2)

    assert y == fx


def test_pairwise():
    """Tests the pairwise function with a simple test of four items."""
    x = "ABCD"
    y = (("A", "B"), ("B", "C"), ("C", "D"))
    fx = ut.pairwise(x)
    fx2 = tuple(fx)
    assert y == fx2


def test_pairwise_circular():
    """Tests pairwise circular with a simple test of four items."""
    x = "ABCD"
    y = (("A", "B"), ("B", "C"), ("C", "D"), ("D", "A"))
    fx = ut.pairwise_circular(x)
    fx2 = tuple(fx)
    assert y == fx2


def test_yaml_to_dict():
    """Tests that the yaml file is converted to the dict."""
    aa = DOWNSCALE_YML
    known_db = {
        "cli_entry_points": ["downscale"],
        "downscale_tolerance": 0.0001,
        "image_dir": None,
        "image_limit_factor": 2.0,
        "image_type": ".tif",
        "out_dir": None,
        "output_stack_type": "padded",
        "padding": {"nx": 1, "ny": 1, "nz": 1},
        "resolution_input": {"dx": 10.0, "dy": 10.0, "dz": 10.0},
        "resolution_output": {"dx": 100.0, "dy": 100.0, "dz": 10.0},
        "save_npy": True,
        "writeVTR": True,
    }
    found_db = ut.yaml_to_dict(aa)

    assert known_db == found_db
