"""The module test the image_stack_to_array module."""

# python standard libraries
from pathlib import Path
import numpy as np
from PIL import Image

# third-party libraries
import pytest

# local libraries
from recon3d import image_stack_to_array


def test_image_stack_to_array():
    """Test the image_stack_to_array on cylinder_machined_grayscale example"""

    cylinder_machined_grayscale_dir = Path(__file__).parent.joinpath(
        "data", "cylinder_machined_grayscale"
    )

    image_file_extension = ".tif"

    out_dir = cylinder_machined_grayscale_dir.parent

    image_stack_to_array.save_image_stack_to_array(
        cylinder_machined_grayscale_dir, image_file_extension, out_dir
    )

    out_filename = out_dir.joinpath(cylinder_machined_grayscale_dir.name + ".npy")
    # check file was written out

    assert Path.exists(out_filename)

    # check file contents is correct
    array_data = np.load(out_filename)

    slice_num = 10

    input_img = Path(__file__).parent.joinpath(
        "data", "cylinder_machined_grayscale", "0010.tif"
    )
    input_img_data = np.array(Image.open(input_img))

    assert np.allclose(array_data[slice_num, :, :], input_img_data)

    out_filename.unlink()  # clean up by removing the just-written file
