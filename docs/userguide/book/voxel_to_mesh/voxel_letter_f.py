"""This modulde demonstrates creating a pixel slice in the form of capital
letter `F`, and then extruding the slide three times to made a solid volume
voxel representation.

This module assume the virtual virtual environment has been loaded:

Example:

    cd ~/recon3d/
    source .venv/bin/activate
    cd ~/recon3d/docs/userguide/src/voxel_to_mesh/
    python voxel_letter_f.py
"""

# standard library
from pathlib import Path
from typing import Final

# third-party libary
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
import numpy as np
from PIL import Image

# module library
# none


def save_image(image, stack_name):
    """Saves an image to the file indicated by the fully pathed stack_name."""
    image.save(stack_name)
    print(f"Saved: {stack_name}")


def main():
    """The main program."""

    # user input begin
    output_dir: Final[str] = "~/scratch"
    # user input end

    # constants
    n_slices: Final[int] = 3  # The number of times "F" matrix is repeated
    n_color_bits: Final[int] = 255  # 8-bit color
    file_base: Final[str] = "/letter_f_slice_"
    file_ext: Final[str] = ".tif"
    dpi: Final[int] = 150  # resolution, dots per inch
    letter_f_pixel_path: Final[Path] = (
        Path(output_dir).expanduser().joinpath("letter_F_slice_pixel.png")
    )
    letter_f_voxel_path: Final[Path] = (
        Path(output_dir).expanduser().joinpath("letter_F_voxel.png")
    )
    # el, az, roll = 20, -15, -95
    el, az, roll = -160, -140, 0
    voxel_alpha: Final[float] = 0.9
    ls = LightSource(180, 45)
    invert_z_axis: Final[bool] = True

    # if the output directory does not already exist, create it
    output_path = Path(output_dir).expanduser()
    if not output_path.exists():
        print(f"Could not find existing output directory: {output_path}")
        Path.mkdir(output_path)
        print(f"Created: {output_path}")
        assert output_path.exists()

    letter_f_pixel = np.array(
        [[1, 1, 1], [1, 0, 0], [1, 1, 0], [1, 0, 0], [1, 0, 0]], dtype=np.uint8
    )

    print(f"The slice has the shape: {letter_f_pixel.shape}")
    fig, ax = plt.subplots()
    ax.set_xlabel("x (pixel)")
    ax.set_ylabel("y (pixel)")
    ax.imshow(letter_f_pixel)
    plt.show()
    fig.savefig(letter_f_pixel_path, dpi=dpi, bbox_inches="tight")
    print(f"Saved: {letter_f_pixel_path}")

    image = Image.fromarray(letter_f_pixel * n_color_bits)
    # end_cap = np.ones(letter_f_pixel.shape, dtype=np.uint8) * n_color_bits
    end_cap = np.ones(letter_f_pixel.shape, dtype=np.uint8)

    output_path = str(output_path)  # cast as string

    for i in range(1, n_slices + 1):  # base 1 index instead of 0
        # save images to the working folder
        stack_name = output_path + file_base + str(i) + file_ext
        save_image(image, stack_name)

    image = Image.fromarray(end_cap * n_color_bits)  # overwrite
    # the base 0 index is the for the first z-layer
    stack_name = output_path + file_base + str(0) + file_ext  # overwrite
    save_image(image, stack_name)

    letter_f_voxel = np.array(
        [
            end_cap,
            letter_f_pixel,
            letter_f_pixel,
            letter_f_pixel,
        ]
    )

    solid = letter_f_voxel == 1
    void = letter_f_voxel == 0

    fig2 = plt.figure()
    # ls = LightSource(45, 45)
    ax0 = fig2.add_subplot(121, projection="3d")
    # ax0.voxels(solid, facecolors="yellow", alpha=voxel_alpha)
    ax0.voxels(solid, facecolors="yellow", alpha=voxel_alpha, lightsource=ls)
    ax0.set_aspect("equal")
    ax0.view_init(elev=el, azim=az, roll=roll)
    ax0.set_title("solid")
    ax0.set_xlabel("z (voxel)")  # x and z are permutated
    ax0.set_ylabel("y (voxel)")
    ax0.set_zlabel("x (voxel)")  # x and z are permutated
    if invert_z_axis:
        print("inverting z-axis")
        ax0.invert_zaxis()
    z_ticks = [0, 1, 2, 3]
    y_ticks = [0, 1, 2, 3, 4]
    x_ticks = [0, 1, 2]
    ax0.set_xticks(z_ticks)  # x and z are permutated
    ax0.set_yticks(y_ticks)
    ax0.set_zticks(x_ticks)  # x and z are permutated

    ax1 = fig2.add_subplot(122, projection="3d")
    # ax1.voxels(void, facecolors="purple", alpha=voxel_alpha)
    ax1.voxels(void, facecolors="purple", alpha=voxel_alpha, lightsource=ls)
    ax1.set_aspect("equal")
    ax1.view_init(elev=el, azim=az, roll=roll)
    ax1.set_title("void")
    ax1.set_xlabel("z (voxel)")  # x and z are permutated
    ax1.set_ylabel("y (voxel)")
    ax1.set_zlabel("x (voxel)")  # x and z are permutated
    if invert_z_axis:
        print("inverting z-axis")
        ax1.invert_zaxis()
    ax1.set_xticks(z_ticks)  # x and z are permutated
    ax1.set_yticks(y_ticks)
    ax1.set_zticks(x_ticks)  # x and z are permutated
    plt.show()
    fig2.savefig(letter_f_voxel_path, dpi=dpi)
    print(f"Saved: {letter_f_voxel_path}")


if __name__ == "__main__":
    main()
