"""Converts a series of images to a stack array."""

import argparse
from pathlib import Path
from typing import NamedTuple

import numpy as np
import yaml

import recon3d.utility as ut


class Recipe(NamedTuple):
    """Defines the recipe to run this module directly from its Python API.

    Attributes
    ----------
    """

    image_dir: Path
    image_type: str
    out_dir: Path


def save_image_stack_to_array(image_dir: Path, file_extension: str, out_dir: Path):
    """
    Reads all images with a specified file extension from a directory,
    stores them in a numpy array, and saves the array to a .npy file in the
    specified directory.

    Parameters:
        image_dir: The directory from which to read images.
        file_extension: The file extension of the images to read.
        out_dir: The directory where the .npy file will be saved.
    """
    try:
        # Assuming 'ut.read_images' is a valid function that reads images and
        # returns a numpy array
        array_data = ut.read_images(image_dir, file_extension)

        # Correctly form the output file path
        out_file_path = out_dir.joinpath(f"{image_dir.name}.npy")

        # Save the numpy array to a .npy file
        np.save(out_file_path, array_data)

        print(f"Images from {image_dir} saved to {out_file_path}")

    except Exception as e:
        print(f"An error occurred: {e}")


def validate_recipe(*, recipe: Recipe) -> bool:
    """Validate the given recipe.

    Ensures that all values in the Recipe NamedTuple are valid.

    Parameters
    ----------
    recipe : Recipe
        The recipe to validate.

    Raises
    ------
    AssertionError
        If any of the recipe attributes are invalid.

    Examples
    --------
    >>> recipe = Recipe(
    ...     input_path=Path("path/to/input"),
    ...     input_file_type=".tif",
    ...     output_path=Path("path/to/output")
    ... )
    >>> validate_recipe(recipe)
    True
    """
    assert isinstance(recipe.image_dir, Path), "image_dir must be a Path object"
    assert recipe.image_dir.is_dir(), "image_dir must be a directory"
    assert isinstance(recipe.image_type, str), "image_type must be a string"
    assert recipe.image_type in [
        ".tif",
        ".tiff",
    ], "image_type must be .tif or .tiff"
    assert isinstance(recipe.out_dir, Path), "out_path must be a Path object"
    assert recipe.out_dir.is_dir(), "out_path must be a directory"

    return True


def image_to_stack_array(*, yml_input_file: Path) -> bool:
    """Converts a series of images to a stack array.

    Parameters
    ----------
    yml_input_file : Path
        The path to the input .yml file containing the image paths and other
        parameters.

    Returns
    -------
    True if the conversion is successful, False otherwise.

    Raises
    ------
    FileNotFoundError
        If the input .yml file is not found.
    TypeError
        If the input file type is not supported.
    OSError
        If there is an error with the yml module.
    """
    print(f"This is {Path(__file__).resolve()}")

    fin = yml_input_file.resolve().expanduser()

    print(f"Processing file: {fin}")

    if not fin.is_file():
        raise FileNotFoundError(f"File not found: {str(fin)}")

    file_type = fin.suffix.casefold()
    supported_types = (".yaml", ".yml")

    if file_type not in supported_types:
        raise TypeError("Only file types .yaml, and .yml are supported.")

    db = []

    try:
        with open(file=fin, mode="r", encoding="utf-8") as stream:
            db = yaml.load(stream, Loader=yaml.SafeLoader)  # overwrite
    except yaml.YAMLError as error:
        print(f"Error with yml module: {error}")
        print(f"Could not open or decode: {fin}")
        raise OSError from error

    print(f"Success: database created from file: {fin}")
    print(db)

    recipe = Recipe(
        image_dir=Path(db["image_dir"]).expanduser(),
        image_type=db["image_type"],
        out_dir=Path(db["out_dir"]).expanduser(),
    )

    validate_recipe(recipe=recipe)

    save_image_stack_to_array(
        image_dir=recipe.image_dir,
        file_extension=recipe.image_type,
        out_dir=recipe.out_dir,
    )

    return True  # success


def main():
    """
    Runs the module from the command line.

    This function serves as the entry point for terminal-based access to the
    module. It uses the argparse library to parse command-line arguments.

    Parameters
    ----------
    None

    Returns
    -------
    None

    Examples
    --------
    To run the module, use the following command in the terminal:
    $ image_to_stack_array path/to/input.yml
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="the .yml input file")
    args = parser.parse_args()
    input_file = args.input_file
    input_file = Path(input_file).expanduser()

    image_to_stack_array(yml_input_file=input_file)


if __name__ == "__main__":
    main()
