"""This module converts a stack of grayscale images into a segmentation
array."""

import argparse
from pathlib import Path
from typing import NamedTuple

import yaml


class Recipe(NamedTuple):
    """Defines the recipe to run this module directly from its Python API.

    Attributes
    ----------
    """
    # To come.


def grayscale_image_stack_to_segmentation(*, yml_input_file: Path) -> bool:
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
    breakpoint()


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

    grayscale_image_stack_to_segmentation(yml_input_file=input_file)


if __name__ == "__main__":
    main()
