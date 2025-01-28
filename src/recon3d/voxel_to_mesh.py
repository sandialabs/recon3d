"""Converts a semantic segmentation into a finite element mesh.
"""

import argparse

# from importlib.metadata import version
from pathlib import Path
from typing import NamedTuple

import numpy as np
import yaml

# import automesh as am
from automesh import Voxels


class AutomeshRecipe(NamedTuple):
    """
    Defines the receipe to run automesh directly from its Python API.

    Attributes
    ----------
    npy_input : Path
        The path to the numpy input file.
    output_file: Path
        The path to the mesh output file.
    remove: List[int]
        Voxel IDs to remove from the mesh [default: [0,]].
    scale_x : float
        The scaling factor along the x-axis [default: 1.0].
    scale_y : float
        The scaling factor along the y-axis [default: 1.0].
    scale_z : float
        The scaling factor along the z-axis [default: 1.0].
    translate_x : float
        The translation along the x-axis [default: 0.0].
    translate_y : float
        The translation along the y-axis [default: 0.0].
    translate_z : float
        The translation along the z-axis [default: 0.0].
    """

    npy_input: Path
    output_file: Path
    remove: list[int]
    scale_x: float = 1.0
    scale_y: float = 1.0
    scale_z: float = 1.0
    translate_x: float = 0.0
    translate_y: float = 0.0
    translate_z: float = 0.0


def validate_recipe(*, recipe: AutomeshRecipe) -> bool:
    """
    Validate the given recipe.

    Ensures that all values in the Recipe NamedTuple are valid.

    Parameters
    ----------
    recipe : AutomeshRecipe
        The Recipe NamedTuple populated with items from the user input
        .yml file.

    Returns
    -------
    bool
        True if the recipe is valid, False otherwise.

    Raises
    ------
    AssertionError
        If any of the validation checks fail.

    Examples
    --------
    >>> recipe = AutomeshRecipe(
    ...     npy_input=Path("/path/to/input.npy"),
    ...     output_file=Path("/path/to/output.inp"),
    ...     remove=[0,],
    ...     scale_x=1.0,
    ...     scale_y=1.0,
    ...     scale_z=1.0,
    ...     translate_x=0.0,
    ...     translate_y=0.0,
    ...     translate_z=0.0,
    ... )
    >>> validate_recipe(recipe=recipe)
    True
    """

    amr = recipe

    # Assure the .npy input file can be found
    assert amr.npy_input.is_file(), f"Cannot find {amr.npy_input}"

    # Assure the folder that will contain the output file can be found
    output_path = amr.output_file.parent
    assert output_path.is_dir(), f"Cannot find output path {output_path}"

    removes = amr.remove
    assert isinstance(removes, list), f"{removes} must be a list"
    all_ints = [isinstance(x, int) for x in removes]
    all_nonneg = [x >= 0 for x in removes]
    assert all(all_ints), f"{removes} must be a list of integers"
    assert all(
        all_nonneg
    ), f"{removes} must be a list of non-negative integers"

    assert amr.scale_x > 0.0, f"{amr.scale_x} must be > 0.0"
    assert amr.scale_y > 0.0, f"{amr.scale_y} must be > 0.0"
    assert amr.scale_z > 0.0, f"{amr.scale_z} must be > 0.0"

    return True


def voxel_to_mesh(*, yml_input_file: Path) -> bool:
    """
    Convert a .npy file specified in a yml recipe to a finite element mesh.

    Parameters
    ----------
    yml_input_file : Path
        The .yml recipe that specifies path variables.

    Returns
    -------
    True if successful, False otherwise.

    Raises
    ------
    FileNotFoundError
        If the input .yml file is not found.
    TypeError
        If the input file type is not supported.
    OSError
        If there is an error with the yml module.

    Examples
    --------
    >>> yml_input_file = Path("/path/to/recipe.yml")
    >>> voxel_to_mesh(yml_input_file=yml_input_file)
    0
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

    recipe = AutomeshRecipe(
        npy_input=Path(db["npy_input"]).expanduser(),
        output_file=Path(db["output_file"]).expanduser(),
        remove=db["remove"],
        scale_x=db["scale_x"],
        scale_y=db["scale_y"],
        scale_z=db["scale_z"],
        translate_x=db["translate_x"],
        translate_y=db["translate_y"],
        translate_z=db["translate_z"],
    )

    validate_recipe(recipe=recipe)

    amr = recipe

    # Run Sculpt
    # cc = [str(recipe.sculpt_binary), "-i", str(path_sculpt_i)]
    # result = subprocess.run(cc)

    # run automesh as a subprocess
    cc = "mesh "
    cc += f"-i {amr.npy_input} "
    cc += f"-o {amr.output_file} "
    for ri in amr.remove:
        cc += f"-r {ri} "
    cc += f"--xscale {amr.scale_x} "
    cc += f"--yscale {amr.scale_y} "
    cc += f"--zscale {amr.scale_z} "
    cc += f"--xtranslate {amr.translate_x} "
    cc += f"--ytranslate {amr.translate_y} "
    cc += f"--ztranslate {amr.translate_z}"

    print("Running automesh with the following command:")
    print(f"{cc}")

    # The automesh library stores .npy arrays as [x[y[z]]], whereas the
    # image processing convention we wish to use is [z[y[x]]].  So, save
    # temporary .npy file that is reordered, and feel that to automesh
    # instead of the original .npy file.
    aa = np.load(str(amr.npy_input))
    bb = amr.npy_input.stem + "_xyz.npy"
    # temp_npy = amr.npy_input.parent.joinpath("temp.npy"
    temp_npy = amr.npy_input.parent.joinpath(bb)
    print(f"Created temporary file in xyz order for automesh: {temp_npy}")
    cc = aa.transpose(2, 1, 0)
    # TODO: Check if the z-order needs to be reveresed based on the CT serial
    # section data.
    np.save(temp_npy, cc)

    # voxels = am.Voxels.from_npy
    # voxels = Voxels.from_npy(str(amr.npy_input))
    voxels = Voxels.from_npy(str(temp_npy))
    elements = voxels.as_finite_elements(
        remove=amr.remove,
        scale=[amr.scale_x, amr.scale_y, amr.scale_z],
        translate=[amr.translate_x, amr.translate_y, amr.translate_z],
    )
    elements.write_exo(str(amr.output_file))
    print(f"Wrote output file: {amr.output_file}")

    # clean up the temporary file:
    try:
        # Delete the file
        temp_npy.unlink()
        print(f"Temporary file successfully deleted: {temp_npy}")
    except FileNotFoundError:
        print(f"Temporary file does not exist: {temp_npy}")
    except PermissionError:
        print(f"Permission denied: Unable to delete tempoary file: {temp_npy}")
    except Exception as e:
        print(f"An error occurred: {e}")

    return True  # success


def main():
    """
    Runs the module from the command line.

    This function serves as the entry point for terminal-based access to the
    module.  It processes a YAML input file to a finite element mesh.

    Parameters
    ----------
    None

    Returns
    -------
    None

    Examples
    --------
    To run the module, use the following command in the terminal:
    $ npy_to_mesh path/to/input.yml
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="the .yml npy to mesh recipe")
    args = parser.parse_args()
    input_file = args.input_file
    input_file = Path(input_file).expanduser()

    voxel_to_mesh(yml_input_file=input_file)


if __name__ == "__main__":
    main()
