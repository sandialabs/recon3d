"""
Command Line Entry Points Module
================================

This module provides command line entry points for various functions
and utilities.  It serves as the interface between the command line
and the underlying functionality of the application.
"""

from typing import Final

# import pkg_resources  # part of setup tools
from importlib.metadata import version
import recon3d.constants as cs

CLI_DOCS: Final[
    str
] = """
--------
recon3d
--------

recon3d
    (this command)

binary_to_semantic <path_to_file>.yml
    Converts binary image stack to semantic image stack in a
    folder specified in the user input .yml file.

    Example:
        # Edit path variables in
        # ~/recon3d/docs/userguide/src/binary_to_semantic/binary_to_semantic.yml
        (.venv) recon3d> binary_to_semantic binary_to_semantic.yml

downscale <path_to_file>.yml
    Downscales images in a folder specified in the user input .yml file.

    Example:
        # Edit path variables in
        # ~/recon3d/docs/userguide/src/downscale/downscale_thunder.yml
        (.venv) recon3d> downscale downscale_thunder.yml

grayscale_image_stack_to_segmentation <path_to_file>.yml
    Converts a series of grayscale images to a segmentation.

    Example:
        # Edit path variables in
        # ~/recon3d/docs/userguide/src/utilities/grayscale_image_stack_to_segmentation.yml
        (.venv) recon3d> grayscale_image_stack_to_segmentation grayscale_image_stack_to_segmentation.yml

hello
    Prints 'Hello world!' to the terminal to illustrate command line
    entry points.

image_to_voxel <path_to_file>.yml
    From a single image (or image stack) in a folder specified in the
    user input .yml file, creates a hdf file in the specified
    output folder.

    Example:
        # Edit path variables in
        # ~/recon3d/docs/userguide/src/voxel_to_image/image_to_voxel.yml
        (.venv) recon3d> image_to_voxel image_to_voxel.yml

image_stack_to_array <path_to_file>.yml
    From a series of images in a folder specified in the user input
    .yml file, creates a npy file in the specified output folder.

    Example:
        # Edit path variables in
        # ~/recon3d/docs/userguide/src/utilities/image_stack_to_array.yml
        (.venv) recon3d> image_stack_to_array image_stack_to_array.yml

instance_analysis <path_to_file>.yml
    Digest a semantic segmentation accessible as a folder containing an image
    stack specified in the user input .yml file.

    Example:
        # Edit path variables in
        # ~/recon3d/docs/userguide/src/instance_analysis/instance_analysis.yml
        (.venv) recon3d> instance_analysis instance_analysis.yml

semantic_to_binary <path_to_file>.yml
    Converts semantic image stack to series of binary image stacks in
    a folder specified in the user input .yml file

    Example:
        # Edit path variables in
        # ~/recon3d/docs/userguide/src/binary_to_semantic/semantic_to_binary.yml
        (.venv) recon3d> semantic_to_binary semantic_to_binary.yml

void_descriptor <path_to_file>.yml
    Work in progress.
    From a pore dataset contained within a hdf file specified by
    the input .yml file, compute the void descriptor attributes
    for the void descriptor function.

    Example:
        # Edit path variables in
        # ~/recon3d/docs/userguide/src/void_descriptor/void_descriptor.yml
        (.venv) recon3d> void_descriptor void_descriptor.yml

voxel_to_image <path_to_file>.yml
    From a dataset contained within a hdf file specified by the input
    .yml file, creates an image stack with the same dataset name in
    the specified parent output folder.

    Example:
        # Edit path variables in
        # ~/recon3d/docs/userguide/src/voxel_to_image/voxel_to_image.yml
        (.venv) recon3d> voxel_to_image voxel_to_image.yml

voxel_to_mesh <path_to_file>.yml
    Converts an instance or semantic segmentation, encoded as a .npy file,
    to an Exodus II finite element mesh using automesh.
    See https://autotwin.github.io/automesh/

    Example:
        # Edit path variables in
        # ~/recon3d/docs/userguide/src/voxel_to_mesh/letter_f_3d.yml
        (.venv) recon3d> voxel_to_mesh letter_f_3d.yml


"""


def recon3d():
    """
    Prints the command line documentation to the command window.

    This function prints the contents of the global variable `CLI_DOCS` to the
    command window. It is assumed that `CLI_DOCS` contains the necessary
    documentation in string format.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    print(CLI_DOCS)


def hello() -> str:
    """
    Simple example of a function hooked to a command line entry point.

    This function serves as an example of how to hook a function to a command
    line entry point. When called, it returns the "Hello world!" string.

    Parameters
    ----------
    None

    Returns
    -------
    str
        The canonical "Hello world!" string.
    """

    return "Hello world!"


def module_version():
    """
    Prints the module version and the yml_schema_version.

    This function retrieves the version of the module and the YAML schema
    version from the `Constants` class in the `cs` module. It prints these
    versions to the command window and returns the module version.

    Parameters
    ----------
    None

    Returns
    -------
    ver : str
        The version of the module.
    """

    # ver = pkg_resources.require("recon3d")[0].version
    module_name = cs.Constants().module_short_name
    ver = version(module_name)
    print(f"module version: {ver}")
    print(f"  yml_schema_version: {cs.Constants().YML_SCHEMA_VERSION}")
    return ver


def work_in_progress():
    """
    Prints the 'Work in Progress (WIP)' warning message.

    This function prints a warning message indicating that the function is a
    work in progress and has not yet been implemented.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    print("Warning: Work in progress (WIP), function not yet implemented.")
