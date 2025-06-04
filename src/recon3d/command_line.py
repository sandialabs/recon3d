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

# ANSI escape codes for formatting
BOLD = "\033[1m"
UNDERLINE = "\033[4m"
RESET = "\033[0m"
BLUE = "\033[34m"
DARK_GRAY = "\033[90m"

CLI_DOCS: Final[
    str
] = f"""
-------
recon3d
-------

{BOLD}{UNDERLINE}{BLUE}recon3d{RESET}

  (this command) Lists the recon3d command line entry points

{BOLD}{UNDERLINE}{BLUE}binary_to_semantic{RESET} <path_to_file>.yml

  Converts binary image stack to semantic image stack in a
  folder specified in the user input .yml file.

  Example:
  {DARK_GRAY}# Edit path variables in
  # ~/recon3d/docs/userguide/src/binary_to_semantic/binary_to_semantic.yml{RESET}
  (.venv) recon3d> binary_to_semantic binary_to_semantic.yml

{BOLD}{UNDERLINE}{BLUE}downscale{RESET} <path_to_file>.yml

  Downscales images in a folder specified in the user input .yml file.

  Example:
  {DARK_GRAY}# Edit path variables in
  # ~/recon3d/docs/userguide/src/downscale/downscale_thunder.yml{RESET}
  (.venv) recon3d> downscale downscale_thunder.yml

{BOLD}{UNDERLINE}{BLUE}grayscale_image_stack_to_segmentation{RESET} <path_to_file>.yml

  Converts a series of grayscale images to a segmentation.

  Example:
  {DARK_GRAY}# Edit path variables in
  # ~/recon3d/docs/userguide/src/utilities/grayscale_image_stack_to_segmentation.yml{RESET}
  (.venv) recon3d> grayscale_image_stack_to_segmentation grayscale_image_stack_to_segmentation.yml

{BOLD}{UNDERLINE}{BLUE}hello{RESET}

  Prints 'Hello world!' to the terminal to illustrate a command line entry point.

{BOLD}{UNDERLINE}{BLUE}hdf_to_image{RESET} <path_to_file>.yml

  From a dataset contained within a .hdf file specified by the input
  .yml file, creates an image stack with the same dataset name in
  the specified parent output folder.

  Example:
  {DARK_GRAY}# Edit path variables in
  # ~/recon3d/docs/userguide/src/hdf_to_image/hdf_to_image.yml{RESET}
  (.venv) recon3d> hdf_to_image hdf_to_image.yml

{BOLD}{UNDERLINE}{BLUE}hdf_to_npy{RESET} <path_to_file>.yml

  From a dataset contained within a .hdf file specified by the input
  .yml file, creates a NumPy .npy file from the segmentation data.

  Example:
  {DARK_GRAY}# Edit path variables in
  # ~/recon3d/docs/userguide/src/to_npy/hdf_to_npy.yml{RESET}
  (.venv) recon3d> hdf_to_npy hdf_to_npy.yml

{BOLD}{UNDERLINE}{BLUE}image_to_hdf{RESET} <path_to_file>.yml

  From a single image (or image stack) in a folder specified in the
  user input .yml file, creates a .hdf file in the specified
  output folder.

  Example:
  {DARK_GRAY}# Edit path variables in
  # ~/recon3d/docs/userguide/src/image_to_hdf/image_to_hdf.yml{RESET}
  (.venv) recon3d> image_to_hdf image_to_hdf.yml

{BOLD}{UNDERLINE}{BLUE}image_to_npy{RESET} <path_to_file>.yml

  From a series of images in a folder specified in the user input
  .yml file, creates a NumPy .npy file in the specified output folder.

  Example:
  {DARK_GRAY}# Edit path variables in
  # ~/recon3d/docs/userguide/src/to_npy/image_to_npy.yml{RESET}
  (.venv) recon3d> image_to_npy image_to_npy.yml

{BOLD}{UNDERLINE}{BLUE}instance_analysis{RESET} <path_to_file>.yml

  Digest a semantic segmentation accessible as a folder containing an image
  stack specified in the user input .yml file.

  Example:
  {DARK_GRAY}# Edit path variables in
  # ~/recon3d/docs/userguide/src/instance_analysis/instance_analysis.yml{RESET}
  (.venv) recon3d> instance_analysis instance_analysis.yml

{BOLD}{UNDERLINE}{BLUE}npy_to_mesh{RESET} <path_to_file>.yml

  Converts an instance or semantic segmentation, encoded as a .npy file,
  to an Exodus II finite element mesh using automesh.
  See https://autotwin.github.io/automesh/

  Example:
  {DARK_GRAY}# Edit path variables in
  # ~/recon3d/docs/userguide/src/npy_to_mesh/letter_f_3d.yml{RESET}
  (.venv) recon3d> npy_to_mesh letter_f_3d.yml

{BOLD}{UNDERLINE}{BLUE}semantic_to_binary{RESET} <path_to_file>.yml

  Converts semantic image stack to series of binary image stacks in
  a folder specified in the user input .yml file

  Example:
  {DARK_GRAY}# Edit path variables in
  # ~/recon3d/docs/userguide/src/binary_to_semantic/semantic_to_binary.yml{RESET}
  (.venv) recon3d> semantic_to_binary semantic_to_binary.yml

{BOLD}{UNDERLINE}{BLUE}void_descriptor{RESET} <path_to_file>.yml

  Work in progress, not yet implemented.
  From a pore dataset contained within a hdf file specified by
  the input .yml file, compute the void descriptor attributes
  for the void descriptor function.

  Example:
  {DARK_GRAY}# Edit path variables in
  # ~/recon3d/docs/userguide/src/void_descriptor/void_descriptor.yml{RESET}
  (.venv) recon3d> void_descriptor void_descriptor.yml
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
