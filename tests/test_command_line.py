""""This module tests the command_line services.

Example:
    To run
    - activate the virutal environment, then
    $ cd ~/recond3d

    pytest (without code coverage)
    pytest tests/test_command_line.py -v  # -v is for verbose

    to run just a single test in this module, for example
    pytest recon3d/tests/test_command_line.py::test_hello_world -v
"""

# from recon3d import command_line as cl
from recon3d import command_line as cl
from recon3d import constants as cs
from importlib.metadata import version


def test_hello_world():
    """Tests that the string 'Hello world!' is returned."""
    print("test_hello_world...")
    known_string_fiducial = "Hello world!"

    found_string = cl.hello()

    assert found_string == known_string_fiducial


def test_version():
    """Tests that the module version is returned."""

    # TODO no package metadata at the moment
    # assert False
    known_version = version(cs.Constants().module_short_name)

    found_version = cl.module_version()
    assert found_version == known_version
