"""The module test the voxel_to_mesh module.
"""

# python standard libraries
from pathlib import Path

# third-party libraries
import pytest

# local libraries
from recon3d import voxel_to_mesh as vm
from recon3d import static_test_paths as test_paths


def test_letter_f_3d():
    """Tests that the letter_f_3d.npy example can be run."""
    path_yml = test_paths.LETTER_F_3D_YML
    path_exo = test_paths.LETTER_F_3D_EXO

    result = vm.voxel_to_mesh(yml_input_file=path_yml)

    assert result  # success

    # clean up temporary Exodus II file that was created in this test.
    try:
        # delete the file
        assert path_exo.is_file()
        path_exo.unlink()
    except FileNotFoundError:
        print(f"Temporary file does not exist: {path_exo}")
    except PermissionError:
        print(f"Permission denied: Unable to delete tempoary file: {path_exo}")
    except Exception as e:
        print(f"An error occurred: {e}")
