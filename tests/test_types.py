"""Tests that the types defined in the module are tested."""

import recon3d.types as tt


def test_origin():
    """Test the origin is at (0.0, 0.0, 0.0)."""
    length_null = tt.Length(value=0.0, unit=tt.Units.MICRON)
    origin = tt.Origin(x0=length_null, y0=length_null, z0=length_null)
    assert type(origin).__name__ == "Origin"
    assert origin.x0.value == 0.0
    assert origin.y0.value == 0.0
    assert origin.z0.value == 0.0


def test_centroid():
    """Tests that the Centroid data type is as expected."""
    centroid = tt.Centroid(
        cx=tt.Length(value=1.0, unit=tt.Units.MICRON),
        cy=tt.Length(value=2.0, unit=tt.Units.MICRON),
        cz=tt.Length(value=3.0, unit=tt.Units.MICRON),
    )
    assert type(centroid).__name__ == "Centroid"
    assert centroid.cx.value == 1.0
    assert centroid.cy.value == 2.0
    assert centroid.cz.value == 3.0


# Below is a history breadcrumb in case we need it for CI/CD.

"""A module for testing types.

Modules, classes, and methods can be
created for testing purposes.

"""

# import unittest
# from recon3d.core import Scan


# class TestTypes(unittest.TestCase):
#     """Class to test types."""

#     def main(self):
#         """Main function for module-level testing functionality."""

#         # TODO: Locate these implementations?  Comment out until found or implemented.
#         # self.test_int()  # where are these implementations?
#         # self.test_float()

#     def test_str(self):
#         """Function to test string representation."""

#         metadata = {"height": 2, "width": 3, "images": None, "xpixel_size": None}
#         self.assertEqual(type(Scan(metadata).__str__()), type("a string"))


# if __name__ == "__main__":
#     unittest.main()
