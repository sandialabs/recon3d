"""
This module stores constant values in a single location for recon3d.

This module defines various constants and conversion factors that are used 
throughout the recon3d application. It provides a centralized location for 
these constants, making the code more maintainable and easier to understand.

Classes
-------
Constants
    Creates all constants used in this module.

Conversions
    Unit conversions.

SkimageConnectivity
    Connectivity level to be used with skimage.measure.label function in 3D.

CC3DConnectivity
    Connectivity level to be used with cc3d.connected_components function in 3D.

Attributes
----------
MODULE_NAME : str
    The name of the module.
"""

from typing import NamedTuple, Final
import numpy as np

MODULE_NAME: Final[str] = "recon3d"


class Constants(NamedTuple):
    """
    Creates all constants used in this module.

    This class defines a set of constants that are used throughout the module.
    It utilizes the `NamedTuple` from the `typing` module to create immutable
    constants.

    Attributes
    ----------
    YML_SCHEMA_VERSION : float
        The version of the YAML schema.
    module_short_name : str
        The short name of the module.
    module_prompt : str
        The prompt string for the module.
    KNUD_THOMSEN_FACTOR : float
        The Knud Thomsen factor for the surface area of a scalene ellipsoid.
    """

    YML_SCHEMA_VERSION: float = 1.0
    module_short_name: str = MODULE_NAME
    module_prompt: str = MODULE_NAME + ">"
    KNUD_THOMSEN_FACTOR: float = 1.6075  # for surface area of scalene ellipsoid

    # def __init__(self) -> None:
    #     self.YML_SCHEMA_VERSION = 1.0  # Good


class Conversions(NamedTuple):
    """
    Unit conversions.

    This class defines constants for converting between radians and degrees.
    It utilizes the `NamedTuple` from the `typing` module to create immutable
    constants.

    Attributes
    ----------
    RAD_TO_DEG : float
        Conversion factor from radians to degrees.
    DEG_TO_RAD : float
        Conversion factor from degrees to radians.
    """

    RAD_TO_DEG: float = 180 / np.pi
    DEG_TO_RAD: float = 1 / RAD_TO_DEG


class SkimageConnectivity(NamedTuple):
    """
    Connectivity level to be used with skimage.measure.label function in 3D.

    This class defines constants for different levels of connectivity that can
    be used with the `skimage.measure.label` function in 3D image processing.
    It utilizes the `NamedTuple` from the `typing` module to create immutable
    constants.

    Attributes
    ----------
    FACE_NEIGHBORS : int
        Connectivity level with 6 face neighbors.
    FACE_EDGE_NEIGHBORS : int
        Connectivity level with 6 face neighbors and 12 edge neighbors.
    FACE_EDGE_VERTEX_NEIGHBORS : int
        Connectivity level with 6 face neighbors, 12 edge neighbors, and 8 vertex neighbors.
    """

    FACE_NEIGHBORS = 1
    FACE_EDGE_NEIGHBORS = 2
    FACE_EDGE_VERTEX_NEIGHBORS = 3


class CC3DConnectivity(NamedTuple):
    """
    Connectivity level to be used with cc3d.connect_components function in 3D.

    This class defines constants for different levels of connectivity that can
    be used with the `cc3d.connected_components` function.
    It utilizes the `NamedTuple` from the `typing` module to create immutable
    constants.

    Attributes
    ----------
    FACE_NEIGHBORS : int
        Connectivity level with 6 face neighbors.
    FACE_EDGE_NEIGHBORS : int
        Connectivity level with 6 face neighbors and 12 edge neighbors.
    FACE_EDGE_VERTEX_NEIGHBORS : int
        Connectivity level with 6 face neighbors, 12 edge neighbors, and 8 vertex neighbors.
    """

    FACE_NEIGHBORS = 6
    FACE_EDGE_NEIGHBORS = 18
    FACE_EDGE_VERTEX_NEIGHBORS = 26
