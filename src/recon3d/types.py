"""
A module to define types used by recon3d.

This module defines various custom types used throughout the recon3d application.
These types provide clear and consistent data structures for different components
of the application.

Types
-----
CartesianAxis2D
    Constant axis ordering in Python for 2D arrays (Y, X) ordering.
CartesianAxis3D
    Constant axis ordering in Python for 3D arrays (Z, Y, X) ordering.
DataVolume
    The integer counts of the data that compose a single image or a volume composed of a stack of images.
Units
    Unit enum for metadata.
UnitVector
    Normalized unit vector.
Angle
    Angle class with associated units.
Area
    Area class with associated units.
InstanceLabels
    Collection of instance labels.
Length
    Length class with associated units.
NthNearestNeighbors
    Nth nearest neighbors for instance analysis.
NVoxel
    Number of voxels, typically for an instance object.
Volume
    Volume class with associated units.
Centroid
    The centroid location for the X, Y, and Z axes, in length units.
Centroids
    Centroids for an instance image stack.
EllipsoidAxis
    Properties of best fit ellipsoid semi-axis.
EquivalentSphereDiameters
    Collection of spherical diameters.
Origin
    The origin in spatial coordinates.
Resolution
    The voxel resolution for the X, Y, and Z dimension in length units.
BestFitEllipsoid
    Three semi-axes of best-fit ellipsoid in descending order of size.
BestFitEllipsoids
    Collection of best-fit ellipsoids.
EllipsoidSurfaceAreas
    Collection of ellipsoid surface area objects.
EllipsoidVolumes
    Collection of ellipsoid volumes.
MetaData
    The data description of the scan.
ImageStack
    The generic image stack.
InstanceImageStack
    The image stack where each pixel in the image is categorized (labeled) into a specific instance of a class or object.
InstanceIndices
    Labels and index locations for each instance in an InstanceImageStack.
SemanticImageStack
    The image stack where each pixel in the image is categorized (labeled) into a class or object.
VoidDescriptorEllipsoidAngles
    Angles for ellipsoid properties to be used in void descriptor function analysis.
InstanceProperties
    Properties for labeled features in an instance image stack.
VoidDescriptorFittingParameters
    Fitting parameters for void descriptor function analysis.
VoidDescriptorProperties
    Properties for calculation of void descriptor function.
VoidDescriptorPoreProperties
    Properties for void descriptor pore analysis.

Examples
--------
To use `CartesianAxis2D` in your code:

    >>> from custom_types_module import CartesianAxis2D
    >>> axis = CartesianAxis2D.X
    >>> print(axis)
    CartesianAxis2D.X
"""

from typing import NamedTuple
from enum import Enum, IntEnum
import numpy as np
import numpy.typing as npt

##### BASE CLASSES #####


class CartesianAxis2D(IntEnum):
    """
    Constant axis ordering in Python for 2D arrays (Y, X) ordering.

    The 2D numpy arrays are ordered Y, X.

    Attributes
    ----------
    X : int
        The X-axis (1).
    Y : int
        The Y-axis (0).

    Examples
    --------
    >>> axis = CartesianAxis2D.X
    >>> print(axis)
    CartesianAxis2D.X
    """

    X = 1
    Y = 0


class CartesianAxis3D(IntEnum):
    """
    Constant axis ordering in Python for 3D arrays (Z, Y, X) ordering.

    The 3D numpy arrays are ordered Z, Y, X.

    Attributes
    ----------
    X : int
        The X-axis (2).
    Y : int
        The Y-axis (1).
    Z : int
        The Z-axis (0).

    Examples
    --------
    >>> axis = CartesianAxis3D.Z
    >>> print(axis)
    CartesianAxis3D.Z
    """

    X = 2
    Y = 1
    Z = 0


class DataVolume(NamedTuple):
    """
    The integer counts of the data that compose a single image or a volume composed of a stack of images.

    For the 2D case, the units are pixels and the z_image_count=1.
    For the 3D case, the units are voxels, and the z_image_count>=2.

    Attributes
    ----------
    x_width : int
        The width of the data volume in the X dimension.
    y_height : int
        The height of the data volume in the Y dimension.
    z_image_count : int
        The number of images in the Z dimension.

    Examples
    --------
    >>> volume = DataVolume(x_width=8, y_height=3, z_image_count=2)
    >>> print(volume)
    DataVolume(x_width=8, y_height=3, z_image_count=2)
    """

    x_width: int
    y_height: int
    z_image_count: int


class Units(Enum):
    """
    Unit enum for metadata.

    Attributes
    ----------
    NANOMETER : str
        Nanometer unit.
    MICRON : str
        Micron unit.
    MILLIMETER : str
        Millimeter unit.
    METER : str
        Meter unit.
    VOXEL : str
        Voxel unit.
    PIXEL : str
        Pixel unit.
    RADIAN : str
        Radian unit.
    DEGREE : str
        Degree unit.

    Examples
    --------
    >>> unit = Units.MICRON
    >>> print(unit)
    <Units.MICRON: 'micron'>
    """

    NANOMETER: str = "nanometer"
    MICRON: str = "micron"
    MILLIMETER: str = "millimeter"
    METER: str = "meter"
    VOXEL: str = "voxel"
    PIXEL: str = "pixel"
    RADIAN: str = "radian"
    DEGREE: str = "degree"


class UnitVector(NamedTuple):
    """
    Normalized unit vector.

    Attributes
    ----------
    u : float
        The U component of the unit vector.
    v : float
        The V component of the unit vector.
    w : float
        The W component of the unit vector.

    Examples
    --------
    >>> vector = UnitVector(u=1.0, v=0.0, w=0.0)
    >>> print(vector)
    UnitVector(u=1.0, v=0.0, w=0.0)
    """

    u: float
    v: float
    w: float


##### DERIVED CLASSES #####


class Angle(NamedTuple):
    """
    Angle class with associated units.

    Attributes
    ----------
    value : float
        The value of the angle.
    unit : Units
        The unit of the angle.

    Examples
    --------
    >>> angle = Angle(value=90.0, unit=Units.DEGREE)
    >>> print(angle)
    Angle(value=90.0, unit=<Units.DEGREE: 'degree'>)
    """

    value: float
    unit: Units


class Area(NamedTuple):
    """
    Area class with associated units.

    Attributes
    ----------
    value : float
        The value of the area.
    unit_squared : Units
        The unit of the area squared.

    Examples
    --------
    >>> area = Area(value=100.0, unit_squared=Units.MICRON)
    >>> print(area)
    Area(value=100.0, unit_squared=<Units.MICRON: 'micron'>)
    """

    value: float
    unit_squared: Units


class InstanceLabels(NamedTuple):
    """
    Collection of instance labels.

    Attributes
    ----------
    data : np.ndarray[int]
        The array of instance labels.

    Examples
    --------
    >>> labels = InstanceLabels(data=np.array([1, 2, 3]))
    >>> print(labels)
    InstanceLabels(data=array([1, 2, 3]))
    """

    data: np.ndarray[int]


class Length(NamedTuple):
    """
    Length class with associated units.

    Attributes
    ----------
    value : float
        The value of the length.
    unit : Units
        The unit of the length.

    Examples
    --------
    >>> length = Length(value=10.0, unit=Units.MICRON)
    >>> print(length)
    Length(value=10.0, unit=<Units.MICRON: 'micron'>)
    """

    value: float
    unit: Units


class NthNearestNeighbors(NamedTuple):
    """
    Nth nearest neighbors for instance analysis.

    1st nearest neighbors are the features themselves.
    Distances are given in the same length units as centroids.

    Attributes
    ----------
    nth_nearest : int
        The nth nearest neighbor.
    distances : list[Length]
        The list of distances to the nth nearest neighbors.
    instance_id : np.ndarray[int]
        The array of instance IDs.

    Examples
    --------
    >>> neighbors = NthNearestNeighbors(nth_nearest=2, distances=[Length(value=5.0, unit=Units.MICRON)], instance_id=np.array([1, 2, 3]))
    >>> print(neighbors)
    NthNearestNeighbors(nth_nearest=2, distances=[Length(value=5.0, unit=<Units.MICRON: 'micron'>)], instance_id=array([1, 2, 3]))
    """

    nth_nearest: int
    distances: list[Length]
    instance_id: np.ndarray[int]


class NVoxel(NamedTuple):
    """
    Number of voxels, typically for an instance object.

    Attributes
    ----------
    value : int
        The number of voxels.
    unit : Units
        The unit of the voxels (default is Units.VOXEL).

    Examples
    --------
    >>> n_voxel = NVoxel(value=100)
    >>> print(n_voxel)
    NVoxel(value=100, unit=<Units.VOXEL: 'voxel'>)
    """

    value: int
    unit: Units = Units.VOXEL


class Volume(NamedTuple):
    """
    Volume class with associated units.

    Attributes
    ----------
    value : float
        The value of the volume.
    unit_cubed : Units
        The unit of the volume cubed.

    Examples
    --------
    >>> volume = Volume(value=1000.0, unit_cubed=Units.MICRON)
    >>> print(volume)
    Volume(value=1000.0, unit_cubed=<Units.MICRON: 'micron'>)
    """

    value: float
    unit_cubed: Units


##### DOUBLY DERIVED CLASSES #####


class Centroid(NamedTuple):
    """
    The centroid location for the X, Y, and Z axes, in length units.

    Attributes
    ----------
    cx : Length
        The centroid location in the X axis.
    cy : Length
        The centroid location in the Y axis.
    cz : Length
        The centroid location in the Z axis.

    Examples
    --------
    >>> centroid = Centroid(cx=Length(value=1.0, unit=Units.MICRON), cy=Length(value=2.0, unit=Units.MICRON), cz=Length(value=3.0, unit=Units.MICRON))
    >>> print(centroid)
    Centroid(cx=Length(value=1.0, unit=<Units.MICRON: 'micron'>), cy=Length(value=2.0, unit=<Units.MICRON: 'micron'>), cz=Length(value=3.0, unit=<Units.MICRON: 'micron'>))
    """

    cx: Length
    cy: Length
    cz: Length


class Centroids(NamedTuple):
    """
    Centroids for an instance image stack.

    Attributes
    ----------
    data : list[Centroid]
        The list of centroids.

    Examples
    --------
    >>> centroids = Centroids(data=[Centroid(cx=Length(value=1.0, unit=Units.MICRON), cy=Length(value=2.0, unit=Units.MICRON), cz=Length(value=3.0, unit=Units.MICRON))])
    >>> print(centroids)
    Centroids(data=[Centroid(cx=Length(value=1.0, unit=<Units.MICRON: 'micron'>), cy=Length(value=2.0, unit=<Units.MICRON: 'micron'>), cz=Length(value=3.0, unit=<Units.MICRON: 'micron'>))])
    """

    data: list[Centroid]


class EllipsoidAxis(NamedTuple):
    """
    Properties of best fit ellipsoid semi-axis (half of length of principal axis).

    Attributes
    ----------
    length : Length
        The length of the semi-axis.
    orientation : UnitVector
        The orientation of the semi-axis.

    Examples
    --------
    >>> axis = EllipsoidAxis(length=Length(value=5.0, unit=Units.MICRON), orientation=UnitVector(u=1.0, v=0.0, w=0.0))
    >>> print(axis)
    EllipsoidAxis(length=Length(value=5.0, unit=<Units.MICRON: 'micron'>), orientation=UnitVector(u=1.0, v=0.0, w=0.0))
    """

    length: Length
    orientation: UnitVector


class EquivalentSphereDiameters(NamedTuple):
    """
    Collection of spherical diameters.

    Attributes
    ----------
    data : list[Length]
        The list of spherical diameters.

    Examples
    --------
    >>> diameters = EquivalentSphereDiameters(data=[Length(value=10.0, unit=Units.MICRON)])
    >>> print(diameters)
    EquivalentSphereDiameters(data=[Length(value=10.0, unit=<Units.MICRON: 'micron'>)])
    """

    data: list[Length]


class Origin(NamedTuple):
    """
    The origin in spatial coordinates.

    Attributes
    ----------
    x0 : Length
        The origin in the X axis.
    y0 : Length
        The origin in the Y axis.
    z0 : Length
        The origin in the Z axis.

    Examples
    --------
    >>> origin = Origin(x0=Length(value=0.0, unit=Units.MICRON), y0=Length(value=0.0, unit=Units.MICRON), z0=Length(value=0.0, unit=Units.MICRON))
    >>> print(origin)
    Origin(x0=Length(value=0.0, unit=<Units.MICRON: 'micron'>), y0=Length(value=0.0, unit=<Units.MICRON: 'micron'>), z0=Length(value=0.0, unit=<Units.MICRON: 'micron'>))
    """

    x0: Length
    y0: Length
    z0: Length


class Resolution(NamedTuple):
    """
    The voxel resolution for the X, Y, and Z dimension in length units.

    Attributes
    ----------
    dx : Length
        The resolution in the X dimension.
    dy : Length
        The resolution in the Y dimension.
    dz : Length
        The resolution in the Z dimension.

    Examples
    --------
    >>> resolution = Resolution(dx=Length(value=1.0, unit=Units.MICRON), dy=Length(value=1.0, unit=Units.MICRON), dz=Length(value=1.0, unit=Units.MICRON))
    >>> print(resolution)
    Resolution(dx=Length(value=1.0, unit=<Units.MICRON: 'micron'>), dy=Length(value=1.0, unit=<Units.MICRON: 'micron'>), dz=Length(value=1.0, unit=<Units.MICRON: 'micron'>))
    """

    dx: Length
    dy: Length
    dz: Length


##### TRIPLY DERIVED CLASSES #####


class BestFitEllipsoid(NamedTuple):
    """
    Three semi-axes of best-fit ellipsoid in descending order of size.

    Attributes
    ----------
    a : EllipsoidAxis
        The largest semi-axis.
    b : EllipsoidAxis
        The middle semi-axis.
    c : EllipsoidAxis
        The smallest semi-axis.
    """

    a: EllipsoidAxis
    b: EllipsoidAxis
    c: EllipsoidAxis


class BestFitEllipsoids(NamedTuple):
    """
    Collection of best-fit ellipsoids.

    Attributes
    ----------
    data : list[BestFitEllipsoid]
        The list of best-fit ellipsoids.
    """

    data: list[BestFitEllipsoid]


class EllipsoidSurfaceAreas(NamedTuple):
    """
    Collection of ellipsoid surface area objects.

    Attributes
    ----------
    data : list[Area]
        The list of ellipsoid surface areas.
    """

    data: list[Area]


class EllipsoidVolumes(NamedTuple):
    """
    Collection of ellipsoid volumes.

    Attributes
    ----------
    data : list[Volume]
        The list of ellipsoid volumes.
    """

    data: list[Volume]


class MetaData(NamedTuple):
    """
    The data description of the scan.

    Attributes
    ----------
    data_volume : DataVolume
        The dimensions of the volume.
    resolution : Resolution
        The voxel resolution for the X, Y, and Z dimensions.
    pixel_units : str
        The units of the pixels (e.g., micrometer).
    origin : Origin
        The origin in spatial coordinates.
    """

    data_volume: DataVolume
    resolution: Resolution
    pixel_units: str
    origin: Origin


class ImageStack(NamedTuple):
    """
    The generic image stack.

    Attributes
    ----------
    name : str
        The name of the image stack, used to create the path in the HDF5 file.
    metadata : MetaData
        The metadata associated with the image stack.
    data : npt.NDArray
        The image data.
    """

    name: str
    metadata: MetaData
    data: npt.NDArray


class InstanceImageStack(NamedTuple):
    """
    The image stack where each pixel in the image is categorized (labeled) into a specific instance of a class or object.

    Attributes
    ----------
    name : str
        The name of the image stack, used to create the path in the HDF5 file.
    metadata : MetaData
        The metadata associated with the image stack.
    data : npt.NDArray[np.int_]
        The labeled image data.
    nlabels : int
        The number of labels.
    min_feature_size : int
        The minimum feature size.
    """

    name: str
    metadata: MetaData
    data: npt.NDArray[np.int_]
    nlabels: int
    min_feature_size: int


class InstanceIndices(NamedTuple):
    """
    Labels and index locations for each instance in an InstanceImageStack.

    Attributes
    ----------
    source_name : str
        The name of the source InstanceImageStack.
    labels : InstanceLabels
        The instance labels.
    indices : list[np.ndarray[int]]
        The list of index locations for each instance.
    """

    source_name: str
    labels: InstanceLabels
    indices: list[np.ndarray[int]]


class SemanticImageStack(NamedTuple):
    """
    The image stack where each pixel in the image is categorized (labeled) into a class or object.

    Attributes
    ----------
    name : str
        The name of the image stack, used to create the path in the HDF5 file.
    metadata : MetaData
        The metadata associated with the image stack.
    data : npt.NDArray[np.int_]
        The segmented image data, of type int.
    """

    name: str
    metadata: MetaData
    data: npt.NDArray[np.int_]


class VoidDescriptorEllipsoidAngles(NamedTuple):
    """
    Angles for ellipsoid properties to be used in void descriptor function analysis.

    Attributes
    ----------
    theta_xy : Angle
        The angle in the XY plane.
    phi_x : Angle
        The angle in the X direction.
    phi_y : Angle
        The angle in the Y direction.
    """

    theta_xy: Angle
    phi_x: Angle
    phi_y: Angle


class InstanceProperties(NamedTuple):
    """
    Properties for labeled features in an instance image stack.

    Attributes
    ----------
    source_name : str
        The name of the source InstanceImageStack.
    labels : InstanceLabels
        The instance labels.
    n_voxels : list[NVoxel]
        The number of voxels for each instance.
    equivalent_sphere_diameters : list[Length]
        The equivalent sphere diameters for each instance.
    centroids : Centroids
        The centroids for each instance.
    ellipsoids : BestFitEllipsoids
        The best-fit ellipsoids for each instance.
    surface_areas : EllipsoidSurfaceAreas
        The surface areas of the ellipsoids for each instance.
    volumes : EllipsoidVolumes
        The volumes of the ellipsoids for each instance.
    """

    source_name: str
    labels: InstanceLabels
    n_voxels: list[NVoxel]
    equivalent_sphere_diameters: list[Length]
    centroids: Centroids
    ellipsoids: BestFitEllipsoids
    surface_areas: EllipsoidSurfaceAreas
    volumes: EllipsoidVolumes


class VoidDescriptorFittingParameters(NamedTuple):
    """
    Fitting parameters for void descriptor function analysis.

    Attributes
    ----------
    alpha : float
        The alpha fitting parameter.
    rho : float
        The rho fitting parameter.
    gamma : float
        The gamma fitting parameter.
    zeta : float
        The zeta fitting parameter.
    """

    alpha: float
    rho: float
    gamma: float
    zeta: float


class VoidDescriptorProperties(NamedTuple):
    """
    Properties for calculation of void descriptor function.

    Attributes
    ----------
    fitting_parameters : VoidDescriptorFittingParameters
        The fitting parameters for the void descriptor function.
    metadata : MetaData
        The metadata associated with the void descriptor properties.
    specimen_volume : Volume
        The volume of the specimen.
    loading_direction : CartesianAxis3D
        The loading direction.
    max_radial_distance : Length
        The maximum radial distance from the center to the free surface.
    gauge_length : Length
        The gauge length.
    """

    fitting_parameters: VoidDescriptorFittingParameters
    metadata: MetaData
    specimen_volume: Volume
    loading_direction: CartesianAxis3D
    max_radial_distance: Length
    gauge_length: Length


class VoidDescriptorPoreProperties(NamedTuple):
    """
    Properties for void descriptor pore analysis.

    Attributes
    ----------
    centroid_px : Centroid
        The centroid of the pore.
    volume : Volume
        The volume of the pore.
    surface_area : Area
        The surface area of the pore.
    r_a : Length
        The semi-major axis length of the pore.
    r_b : Length
        The semi-intermediate axis length of the pore.
    r_c : Length
        The semi-minor axis length of the pore.
    vdf_angles : VoidDescriptorEllipsoidAngles
        The angles for the ellipsoid properties.
    """

    centroid_px: Centroid
    volume: Volume
    surface_area: Area
    r_a: Length
    r_b: Length
    r_c: Length
    vdf_angles: VoidDescriptorEllipsoidAngles
