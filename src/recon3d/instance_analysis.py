"""
Basic workflow for processing grayscale image stacks to labeled images for feature property measurement.

This module supports multiple phases (>2 phases) and includes simple segmentation for demonstration purposes. 
It writes outputs as image directories, but can be integrated with HDF for internal dataset storage.
"""

# TODO Integrate with HDF writing
# TODO Integrate with feature measurements
# TODO Use constants for Z axis in other module
# TODO Simplify/Impose structure for image directory (would be dataset directory in HDF)
# to remove magic numbers and indexing of phase_paths and label_paths
# TODO Devise better solution for multiphase support without if/else blocks


import argparse
import glob
from typing import Final, NamedTuple, Tuple, Type, List
import math

import h5py
import numpy as np
import numpy.typing as npt
from pathlib import Path
from PIL import Image
import skimage.io as skio
from sklearn.neighbors import NearestNeighbors as sklearnNN

from scipy import ndimage
from skimage.measure import label, regionprops
import cc3d
import fastremap

import recon3d.constants as cs
import recon3d.types
import recon3d.utility as ut
import recon3d.hdf_io as hio
from recon3d.types import *


# # ----------------------------------
# # USER INPUT DATA HARD CODES - begin
# # ----------------------------------
# # NEED TO PUT INTO A .yml file instead of hard codes
# # input_dir: str = "/Users/apolon/Desktop/test_recon/input/"
# # save_dir: str = "/Users/apolon/Desktop/test_recon/output/"

# input_dir: Final[str] = "/Users/chovey/temp/test_recon/input/"
# save_dir: Final[str] = "/Users/chovey/temp/test_recon/output/"

# file_type: Final[str] = ".tif"
# gray_path: Final[str] = f"{save_dir}/gray/"
# threshold: Final[int] = 185
# segmented_path: Final[str] = f"{save_dir}/segmented/"
# # first phase is always assumed to be the sample (1 or 255 in segmented images)
# phase_paths: np.ndarray = np.array(
#     [f"{save_dir}/sample_mask/", f"{save_dir}/void_mask/"]
# )
# # first label path is always assumed to the first phase after the sample,
# # same ordering as phase_paths
# label_paths: np.ndarray = np.array([f"{save_dir}/void_ids/"])

# # --------------------------------
# # USER INPUT DATA HARD CODES - end
# # --------------------------------


def calc_moment(
    indices: np.ndarray[np.int32],
    p: int,
    q: int,
    r: int,
    centroid_px: Centroid,
) -> float:
    """
    Calculate 3D moment invariants of a feature based on its voxel locations.
    Mask is given in index locations, so math here is all done in pixel units.

    Parameters
    ----------
    indices : np.ndarray[np.int32]
        Array of [Z,Y,X] indices corresponding to the feature's voxel locations.
        The shape of the array should be (n, 3), where n is the number of voxels of the feature.
        Example: np.array([[Z1, Y1, X1], [Z2, Y2, X2], ..., [Zn, Yn, Xn]]).
    p : int
        Order of the moment for the X dimension.
    q : int
        Order of the moment for the Y dimension.
    r : int
        Order of the moment for the Z dimension.
    centroid_px : Centroid
        Centroid location for X, Y, and Z axis, respectively, in pixel units.

    Returns
    -------
    float
        The 3D central moment (translation invariant) in pixel units normalized by object volume.

    Raises
    ------
    ValueError
        If the units of the centroid are not in pixel/voxel units.

    Examples
    --------
    >>> indices = np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2]], dtype=np.int32)
    >>> centroid_px = Centroid(
    ...     cx=Length(1.0, Units.VOXEL),
    ...     cy=Length(1.0, Units.VOXEL),
    ...     cz=Length(1.0, Units.VOXEL)
    ... )
    >>> calc_moment(indices, 2, 2, 2, centroid_px)
    0.66
    """

    unit = centroid_px.cx.unit
    if unit != Units.VOXEL:
        raise ValueError(
            f"Moment invariants are calculated using pixel/voxel units, but {unit} units were given"
        )

    mu = 0
    for i in range(0, indices.shape[0]):
        mu += (
            (indices[i][2] - centroid_px.cx.value)
            ** p  # x-location of voxel w.r.t x-centroid
            * (indices[i][1] - centroid_px.cy.value)
            ** q  # y-location of voxel w.r.t y-centroid
            * (indices[i][0] - centroid_px.cz.value)
            ** r  # z-location of voxel w.r.t z-centroid
        )
    return mu / len(indices)


def center_of_mass(
    indices: np.ndarray[np.int32], resolution: Resolution, origin: Origin
) -> Centroid:
    """
    Calculates the center of mass (centroid) of an instance object.
    Requires that the instance object has uniform density (is homogeneous).

    Parameters
    ----------
    indices : np.ndarray[np.int32]
        Array of [Z,Y,X] indices corresponding to the feature's voxel locations.
        The shape of the array should be (n, 3), where n is the number of voxels of the feature.
        Example: np.array([[Z1, Y1, X1], [Z2, Y2, X2], ..., [Zn, Yn, Xn]]).
    resolution : Resolution
        Resolution of the voxel grid in each dimension (dx, dy, dz).
    origin : Origin
        Origin of the voxel grid in each dimension (x0, y0, z0).

    Returns
    -------
    Centroid
        The centroid location for X, Y, and Z axis, respectively, in microns.

    Examples
    --------
    >>> indices = np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2]], dtype=np.int32)
    >>> resolution = Resolution(dx=Length(1.0, Units.MICRON), dy=Length(1.0, Units.MICRON), dz=Length(1.0, Units.MICRON))
    >>> origin = Origin(x0=Length(0.0, Units.MICRON), y0=Length(0.0, Units.MICRON), z0=Length(0.0, Units.MICRON))
    >>> center_of_mass(indices, resolution, origin)
    Centroid(cx=Length(value=1.0, unit=Units.MICRON), cy=Length(value=1.0, unit=Units.MICRON), cz=Length(value=1.0, unit=Units.MICRON))
    """

    centroid_zyx = np.add(
        np.multiply(
            np.mean(indices, axis=0),
            np.array(
                [
                    resolution.dz.value,
                    resolution.dy.value,
                    resolution.dx.value,
                ]
            ),
        ),
        np.array(
            [
                origin.z0.value,
                origin.y0.value,
                origin.x0.value,
            ]
        ),
    )
    cx = centroid_zyx[2]
    cy = centroid_zyx[1]
    cz = centroid_zyx[0]
    unit = resolution.dx.unit
    return Centroid(
        cx=Length(cx, unit=unit),
        cy=Length(cy, unit=unit),
        cz=Length(cz, unit=unit),
    )


def centroid_pix(
    centroid: Centroid, resolution: Resolution, origin: Origin
) -> Centroid:
    """
    Convert between real space units and pixel units.

    Parameters
    ----------
    centroid : Centroid
        The centroid location in real space units.
    resolution : Resolution
        Resolution of the voxel grid in each dimension (dx, dy, dz).
    origin : Origin
        Origin of the voxel grid in each dimension (x0, y0, z0).

    Returns
    -------
    Centroid
        The centroid location in pixel units.

    Examples
    --------
    >>> centroid = Centroid(cx=Length(10.0, Units.MICRON), cy=Length(20.0, Units.MICRON), cz=Length(30.0, Units.MICRON))
    >>> resolution = Resolution(dx=Length(1.0, Units.MICRON), dy=Length(1.0, Units.MICRON), dz=Length(1.0, Units.MICRON))
    >>> origin = Origin(x0=Length(0.0, Units.MICRON), y0=Length(0.0, Units.MICRON), z0=Length(0.0, Units.MICRON))
    >>> centroid_pix(centroid, resolution, origin)
    Centroid(cx=Length(value=10.0, unit=Units.VOXEL), cy=Length(value=20.0, unit=Units.VOXEL), cz=Length(value=30.0, unit=Units.VOXEL))
    """

    # convert centroids to normalized pixel units
    cx = (centroid.cx.value - origin.x0.value) / resolution.dx.value
    cy = (centroid.cy.value - origin.y0.value) / resolution.dy.value
    cz = (centroid.cz.value - origin.z0.value) / resolution.dz.value

    unit = Units.VOXEL

    return Centroid(
        cx=Length(cx, unit=unit),
        cy=Length(cy, unit=unit),
        cz=Length(cz, unit=unit),
    )


def ellipsoid_surface_area(ellipsoid: BestFitEllipsoid) -> Area:
    """
    Calculate the surface area of the best fit ellipsoid using the Knud Thomsen formula.
    Applicable to scalene ellipsoids (semi-axes a > b > c), with a relative error of at most 1.061%.

    Parameters
    ----------
    ellipsoid : BestFitEllipsoid
        The best fit ellipsoid with semi-axes lengths a, b, and c.

    Returns
    -------
    Area
        The surface area of the ellipsoid in microns squared.

    Examples
    --------
    >>> ellipsoid = BestFitEllipsoid(
    ...     a=EllipsoidAxis(length=Length(5.0, Units.MICRON), orientation=UnitVector(u=1, v=0, w=0)),
    ...     b=EllipsoidAxis(length=Length(3.0, Units.MICRON), orientation=UnitVector(u=0, v=1, w=0)),
    ...     c=EllipsoidAxis(length=Length(2.0, Units.MICRON), orientation=UnitVector(u=0, v=0, w=1))
    ... )
    >>> ellipsoid_surface_area(ellipsoid)
    Area(value=134.8149867941625, unit_squared=<Units.MICRON: 'micron'>)
    """
    unit = ellipsoid.a.length.unit
    p = cs.Constants().KNUD_THOMSEN_FACTOR
    a, b, c = (
        ellipsoid.a.length.value,
        ellipsoid.b.length.value,
        ellipsoid.c.length.value,
    )

    surface_area = (
        4 * math.pi * ((((a**p * b**p) + (a**p * c**p) + (b**p * c**p)) / 3) ** (1 / p))
    )
    return Area(value=surface_area, unit_squared=unit)


def ellipsoid_volume(ellipsoid: BestFitEllipsoid) -> Volume:
    """
    Calculate the volume of the best fit ellipsoid using semi-axes lengths.

    Parameters
    ----------
    ellipsoid : BestFitEllipsoid
        The best fit ellipsoid with semi-axes lengths a, b, and c.

    Returns
    -------
    Volume
        The volume of the ellipsoid in microns cubed.

    Examples
    --------
    >>> ellipsoid = BestFitEllipsoid(
    ...     a=EllipsoidAxis(length=Length(5.0, Units.MICRON), orientation=UnitVector(u=1, v=0, w=0)),
    ...     b=EllipsoidAxis(length=Length(3.0, Units.MICRON), orientation=UnitVector(u=0, v=1, w=0)),
    ...     c=EllipsoidAxis(length=Length(2.0, Units.MICRON), orientation=UnitVector(u=0, v=0, w=1))
    ... )
    >>> ellipsoid_volume(ellipsoid)
    Volume(value=125.66370614359171, unit_cubed=<Units.MICRON: 'micron'>)
    """

    unit = ellipsoid.a.length.unit
    a, b, c = (
        ellipsoid.a.length.value,
        ellipsoid.b.length.value,
        ellipsoid.c.length.value,
    )

    volume = 4 / 3 * math.pi * a * b * c

    return Volume(value=volume, unit_cubed=unit)


def equivalent_spherical_diameter(
    n_voxels: int,
    resolution: Resolution,
) -> Length:
    """
    Calculate the equivalent spherical diameter of a feature.

    Parameters
    ----------
    n_voxels : int
        The number of voxels in the feature.
    resolution : Resolution
        The resolution of the voxel grid in each dimension (dx, dy, dz).

    Returns
    -------
    Length
        The equivalent spherical diameter of the feature in the same units as the resolution.

    Examples
    --------
    >>> resolution = Resolution(
    ...     dx=Length(1.0, Units.MICRON),
    ...     dy=Length(1.0, Units.MICRON),
    ...     dz=Length(1.0, Units.MICRON)
    ... )
    >>> equivalent_spherical_diameter(1000, resolution)
    Length(value=12.407009817988, unit=<Units.MICRON: 'micron'>)
    """

    value = 2 * np.power(
        (
            (3 / (4 * np.pi))
            * (
                n_voxels
                * resolution.dx.value
                * resolution.dy.value
                * resolution.dz.value
            )
        ),
        (1 / 3),
    )
    return Length(value=value, unit=resolution.dx.unit)


def fit_ellipsoid(
    indices: np.ndarray,
    centroid: Centroid,
    resolution: Resolution,
    origin: Origin,
) -> BestFitEllipsoid:
    """
    Fit a 3D ellipsoid to a set of voxel indices.

    Parameters
    ----------
    indices : np.ndarray
        Array of [Z,Y,X] indices corresponding to the feature's voxel locations.
        The shape of the array should be (n, 3), where n is the number of voxels of the feature.
    centroid : Centroid
        Centroid location for X, Y, and Z axis, respectively, in spatial coordinates (e.g., microns).
    resolution : Resolution
        Voxel resolution for X, Y, and Z dimensions, respectively, in length units (e.g., microns).
    origin : Origin
        Optional origin in spatial coordinates if the volume has been previously transformed.

    Returns
    -------
    BestFitEllipsoid
        The best fit ellipsoid with semi-axes lengths and orientations.

    Examples
    --------
    >>> indices = np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2]], dtype=np.int32)
    >>> centroid = Centroid(cx=Length(1.0, Units.MICRON), cy=Length(1.0, Units.MICRON), cz=Length(1.0, Units.MICRON))
    >>> resolution = Resolution(dx=Length(1.0, Units.MICRON), dy=Length(1.0, Units.MICRON), dz=Length(1.0, Units.MICRON))
    >>> origin = Origin(x0=Length(0.0, Units.MICRON), y0=Length(0.0, Units.MICRON), z0=Length(0.0, Units.MICRON))
    >>> fit_ellipsoid(indices, centroid, resolution, origin)
    BestFitEllipsoid(a=EllipsoidAxis(length=Length(value=80.146386085398, unit=<Units.MICRON: 'micron'>), orientation=UnitVector(u=-0.25175905655605574, v=-0.5305735687527386, w=-0.8093880809494217)), b=EllipsoidAxis(length=Length(value=1.2477168951003976, unit=<Units.MICRON: 'micron'>), orientation=UnitVector(u=-0.8774683531474325, v=-0.227651095921704, w=0.4221661613042665)), c=EllipsoidAxis(length=Length(value=2.380014006104306e-07, unit=<Units.MICRON: 'micron'>), orientation=UnitVector(u=0.4082482904637247, v=-0.816496580927762, w=0.4082482904639296)))
    """

    if not resolution.dx == resolution.dy == resolution.dz:
        raise NotImplementedError(
            "Only isotropic voxels are currently supported for axis length calculation"
        )

    unit = resolution.dx.unit

    # convert centroids to normalized pixel units
    cen_px = centroid_pix(centroid, resolution, origin)

    # Get moments
    s_xx = calc_moment(indices, 2, 0, 0, cen_px)
    s_yy = calc_moment(indices, 0, 2, 0, cen_px)
    s_zz = calc_moment(indices, 0, 0, 2, cen_px)
    s_xy = calc_moment(indices, 1, 1, 0, cen_px)
    s_xz = calc_moment(indices, 1, 0, 1, cen_px)
    s_yz = calc_moment(indices, 0, 1, 1, cen_px)

    # Create matrix of 3D moments
    s_ij = np.array(
        [
            [s_xx, s_xy, s_xz],
            [s_xy, s_yy, s_yz],
            [s_xz, s_yz, s_zz],
        ]
    )
    eig_val, eig_vec = np.linalg.eig(s_ij)

    semi_ax_lengths_px = np.sqrt(5 * eig_val)  # still in pixel units
    sort = np.argsort(semi_ax_lengths_px)[::-1]  # sort descending axis lengths

    # TODO: math to account for change in length as function
    #       as function of axis orientation in anisotropic voxels
    #       remove NotImplementedError above
    semi_ax_lengths = np.array(
        [
            semi_ax_lengths_px[sort[0]] * resolution.dx.value,
            semi_ax_lengths_px[sort[1]] * resolution.dx.value,
            semi_ax_lengths_px[sort[2]] * resolution.dx.value,
        ]
    )
    ## some axes lengths are NAN, should overwrite to 0.0
    semi_ax_lengths = np.nan_to_num(semi_ax_lengths, nan=0.0)

    ## from docs on np.linalg.eig: The normalized (unit “length”) eigenvectors,
    ## such that the column eigenvectors[:,i] is the eigenvector corresponding to the eigenvalue eigenvalues[i].
    ax_vectors = np.array(
        [
            eig_vec[:, sort[0]],
            eig_vec[:, sort[1]],
            eig_vec[:, sort[2]],
        ]
    )
    ax_vectors = np.reshape(ax_vectors, (9))

    a_vector = UnitVector(
        u=ax_vectors[0],
        v=ax_vectors[1],
        w=ax_vectors[2],
    )
    a = EllipsoidAxis(
        length=Length(semi_ax_lengths[0], unit=unit),
        orientation=a_vector,
    )
    b_vector = UnitVector(
        u=ax_vectors[3],
        v=ax_vectors[4],
        w=ax_vectors[5],
    )
    b = EllipsoidAxis(
        length=Length(semi_ax_lengths[1], unit=unit),
        orientation=b_vector,
    )
    c_vector = UnitVector(
        u=ax_vectors[6],
        v=ax_vectors[7],
        w=ax_vectors[8],
    )
    c = EllipsoidAxis(
        length=Length(semi_ax_lengths[2], unit=unit),
        orientation=c_vector,
    )

    ellipsoid = BestFitEllipsoid(a=a, b=b, c=c)
    return ellipsoid


def instance_indices(
    instance_stack: InstanceImageStack,
) -> InstanceIndices:
    """
    Returns an array of instances and a list of N-dimensional arrays of indices
    of the unique values in the instance stack.

    Parameters
    ----------
    instance_stack : InstanceImageStack
        Array with arbitrary dimensions.

    Returns
    -------
    InstanceIndices
        - 1D-array of sorted unique values.
        - List of arrays. Each array contains the indices where a given value in x is found.

    Examples
    --------
    >>> instance_stack = InstanceImageStack(
    ...     name="example_stack",
    ...     metadata=MetaData(
    ...         data_volume=DataVolume(x_width=2, y_height=2, z_image_count=1),
    ...         resolution=Resolution(
    ...             dx=Length(1.0, Units.MICRON),
    ...             dy=Length(1.0, Units.MICRON),
    ...             dz=Length(1.0, Units.MICRON)
    ...         ),
    ...         pixel_units=Units.MICRON,
    ...         origin=Origin(
    ...             x0=Length(0.0, Units.MICRON),
    ...             y0=Length(0.0, Units.MICRON),
    ...             z0=Length(0.0, Units.MICRON)
    ...         )
    ...     ),
    ...     data=np.array([[1, 2, 2], [3, 1, 1]]),
    ...     nlabels=3,
    ...     min_feature_size=1
    ... )
    >>> instance_indices(instance_stack)
    InstanceIndices(
        source_name='example_stack',
        indices=[array([[0, 0], [1, 1], [1, 2]]), array([[0, 1], [0, 2]]), array([[1, 0]])],
        labels=InstanceLabels(data=np.array([1, 2, 3]))
    )
    """

    x = instance_stack.data
    x_flat = x.ravel()
    ix_flat = np.argsort(x_flat)
    u, ix_u = fastremap.unique(x_flat[ix_flat], return_index=True)
    ix_ndim = np.unravel_index(ix_flat, x.shape)
    ix_ndim = np.c_[ix_ndim] if x.ndim > 1 else ix_flat
    indices = np.split(ix_ndim, ix_u[1:])
    labels = u
    instance_ids = InstanceIndices(
        source_name=instance_stack.name,
        indices=indices,
        labels=InstanceLabels(data=labels),
    )
    return instance_ids


def minimum_size_filter(
    initial_stack: InstanceImageStack,
    initial_indices: InstanceIndices,
) -> Tuple[InstanceImageStack, InstanceIndices]:
    """
    Removes features below minimum size from InstanceIndices and corresponding
    InstanceImageStack. Should be done before shape metrics are determined.

    Parameters
    ----------
    initial_stack : InstanceImageStack
        The initial instance image stack.
    initial_indices : InstanceIndices
        The initial instance indices.

    Returns
    -------
    Tuple[InstanceImageStack, InstanceIndices]
        The filtered instance image stack and the filtered instance indices.

    Raises
    ------
    ValueError
        If the names of the initial stack and indices do not match.

    Examples
    --------
    >>> initial_stack = InstanceImageStack(
    ...     name="example_stack",
    ...     metadata=MetaData(
    ...         data_volume=DataVolume(x_width=4, y_height=2, z_image_count=1),
    ...         resolution=Resolution(
    ...             dx=Length(1.0, Units.MICRON),
    ...             dy=Length(1.0, Units.MICRON),
    ...             dz=Length(1.0, Units.MICRON)
    ...         ),
    ...         pixel_units=Units.MICRON,
    ...         origin=Origin(
    ...             x0=Length(0.0, Units.MICRON),
    ...             y0=Length(0.0, Units.MICRON),
    ...             z0=Length(0.0, Units.MICRON)
    ...         )
    ...     ),
    ...     data=np.array([[1, 2, 2, 1], [3, 1, 1, 1]]),
    ...     nlabels=3,
    ...     min_feature_size=2
    ... )
    >>> initial_indices = instance_indices(initial_stack)
    >>> minimum_size_filter(initial_stack, initial_indices)
    2 connected components remaining after minimum size filter of 2 voxels
    (InstanceImageStack(name='example_stack', metadata=MetaData(data_volume=DataVolume(x_width=4, y_height=2, z_image_count=1), resolution=Resolution(dx=Length(value=1.0, unit=<Units.MICRON: 'micron'>), dy=Length(value=1.0, unit=<Units.MICRON: 'micron'>), dz=Length(value=1.0, unit=<Units.MICRON: 'micron'>)), pixel_units=<Units.MICRON: 'micron'>, origin=Origin(x0=Length(value=0.0, unit=<Units.MICRON: 'micron'>), y0=Length(value=0.0, unit=<Units.MICRON: 'micron'>), z0=Length(value=0.0, unit=<Units.MICRON: 'micron'>))), data=array([[1, 0, 0, 1],[2, 1, 1, 1]], dtype=int8), nlabels=2, min_feature_size=2), InstanceIndices(source_name='example_stack', labels=InstanceLabels(data=array([0, 1, 2], dtype=int8)), indices=[array([[0, 1],[0, 2]]), array([[0, 0],[0, 3],[1, 1],[1, 2],[1, 3]]), array([[1, 0]])]))
    """

    if not initial_stack.name == initial_indices.source_name:
        raise ValueError(
            f"""Initial InstanceImageStack name ({initial_stack.name}) and 
            InstanceIndices name ({initial_indices.source_name}) do not match.
            Function inputs must be from same dataset."""
        )

    n_voxels = num_voxels(initial_indices)
    n_voxels_vals = np.asarray([i.value for i in n_voxels])
    invalid_instance_ids = np.nonzero(n_voxels_vals < initial_stack.min_feature_size)[0]
    masked_instance_ids = fastremap.mask(initial_stack.data, invalid_instance_ids)
    remapped_instance_ids, _ = fastremap.renumber(
        masked_instance_ids,
        preserve_zero=True,
    )
    nlabels = np.max(np.unique(remapped_instance_ids))

    filtered_stack = InstanceImageStack(
        name=initial_stack.name,
        metadata=initial_stack.metadata,
        data=remapped_instance_ids,
        nlabels=nlabels,
        min_feature_size=initial_stack.min_feature_size,
    )

    filtered_indices = instance_indices(instance_stack=filtered_stack)

    print(
        f"\t\t{nlabels} connected components remaining after minimum size filter of {initial_stack.min_feature_size} voxels"
    )
    return filtered_stack, filtered_indices


# TODO
def map_features_to_voxels():
    """To come."""


def nearest_neighbor_distance(
    centroids: Centroids,
) -> NthNearestNeighbors:
    """
    Calculate the nearest neighbor distance for a set of centroids.
    With nth_nearest_neighbor=1 returning the features themselves
    (distances=0, neighbor is themselves).

    Parameters
    ----------
    centroids : Centroids
        The centroids for which to calculate the nearest neighbor distance.

    Returns
    -------
    NthNearestNeighbors
        The nth nearest neighbors, distances, and instance IDs.

    Examples
    --------
    >>> centroids = Centroids(data=[
    ...     Centroid(cx=Length(0.0, Units.MICRON), cy=Length(0.0, Units.MICRON), cz=Length(0.0, Units.MICRON)),
    ...     Centroid(cx=Length(1.0, Units.MICRON), cy=Length(1.0, Units.MICRON), cz=Length(1.0, Units.MICRON)),
    ...     Centroid(cx=Length(2.0, Units.MICRON), cy=Length(2.0, Units.MICRON), cz=Length(2.0, Units.MICRON))
    ... ])
    >>> nearest_neighbor_distance(centroids)
    NthNearestNeighbors(nth_nearest=2, distances=[Length(value=1.7320508075688772, unit=<Units.MICRON: 'micron'>), Length(value=1.7320508075688772, unit=<Units.MICRON: 'micron'>), Length(value=1.7320508075688772, unit=<Units.MICRON: 'micron'>)], instance_id=array([0, 0, 1]))
    """

    neighbor = sklearnNN(n_neighbors=2)
    centroid_array = ut.centroids_to_ndarray(centroids=centroids)
    neighbor.fit(centroid_array)
    neighbor_dist, neighbor_index = neighbor.kneighbors(
        centroid_array,
        return_distance=True,
    )
    distances = neighbor_dist[:, 1]
    instance_id = neighbor_index[:, 1]
    distances[0], instance_id[0] = 0, 0  # set feature 0 to 0 distance and neighbor 0

    unit = centroids.data[0].cx.unit
    distances_list = [Length(value=i, unit=unit) for i in distances.tolist()]
    nth_nearest = NthNearestNeighbors(
        nth_nearest=2,
        distances=distances_list,
        instance_id=instance_id,
    )
    return nth_nearest


def num_voxels(indices: InstanceIndices) -> list[NVoxel]:
    """
    Returns the number of voxels for each instance ID from an InstanceImageStack.

    Parameters
    ----------
    indices : InstanceIndices
        The instance indices containing the indices of each instance.

    Returns
    -------
    List[NVoxel]
        A list of NVoxel objects, each representing the number of voxels for an instance.

    Examples
    --------
    >>> indices = InstanceIndices(
    ...     source_name="example_stack",
    ...     indices=[np.array([[0, 0], [1, 1]]), np.array([[0, 1], [1, 2]])],
    ...     labels=InstanceLabels(data=np.array([1, 2]))
    ... )
    >>> num_voxels(indices)
    [NVoxel(value=2, unit=<Units.VOXEL: 'voxel'>), NVoxel(value=2, unit=<Units.VOXEL: 'voxel'>)]
    """

    indices = indices.indices
    data = [NVoxel(value=len(i)) for i in indices]
    # n_voxels = np.asarray(data, dtype=NVoxel)
    n_voxels = data
    return n_voxels


def process_image_stack(yml_input_file: Path) -> SemanticImageStack:
    """
    Given a yml input file, converts it to a SemanticImageStack.

    Parameters
    ----------
    yml_input_file : Path
        The fully pathed yml file that specifies settings.

    Returns
    -------
    SemanticImageStack
        The SemanticImageStack object.

    Examples
    --------
    >>> yml_input_file = Path("path/to/input.yml")
    >>> semantic_image_stack = process_image_stack(yml_input_file)
    """

    print(f"Processing specification file: {yml_input_file}")
    db = ut.yaml_to_dict(yml_input_file)

    file_type = db["semantic_images_type"]
    valid_extensions = (".tif", ".tiff")

    if file_type not in valid_extensions:
        raise ValueError("Error: 'image_type' must be '.tif' or '.tiff'.")

    path_input = Path(db["semantic_images_dir"]).expanduser()
    if not path_input.is_dir():
        raise ValueError(f"Error, 'image_dir' of {path_input} not found.")
    print(f"Input path: {path_input}")

    path_output = Path(db["out_dir"]).expanduser()
    if not path_output.is_dir():
        Path(path_output).mkdir(
            parents=True, exist_ok=True
        )  # makes directory if it doesn't exist
    print(f"Output path: {path_output}")

    # semantic_seg_dict = db["class_labels"]

    # # load images
    # imageStack = skio.imread_collection(tiff_files, conserve_memory=True)
    data = ut.read_images(path_input, file_type)
    semantic_seg_stack_name = db["semantic_imagestack_name"]

    # create meta data
    (nz, ny, nx) = data.shape
    data_volume = DataVolume(z_image_count=nz, y_height=ny, x_width=nx)
    print(
        f"image stack, {semantic_seg_stack_name}, has dimensions (num_images, row, col): {data.shape}"
    )

    pixel_units = db["pixel_units"]
    valid_units = set(item.value for item in Units)
    if pixel_units not in valid_units:
        raise ValueError(
            f"Error, '{pixel_units}' is not a valid unit, accepted units are: {valid_units}"
        )
    pixel_units = Units(db["pixel_units"])

    dx, dy, dz = (
        db["voxel_size"]["dx"],
        db["voxel_size"]["dy"],
        db["voxel_size"]["dz"],
    )
    resolution = Resolution(
        dx=Length(dx, unit=pixel_units),
        dy=Length(dy, unit=pixel_units),
        dz=Length(dz, unit=pixel_units),
    )

    x0, y0, z0 = (
        db["origin"]["x0"],
        db["origin"]["y0"],
        db["origin"]["z0"],
    )
    origin = Origin(
        x0=Length(x0, unit=pixel_units),
        y0=Length(y0, unit=pixel_units),
        z0=Length(z0, unit=pixel_units),
    )

    meta = MetaData(
        data_volume=data_volume,
        resolution=resolution,
        pixel_units=pixel_units,
        origin=origin,
    )

    # # Optional, grey image stack
    # # TODO: test functionality
    # if db["grey_image_dir"]:
    #     path_grey_image_dir = Path(db["path_grey_image_dir"]).expanduser()
    #     assert path_grey_image_dir.is_dir(), "Error, 'image_dir' not found."
    #     grey_image_name = "Raw greyscale"
    #     grey_data = ut.read_images(path_grey_image_dir, file_type)
    #     assert data.shape == grey_data.shape
    #     semantic_seg_image_stack = ImageStack(
    #         name=grey_image_name, metadata=meta, data=grey_data
    #     )

    # create SemanticImageStack
    semantic_seg_image_stack = SemanticImageStack(
        name=semantic_seg_stack_name, metadata=meta, data=data
    )

    return semantic_seg_image_stack  # , path_file_output)


def save_instance_images(
    input_volumes: list[InstanceImageStack], save_path: Path
) -> bool:
    """
    Save InstanceImageStack objects to image files.

    This function takes a list of InstanceImageStack objects and saves each one as a series of image files.
    The images are saved in the specified directory, with each stack saved in its own subdirectory.

    Parameters
    ----------
    input_volumes : list of InstanceImageStack
        List of InstanceImageStack objects to be saved as images.
    save_path : Path
        The directory where the images will be saved.

    Returns
    -------
    bool
        True if the images were saved successfully, False otherwise.

    Examples
    --------
    >>> input_volumes = [
    ...     InstanceImageStack(
    ...         name="example_stack_1",
    ...         metadata=MetaData(
    ...             data_volume=DataVolume(x_width=256, y_height=256, z_image_count=10),
    ...             resolution=Resolution(
    ...                 dx=Length(1.0, Units.MICRON),
    ...                 dy=Length(1.0, Units.MICRON),
    ...                 dz=Length(1.0, Units.MICRON)
    ...             ),
    ...             pixel_units=Units.MICRON,
    ...             origin=Origin(
    ...                 x0=Length(0.0, Units.MICRON),
    ...                 y0=Length(0.0, Units.MICRON),
    ...                 z0=Length(0.0, Units.MICRON)
    ...             )
    ...         ),
    ...         data=np.random.rand(10, 256, 256),
    ...         nlabels=3,
    ...         min_feature_size=2
    ...     ),
    ...     InstanceImageStack(
    ...         name="example_stack_2",
    ...         metadata=MetaData(
    ...             data_volume=DataVolume(x_width=256, y_height=256, z_image_count=10),
    ...             resolution=Resolution(
    ...                 dx=Length(1.0, Units.MICRON),
    ...                 dy=Length(1.0, Units.MICRON),
    ...                 dz=Length(1.0, Units.MICRON)
    ...             ),
    ...             pixel_units=Units.MICRON,
    ...             origin=Origin(
    ...                 x0=Length(0.0, Units.MICRON),
    ...                 y0=Length(0.0, Units.MICRON),
    ...                 z0=Length(0.0, Units.MICRON)
    ...             )
    ...         ),
    ...         data=np.random.rand(10, 256, 256),
    ...         nlabels=3,
    ...         min_feature_size=2
    ...     )
    ... ]
    >>> save_path = Path("path/to/save")
    >>> save_instance_images(input_volumes, save_path)
    True
    """

    for instance_volumes in input_volumes:
        ut.ndarray_to_img(
            data=instance_volumes.data,
            slice_axis=recon3d.types.CartesianAxis3D.Z,
            parent_dir=save_path,
            folder_name=instance_volumes.name,
            file_type=".tif",
            pad_length=4,
        )


def semantic_to_instance(
    semantic_stack: SemanticImageStack,
    instance_name: str,
    instance_value: int,
    min_feature_size: int,
) -> InstanceImageStack:
    """
    Isolate each unique class in the semantic segmentation (e.g., Cats or Dogs)
    and create the Instance Image Stacks.

    Parameters
    ----------
    semantic_stack : SemanticImageStack
        The semantic image stack containing the segmented data.
    instance_name : str
        The name for the instance image stack.
    instance_value : int
        The value in the semantic stack representing the class to isolate.
    min_feature_size : int
        The minimum feature size to retain in the instance image stack.

    Returns
    -------
    InstanceImageStack
        The instance image stack with isolated features.

    Examples
    --------
    >>> semantic_stack = SemanticImageStack(
    ...     name="example_semantic_stack",
    ...     metadata=MetaData(
    ...         data_volume=DataVolume(z_image_count=10, y_height=256, x_width=256),
    ...         resolution=Resolution(
    ...             dx=Length(1.0, Units.MICRON),
    ...             dy=Length(1.0, Units.MICRON),
    ...             dz=Length(1.0, Units.MICRON)
    ...         ),
    ...         pixel_units=Units.MICRON,
    ...         origin=Origin(
    ...             x0=Length(0.0, Units.MICRON),
    ...             y0=Length(0.0, Units.MICRON),
    ...             z0=Length(0.0, Units.MICRON)
    ...         )
    ...     ),
    ...     data=np.random.randint(0, 2, (10, 256, 256))
    ... )
    >>> instance_stack = semantic_to_instance(
    ...     semantic_stack,
    ...     instance_name="example_instance_stack",
    ...     instance_value=1,
    ...     min_feature_size=10
    ... )
    >>> print(instance_stack)
    InstanceImageStack(name='example_instance_stack', metadata=MetaData(data_volume=DataVolume(z_image_count=10, y_height=256, x_width=256), resolution=Resolution(dx=Length(value=1.0, unit=<Units.MICRON: 'micron'>), dy=Length(value=1.0, unit=<Units.MICRON: 'micron'>), dz=Length(value=1.0, unit=<Units.MICRON: 'micron'>)), pixel_units=<Units.MICRON: 'micron'>, origin=Origin(x0=Length(value=0.0, unit=<Units.MICRON: 'micron'>), y0=Length(value=0.0, unit=<Units.MICRON: 'micron'>), z0=Length(value=0.0, unit=<Units.MICRON: 'micron'>))), data=array(...), nlabels=..., min_feature_size=10)
    """

    # Create a boolean mask of the array
    masked_stack = semantic_stack.data == instance_value

    print(f"\tlabelling connected components in '{instance_name}'")

    # Use cc3d for connected components labeling
    cc3d_instance_stack, cc3d_nlabels = cc3d.connected_components(
        masked_stack,
        connectivity=cs.CC3DConnectivity.FACE_NEIGHBORS,
        return_N=True,
    )
    print(
        f"\t\twith cc3d package, found {cc3d_nlabels} connected components in '{instance_name}'"
    )
    return InstanceImageStack(
        name=instance_name,
        data=cc3d_instance_stack,
        metadata=semantic_stack.metadata,
        nlabels=cc3d_nlabels,
        min_feature_size=min_feature_size,
    )


##### MAIN FUNCTIONS


# TODO remove hardcoding here with yml config type
def instance_analysis_included(settings: dict, key: str) -> bool:
    """
    Checks YAML settings to determine if instance analysis should be performed.

    This function checks the provided settings dictionary to see if instance analysis
    is included for the specified key.

    Parameters
    ----------
    settings : dict
        The settings dictionary, typically loaded from a YAML file.
    key : str
        The key representing the class label to check for instance analysis inclusion.

    Returns
    -------
    bool
        True if instance analysis is included for the specified key, False otherwise.

    Examples
    --------
    >>> settings = {
    ...     "class_labels": {
    ...         "cat": {
    ...             "instance_analysis": {
    ...                 "include": True
    ...             }
    ...         },
    ...         "dog": {
    ...             "instance_analysis": {
    ...                 "include": False
    ...             }
    ...         }
    ...     }
    ... }
    >>> instance_analysis_included(settings, "cat")
    True
    >>> instance_analysis_included(settings, "dog")
    False
    >>> instance_analysis_included(settings, "bird")
    False
    """

    labels = settings["class_labels"]
    try:
        return labels[key]["instance_analysis"]["include"]
    except KeyError:
        return False


def instance_properties(
    instance_stack: InstanceImageStack,
    inst_indices: InstanceIndices,
) -> InstanceProperties:
    """
    Calculate various properties for each instance in an InstanceImageStack.

    Parameters
    ----------
    instance_stack : InstanceImageStack
        The instance image stack containing the segmented data.
    inst_indices : InstanceIndices
        The instance indices containing the indices of each instance.

    Returns
    -------
    InstanceProperties
        The properties of each instance, including number of voxels, equivalent spherical diameters,
        centroids, ellipsoids, surface areas, and volumes.

    Examples
    --------
    >>> instance_stack = InstanceImageStack(
    ...     name="example_instance_stack",
    ...     metadata=MetaData(
    ...         data_volume=DataVolume(z_image_count=10, y_height=256, x_width=256),
    ...         resolution=Resolution(
    ...             dx=Length(1.0, Units.MICRON),
    ...             dy=Length(1.0, Units.MICRON),
    ...             dz=Length(1.0, Units.MICRON)
    ...         ),
    ...         pixel_units=Units.MICRON,
    ...         origin=Origin(
    ...             x0=Length(0.0, Units.MICRON),
    ...             y0=Length(0.0, Units.MICRON),
    ...             z0=Length(0.0, Units.MICRON)
    ...         )
    ...     ),
    ...     data=np.random.randint(0, 2, (10, 256, 256)),
    ...     nlabels=3,
    ...     min_feature_size=2
    ... )
    >>> inst_indices = InstanceIndices(
    ...     source_name="example_instance_stack",
    ...     indices=[np.array([[0, 0], [1, 1]]), np.array([[0, 1], [1, 2]])],
    ...     labels=InstanceLabels(data=np.array([1, 2]))
    ... )
    >>> properties = instance_properties(instance_stack, inst_indices)
    >>> print(properties)
    InstanceProperties(source_name='example_instance_stack', labels=InstanceLabels(data=array([1, 2])), n_voxels=[NVoxel(value=0), NVoxel(value=2), NVoxel(value=2)], equivalent_sphere_diameters=[Length(value=0.0, unit=<Units.MICRON: 'micron'>), Length(value=1.240700981798799, unit=<Units.MICRON: 'micron'>), Length(value=1.240700981798799, unit=<Units.MICRON: 'micron'>)], centroids=Centroids(data=[Centroid(cx=Length(value=0.0, unit=<Units.MICRON: 'micron'>), cy=Length(value=0.0, unit=<Units.MICRON: 'micron'>), cz=Length(value=0.0, unit=<Units.MICRON: 'micron'>)), Centroid(cx=Length(value=0.5, unit=<Units.MICRON: 'micron'>), cy=Length(value=0.5, unit=<Units.MICRON: 'micron'>), cz=Length(value=0.5, unit=<Units.MICRON: 'micron'>)), Centroid(cx=Length(value=1.5, unit=<Units.MICRON: 'micron'>), cy=Length(value=1.5, unit=<Units.MICRON: 'micron'>), cz=Length(value=1.5, unit=<Units.MICRON: 'micron'>))]), ellipsoids=BestFitEllipsoids(data=[BestFitEllipsoid(a=EllipsoidAxis(length=Length(value=0.0, unit=<Units.MICRON: 'micron'>), orientation=UnitVector(u=0.0, v=0.0, w=0.0)), b=EllipsoidAxis(length=Length(value=0.0, unit=<Units.MICRON: 'micron'>), orientation=UnitVector(u=0.0, v=0.0, w=0.0)), c=EllipsoidAxis(length=Length(value=0.0, unit=<Units.MICRON: 'micron'>), orientation=UnitVector(u=0.0, v=0.0, w=0.0))]), BestFitEllipsoid(a=EllipsoidAxis(length=Length(value=1.240700981798799, unit=<Units.MICRON: 'micron'>), orientation=UnitVector(u=1.0, v=0.0, w=0.0)), b=EllipsoidAxis(length=Length(value=0.0, unit=<Units.MICRON: 'micron'>), orientation=UnitVector(u=0.0, v=1.0, w=0.0)), c=EllipsoidAxis(length=Length(value=0.0, unit=<Units.MICRON: 'micron'>), orientation=UnitVector(u=0.0, v=0.0, w=1.0))]), BestFitEllipsoid(a=EllipsoidAxis(length=Length(value=1.240700981798799, unit=<Units.MICRON: 'micron'>), orientation=UnitVector(u=1.0, v=0.0, w=0.0)), b=EllipsoidAxis(length=Length(value=0.0, unit=<Units.MICRON: 'micron'>), orientation=UnitVector(u=0.0, v=1.0, w=0.0)), c=EllipsoidAxis(length=Length(value=0.0, unit=<Units.MICRON: 'micron'>), orientation=UnitVector(u=0.0, v=0.0, w=1.0))])]), surface_areas=EllipsoidSurfaceAreas(data=[Area(value=0.0, unit_squared=<Units.MICRON: 'micron'>), Area(value=19.739208802178716, unit_squared=<Units.MICRON: 'micron'>), Area(value=19.739208802178716, unit_squared=<Units.MICRON: 'micron'>)]), volumes=EllipsoidVolumes(data=[Volume(value=0.0, unit_cubed=<Units.MICRON: 'micron'>), Volume(value=4.1887902047863905, unit_cubed=<Units.MICRON: 'micron'>), Volume(value=4.1887902047863905, unit_cubed=<Units.MICRON: 'micron'>)]))
    """

    num_labels = instance_stack.nlabels
    # feature 0 is the background, and we want to include the last label
    stack_n_voxels = num_voxels(indices=inst_indices)
    stack_eq_diameters = np.zeros(num_labels + 1, dtype=object)
    stack_centroids = np.zeros(num_labels + 1, dtype=object)
    stack_ellipsoids = np.zeros(num_labels + 1, dtype=object)
    stack_surface_areas = np.zeros(num_labels + 1, dtype=object)
    stack_volumes = np.zeros(num_labels + 1, dtype=object)

    resolution = instance_stack.metadata.resolution
    origin = instance_stack.metadata.origin

    # null types
    null_length = Length(value=0.0, unit=resolution.dx.unit)
    null_area = Area(value=0.0, unit_squared=resolution.dx.unit)
    null_volume = Volume(value=0.0, unit_cubed=resolution.dx.unit)
    null_axis = EllipsoidAxis(
        length=null_length,
        orientation=UnitVector(
            u=0.0,
            v=0.0,
            w=0.0,
        ),
    )
    # set first entry to correct type of zero for feature 0
    stack_n_voxels[0] = NVoxel(value=0)
    stack_eq_diameters[0] = null_length
    stack_centroids[0] = Centroid(
        cx=null_length,
        cy=null_length,
        cz=null_length,
    )
    stack_ellipsoids[0] = BestFitEllipsoid(
        a=null_axis,
        b=null_axis,
        c=null_axis,
    )
    stack_surface_areas[0] = null_area
    stack_volumes[0] = null_volume

    for instance in range(1, num_labels + 1):
        if ((instance - 1) % 1000) == 0:
            print(f"\t\tAnalyzing feature {instance} of {num_labels}")
        indices = inst_indices.indices[instance]
        n_voxels = stack_n_voxels[instance].value

        # Equivalent spherical diameter
        stack_eq_diameters[instance] = equivalent_spherical_diameter(
            n_voxels=n_voxels,
            resolution=resolution,
        )

        # Centroid (real-space coordinates)
        centroid = center_of_mass(
            indices=indices,
            resolution=resolution,
            origin=origin,
        )
        stack_centroids[instance] = centroid

        # Ellipsoid properties
        ellipsoid = fit_ellipsoid(
            indices=indices,
            centroid=centroid,
            resolution=resolution,
            origin=origin,
        )
        stack_ellipsoids[instance] = ellipsoid

        surface_area = ellipsoid_surface_area(ellipsoid=ellipsoid)
        stack_surface_areas[instance] = surface_area

        volume = ellipsoid_volume(ellipsoid=ellipsoid)
        stack_volumes[instance] = volume

    properties = InstanceProperties(
        source_name=instance_stack.name,
        labels=inst_indices.labels,
        n_voxels=stack_n_voxels,
        equivalent_sphere_diameters=stack_eq_diameters.tolist(),
        centroids=Centroids(data=stack_centroids.tolist()),
        ellipsoids=BestFitEllipsoids(data=stack_ellipsoids.tolist()),
        surface_areas=EllipsoidSurfaceAreas(data=stack_surface_areas.tolist()),
        volumes=EllipsoidVolumes(data=stack_volumes.tolist()),
    )

    return properties


def process(yml_file: Path) -> bool:
    """
    Processes a directory of images specified in a YAML file, creates a semantic image stack,
    and performs instance analysis on all semantic label classes within the stack.
    Saves all instance stacks and properties into an HDF5 file containing the original semantic image stack.

    Parameters
    ----------
    yml_file : Path
        The fully pathed YAML file that specifies settings.

    Returns
    -------
    bool
        True if the process was successful, False otherwise.

    Examples
    --------
    >>> yml_file = Path("path/to/input.yml")
    >>> process(yml_file)
    True
    """

    semantic_stack = process_image_stack(yml_file)
    n_semantic_labels = np.unique(semantic_stack.data)

    # create h5
    h5_path = hio.image_to_voxel(yml_file)

    # start instance analysis
    db = ut.yaml_to_dict(yml_file)
    for instance_name in db["class_labels"]:

        # # check yml values match the labels in the semantic image stack
        # for semantic_label in n_semantic_labels:

        if not instance_analysis_included(settings=db, key=instance_name):
            continue

        print(f"Analyzing '{instance_name}' semantic label")

        # create the instance_stacks
        inst_stack = semantic_to_instance(
            semantic_stack=semantic_stack,
            instance_name=instance_name,
            instance_value=db["class_labels"][instance_name]["value"],
            min_feature_size=db["class_labels"][instance_name]["instance_analysis"][
                "min_feature_size"
            ],
        )

        # get instance_indices
        inst_indices = instance_indices(
            instance_stack=inst_stack,
        )
        if inst_stack.min_feature_size > 0:
            # overwrites instance_stack and inst_indices
            inst_stack, inst_indices = minimum_size_filter(
                initial_stack=inst_stack,
                initial_indices=inst_indices,
            )

        # add to h5
        hio.add_to_h5(
            inst_stack,
            h5_path=h5_path,
            h5_group="VoxelData",
        )
        print(f"\tlabelled '{instance_name}' voxel data added")

        # TODO not implemented instance indices yet
        # hio.add_to_h5(
        #     inst_indices,
        #     h5_path=h5_path,
        #     h5_group=f"{instance_name}",
        # )

        print(f"\tcalculating '{instance_name}' instance properties...")
        # calculate indepedendent feature properties
        inst_properties = instance_properties(
            instance_stack=inst_stack,
            inst_indices=inst_indices,
        )

        # add to h5
        hio.add_to_h5(
            inst_properties,
            h5_path=h5_path,
            h5_group=f"{instance_name}",
        )
        print(f"\t'{instance_name}' instance property data added")

        nearest_neighbors = nearest_neighbor_distance(
            centroids=inst_properties.centroids
        )
        hio.add_to_h5(
            nearest_neighbors,
            h5_path=h5_path,
            h5_group=f"{instance_name}",
        )
        print(f"\t'{instance_name}' neighborhood property data added")

    return True


def main():
    """
    Command-line interface for running instance analysis on image stacks.

    This function serves as the entry point for terminal-based access to the instance analysis module.
    It is invoked from the command line using the 'instance_analysis' command specified in pyproject.toml.
    The function processes a YAML input file to create a semantic image stack and performs instance analysis
    on all semantic label classes within the stack. The results, including instance stacks and properties,
    are saved into an HDF5 file.

    Parameters
    ----------
    None

    Returns
    -------
    None

    Examples
    --------
    To run the instance analysis, use the following command in the terminal:
    $ instance_analysis path/to/input.yml
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="the .yml user input file")
    args = parser.parse_args()
    input_file = args.input_file
    process(input_file)

    print(f"\n{input_file} processed!")


if __name__ == "__main__":
    main()
