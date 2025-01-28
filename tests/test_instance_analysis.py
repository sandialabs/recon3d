"""This module tests individual functions for the core module.

Example:
    To run
    - activate the virutal environment, then
    $ cd ~/recond3d

    pytest (without code coverage)
    pytest recon3d/tests/test_instance_analysis.py -v  # -v is for verbose

    to run just a single test in this module, for example
    pytest recon3d/tests/test_instance_analysis.py::test_read_images -v

"""

# Python standard libarires
import glob
from pathlib import Path

# Third party libraries
import numpy as np
import pytest

# local libraries
# import recon3d.feature_analysis as fa
import recon3d.instance_analysis as ia
import recon3d.utility as ut
import recon3d.constants as cs
from recon3d.types import *

from recon3d.static_test_paths import *


def create_instance_stack(instance_name):
    """Generate object required for further tests"""

    db = ut.yaml_to_dict(INSTANCE_ANALYSIS_YML)

    # overwrite so that this test runs on any machine
    db["semantic_image_dir"] = str(
        Path(__file__).parent.joinpath("data", "cylinder_machined_semantic")
    )
    temp_path = Path(__file__).parent.joinpath("data", "temp.yml")

    # make sure the recon3d/examples/temp.yml doesn't already
    # exist from a previous test run or from a test paused
    # midway through
    try:
        temp_path.unlink()  # delete the temp.yml if it exists aready
    except FileNotFoundError:
        pass

    # make temp yml
    bb = ut.dict_to_yaml(db, str(temp_path))

    found_stack = ia.process_image_stack(yml_input_file=bb)
    print(db)
    # instance_name = "pore"
    instance_value = db["class_labels"][instance_name]["value"]
    min_feat_size = db["class_labels"][instance_name]["instance_analysis"][
        "min_feature_size"
    ]

    found_instance = ia.semantic_to_instance(
        semantic_stack=found_stack,
        instance_name=instance_name,
        instance_value=instance_value,
        min_feature_size=min_feat_size,
    )

    temp_path.unlink()

    return found_instance


# @pytest.fixture
# def test_file() -> Path:
#     """The relative path and file string locating the default yml test file."""

#     return Path(__file__).parent.joinpath("examples", "instance_analysis.yml")


def test_process_image_stack():
    """Test that a yml input file that specifies an image stack can be read in
    and used to create a h5 file.
    """

    aa = INSTANCE_ANALYSIS_YML
    db = ut.yaml_to_dict(aa)
    # overwrite so that this test runs on any machine
    db["semantic_image_dir"] = str(
        Path(__file__).parent.joinpath("data", "cylinder_machined_semantic")
    )

    temp_path = Path(__file__).parent.joinpath("data", "temp.yml")

    try:
        temp_path.unlink()  # delete the temp.yml if it exists aready
    except FileNotFoundError:
        pass

    bb = ut.dict_to_yaml(db, str(temp_path))

    found_stack = ia.process_image_stack(yml_input_file=bb)

    images_path = Path(__file__).parent.joinpath("data", "cylinder_machined_semantic")

    known_stack = SemanticImageStack(
        name="semantic_output",
        data=ut.read_images(images_path, "tif"),
        metadata=MetaData(
            data_volume=DataVolume(
                x_width=923,
                y_height=925,
                z_image_count=21,
            ),
            resolution=Resolution(
                dx=Length(1, Units.MICRON),
                dy=Length(1, Units.MICRON),
                dz=Length(1, Units.MICRON),
            ),
            pixel_units="micron",
            origin=Origin(
                x0=Length(0.0, Units.MICRON),
                y0=Length(0.0, Units.MICRON),
                z0=Length(0.0, Units.MICRON),
            ),
        ),
    )

    assert found_stack.data.all() == known_stack.data.all()

    temp_path.unlink()


def test_read_images():
    """Tests that either images exist and can be read, or that None is returned."""

    dataset_path = Path(__file__).parent.joinpath("data", "cylinder_machined_semantic")
    image_type = ".tif"

    output = ut.read_images(file_dir=dataset_path, file_type=image_type)

    assert isinstance(output, np.ndarray)

    # now specify some non-existed files
    image_type = (
        ".gif"  # overwrite, this file type does not exist in the dataset_path folder.
    )
    with pytest.raises(FileNotFoundError):
        output = ut.read_images(file_dir=dataset_path, file_type=image_type)


def test_semantic_to_instance():
    """Test the instance image stack is properly extracted from the semantic image stack."""

    found_instance = create_instance_stack("pore")

    assert found_instance.nlabels == 103


def test_instance_mask():
    """Test the instance mask is properly created"""
    instance_stack = create_instance_stack("pore")
    instance_ids = ia.instance_indices(instance_stack)
    indices = instance_ids.indices
    labels = instance_ids.labels.data
    found_idxs = indices[1]

    known_idxs = np.array([[1, 14, 52], [1, 13, 52]])
    assert len(labels) == len(np.unique(instance_stack.data))
    assert found_idxs == pytest.approx(known_idxs)


def test_centroid():
    """Test the centroid is properly calcualted"""
    instance_stack = create_instance_stack("pore")
    instance_ids = ia.instance_indices(instance_stack)
    indices = instance_ids.indices
    found_idxs = indices[1]
    resolution = Resolution(
        dx=Length(2.0, unit=Units.MICRON),
        dy=Length(2.0, unit=Units.MICRON),
        dz=Length(2.0, unit=Units.MICRON),
    )
    origin = Origin(
        x0=Length(1.0, unit=Units.MICRON),
        y0=Length(1.0, unit=Units.MICRON),
        z0=Length(1.0, unit=Units.MICRON),
    )
    found_centroid = ia.center_of_mass(
        indices=found_idxs, resolution=resolution, origin=origin
    )

    known_centroid = Centroid(
        cx=Length((52 * 2.0) + 1.0, unit=Units.MICRON),
        cy=Length((13.5 * 2.0) + 1.0, unit=Units.MICRON),
        cz=Length((1 * 2.0) + 1.0, unit=Units.MICRON),
    )
    assert found_centroid == pytest.approx(known_centroid)


def test_centroid_pix():
    "test the centroid in pixels/voxels units is properly calculated"
    centroid = Centroid(
        cx=Length(5.2, unit=Units.MICRON),
        cy=Length(1.1, unit=Units.MICRON),
        cz=Length(3.5, unit=Units.MICRON),
    )
    resolution = Resolution(
        dx=Length(2.0, unit=Units.MICRON),
        dy=Length(2.0, unit=Units.MICRON),
        dz=Length(2.0, unit=Units.MICRON),
    )
    origin = Origin(
        x0=Length(1.0, unit=Units.MICRON),
        y0=Length(1.0, unit=Units.MICRON),
        z0=Length(1.0, unit=Units.MICRON),
    )
    found_centroid_px = ia.centroid_pix(
        centroid=centroid, resolution=resolution, origin=origin
    )

    known_centroid_px = Centroid(
        cx=Length(2.1, Units.VOXEL),
        cy=Length(0.05, Units.VOXEL),
        cz=Length(1.25, Units.VOXEL),
    )

    assert found_centroid_px.cx.value == pytest.approx(known_centroid_px.cx.value)
    assert found_centroid_px.cy.value == pytest.approx(known_centroid_px.cy.value)
    assert found_centroid_px.cz.value == pytest.approx(known_centroid_px.cz.value)


def test_calc_moment():
    "test the calc_moment properly calculated"
    instance_stack = create_instance_stack("pore")
    instance_ids = ia.instance_indices(instance_stack)
    indices = instance_ids.indices
    found_idxs = indices[73]

    resolution = Resolution(
        dx=Length(2.0, Units.MICRON),
        dy=Length(2.0, Units.MICRON),
        dz=Length(2.0, Units.MICRON),
    )
    origin = Origin(
        x0=Length(1.0, Units.MICRON),
        y0=Length(1.0, Units.MICRON),
        z0=Length(1.0, Units.MICRON),
    )
    found_centroid = ia.center_of_mass(
        indices=found_idxs, resolution=resolution, origin=origin
    )
    found_centroid_px = ia.centroid_pix(
        centroid=found_centroid, resolution=resolution, origin=origin
    )
    found_moment = ia.calc_moment(
        indices=found_idxs, p=2, q=2, r=2, centroid_px=found_centroid_px
    )
    known_moment = (0.5**2) ** 3
    assert found_moment == pytest.approx(known_moment)


def test_fit_ellipsoid():
    "test the fit_ellipsoid properly calculated"
    instance_stack = create_instance_stack("pore")
    instance_ids = ia.instance_indices(instance_stack)
    indices = instance_ids.indices
    found_idxs = indices[1]
    resolution = Resolution(
        dx=Length(2.0, unit=Units.MICRON),
        dy=Length(2.0, unit=Units.MICRON),
        dz=Length(2.0, unit=Units.MICRON),
    )
    origin = Origin(
        x0=Length(1.0, unit=Units.MICRON),
        y0=Length(1.0, unit=Units.MICRON),
        z0=Length(1.0, unit=Units.MICRON),
    )
    found_centroid = ia.center_of_mass(
        indices=found_idxs, resolution=resolution, origin=origin
    )

    found_ellipsoid = ia.fit_ellipsoid(
        indices=found_idxs,
        centroid=found_centroid,
        resolution=resolution,
        origin=origin,
    )

    a_vector = UnitVector(u=0.0, v=1.0, w=0.0)
    a = EllipsoidAxis(
        length=Length(2.23607, unit=Units.MICRON),
        orientation=a_vector,
    )
    b_vector = UnitVector(
        u=0.0,
        v=0.0,
        w=1.0,
    )
    b = EllipsoidAxis(
        length=Length(0.0, unit=Units.MICRON),
        orientation=b_vector,
    )
    c_vector = UnitVector(
        u=1.0,
        v=0.0,
        w=0.0,
    )
    c = EllipsoidAxis(
        length=Length(0.0, unit=Units.MICRON),
        orientation=c_vector,
    )
    known_ellipsoid = BestFitEllipsoid(a=a, b=b, c=c)

    # TODO look into __repr__ method to directly compare BestFitEllipsoid obejcts
    assert found_ellipsoid.a.length.value == pytest.approx(
        known_ellipsoid.a.length.value
    )
    assert found_ellipsoid.b.length.value == pytest.approx(
        known_ellipsoid.b.length.value
    )
    assert found_ellipsoid.c.length.value == pytest.approx(
        known_ellipsoid.c.length.value
    )


def test_equivalent_spherical_diameter():
    """test calculation of equivalent diameter"""
    num_vox = 4
    res = Resolution(
        dx=Length(2.0, unit=Units.MICRON),
        dy=Length(2.0, unit=Units.MICRON),
        dz=Length(2.0, unit=Units.MICRON),
    )
    known_diam = 3.93898
    found_diam = ia.equivalent_spherical_diameter(
        n_voxels=num_vox,
        resolution=res,
    )
    assert known_diam == pytest.approx(found_diam.value)

    num_vox2 = 27
    res2 = Resolution(
        dx=Length(1.1, unit=Units.MICRON),
        dy=Length(1.1, unit=Units.MICRON),
        dz=Length(2.0, unit=Units.MICRON),
    )
    known_diam2 = 4.9972
    found_diam2 = ia.equivalent_spherical_diameter(
        n_voxels=num_vox2,
        resolution=res2,
    )
    assert known_diam2 == pytest.approx(found_diam2.value)


def test_num_voxels():
    """test calculation of number of voxels of instance IDs"""
    instance_stack = create_instance_stack("pore")
    instance_ids = ia.instance_indices(instance_stack)
    n_voxels = ia.num_voxels(indices=instance_ids)
    data_vol = instance_stack.metadata.data_volume
    total_voxels = data_vol.x_width * data_vol.y_height * data_vol.z_image_count

    found_n_vox = [i.value for i in n_voxels]
    assert sum(found_n_vox) == total_voxels
    assert n_voxels[3].value == 4


def test_minimum_size_filter():
    """test application of minimum size filter"""
    instance_stack = create_instance_stack("pore")
    instance_ids = ia.instance_indices(instance_stack)
    filtered_stack, filtered_ids = ia.minimum_size_filter(
        instance_stack,
        instance_ids,
    )
    assert len(filtered_ids.labels.data) == 6
    assert filtered_stack.nlabels == 5
    assert len(filtered_ids.indices[1]) == 21

    data_vol = filtered_stack.metadata.data_volume
    known_voxels = data_vol.x_width * data_vol.y_height * data_vol.z_image_count

    idx = filtered_ids.indices

    # list comprehension on list of arrays doesn't seem to work
    # found_voxels = [len(idx[i]) for i in idx]

    found_voxels = 0
    for i in idx:
        found_voxels += len(i)
    assert known_voxels == found_voxels


def test_instance_analysis_include():
    """test whether instance analysis should be ran"""
    aa = INSTANCE_ANALYSIS_YML
    db = ut.yaml_to_dict(aa)
    bb = ia.instance_analysis_included(settings=db, key="air")
    assert bb is False

    cc = ia.instance_analysis_included(settings=db, key="metal")
    assert cc is False

    dd = ia.instance_analysis_included(settings=db, key="pore")
    assert dd is True

    ee = ia.instance_analysis_included(settings=db, key="rand")
    assert ee is False


def test_instance_properties():
    instance_stack = create_instance_stack("pore")
    instance_ids = ia.instance_indices(instance_stack)

    found_properties = ia.instance_properties(
        instance_stack=instance_stack, inst_indices=instance_ids
    )
    assert found_properties.source_name == "pore"
    assert len(found_properties.labels.data) == 104
    assert found_properties.n_voxels[55].value == 4
    assert found_properties.ellipsoids.data[53].a.length.value == pytest.approx(
        1.118034
    )

    known_centroid = Centroid(
        cx=Length(12.5, Units.MICRON),
        cy=Length(81.5, Units.MICRON),
        cz=Length(14.0, Units.MICRON),
    )
    found_centroid = found_properties.centroids.data[68]

    assert found_centroid == known_centroid

    # rod_demo_file = Path(__file__).parent.joinpath(
    #     "examples", "rods_instance_analysis.yml"
    # )
    # instance_stack = create_instance_stack(rod_demo_file, "rods")
    # instance_ids = ia.instance_indices(instance_stack)
    # found_properties = ia.instance_properties(
    #     instance_stack=instance_stack, inst_indices=instance_ids
    # )

    # assert found_properties.ellipsoids.data[4].a.orientation.u == pytest.approx(
    #     -0.584283594756
    # )


def test_nth_nearest_neighbor_distance():
    """test calculation of nth nearest neighbors"""
    instance_stack = create_instance_stack("pore")
    instance_indices = ia.instance_indices(instance_stack)
    instance_props = ia.instance_properties(
        instance_stack=instance_stack,
        inst_indices=instance_indices,
    )

    nearest_neighbors = ia.nearest_neighbor_distance(centroids=instance_props.centroids)

    assert nearest_neighbors.instance_id[7] == 12
    assert nearest_neighbors.distances[7].value == pytest.approx(6.576473218982953)


def test_ellipsoid_surface_area():
    """test calculation of ellipsoid surface area"""
    ellipsoid = BestFitEllipsoid(
        a=EllipsoidAxis(
            length=Length(value=2, unit=Units.METER),
            orientation=UnitVector(u=1.0, v=0.0, w=0.0),
        ),
        b=EllipsoidAxis(
            length=Length(value=3, unit=Units.METER),
            orientation=UnitVector(u=0.0, v=1.0, w=0.0),
        ),
        c=EllipsoidAxis(
            length=Length(value=4, unit=Units.METER),
            orientation=UnitVector(u=0.0, v=0.0, w=1.0),
        ),
    )

    found_surface_area = ia.ellipsoid_surface_area(ellipsoid=ellipsoid)
    known_surface_area = 111.604
    assert found_surface_area.value == pytest.approx(known_surface_area)
    assert found_surface_area.unit_squared == Units.METER
