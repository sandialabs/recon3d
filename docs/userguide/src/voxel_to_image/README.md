# `image_to_voxel` (and reverse)

---

The `image_to_voxel` function converts image data into a voxel representation suitable for 3D analysis, while the `voxel_to_image` function reverses this process, extracting 2D image slices from a 3D voxel dataset.

---

Some analyses may require the extraction of voxel data contained within an hdf file and conversion to an internal python type. One extension of this application is the direct conversion of a voxel array into an image stack, which may be utilized for subsequent downscaling and meshing activities.

With ``recon3d`` installed in a virtual environment called `.venv`, the `image_to_voxel` and `voxel_to_image` functionality is provided as a command line interface. Providing a HDF file with voxel dataset of interest, image data can be extracted to a separate folder (and vice versa).

Contents of `image_to_voxel.yml`:

```yml
<!-- cmdrun cat image_to_voxel.yml -->
```

Contents of `voxel_to_image.yml`:

```yml
<!-- cmdrun cat voxel_to_image.yml -->
```

`image_to_voxel image_to_voxel.yml` produces:

```sh
<!-- cmdrun image_to_voxel image_to_voxel.yml -->
```

`voxel_to_image voxel_to_image.yml` produces:

```sh
<!-- cmdrun voxel_to_image voxel_to_image.yml -->
```
