# `hdf_to_image`

---

The `hdf_to_image` function  converts a hdf representation into a stack of images.

---

Some analyses may require the extraction of voxel data contained within an hdf file and conversion to an internal Python type. One extension of this application is the direct conversion of a voxel array into an image stack, which may be utilized for subsequent downscaling and meshing activities.

The `hdf_to_image` functionality is provided as a command line interface. Providing a HDF file with voxel dataset of interest, image data can be extracted to a separate folder.

Contents of `hdf_to_image.yml`:

```yml
<!-- cmdrun cat hdf_to_image.yml -->
```

`hdf_to_image hdf_to_image.yml` produces:

```sh
<!-- cmdrun hdf_to_image hdf_to_image.yml -->
```
