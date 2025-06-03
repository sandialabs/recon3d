# `hdf_to_npy` and `image_to_npy`

---

The `hdf_to_npy` function converts a hdf representation into a NumPy .npy file.

---

Some analyses may require the extraction of voxel data contained within an hdf file and conversion to an internal Python type. 

The `hdf_to_npy` functionality is provided as a command line interface. Providing a HDF file with voxel dataset of interest, image data can be extracted to a separate folder.

Contents of `hdf_to_npy.yml`:

```yml
<!-- cmdrun cat hdf_to_npy.yml -->
```

`hdf_to_image hdf_to_npy.yml` produces:

```sh
<!-- cmdrun hdf_to_npy hdf_to_npy.yml -->
```

---

The `image_to_npy` function converts an image stack into a NumPy .npy file.

---

Some analyses may require the extraction of image stacks to an internal Python type. 

The `image_to_npy` functionality is provided as a command line interface. Providing an image stack path with voxel dataset of interest, image data can be extracted to a separate folder.

Contents of `image_to_npy.yml`:

```yml
<!-- cmdrun cat image_to_npy.yml -->
```

`image_to_npy image_to_npy.yml` produces:

```sh
<!-- cmdrun image_to_npy image_to_npy.yml -->
```
