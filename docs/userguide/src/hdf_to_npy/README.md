# `hdf_to_npy`

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
