# `image_to_voxel` (and reverse)

---

The `image_to_voxel` function converts image data into a voxel representation suitable for 3D analysis, while the `voxel_to_image` function reverses this process, extracting 2D image slices from a 3D voxel dataset.

---

Some analyses may require the extraction of voxel data contained within an hdf file and conversion to an internal python type. One extension of this application is the direct conversion of a voxel array into an image stack, which may be utilized for subseqeunt downscaling and meshing activities.

Provided a HDF file with voxel dataset of interest, image data can be extracted to a separate folder (and vice versa) with the following command:

For macOS/Linux:

```sh
image_to_voxel .venv/lib/python3.11/site-packages/recon3d/examples/image_to_voxel.yml
voxel_to_image .venv/lib/python3.11/site-packages/recon3d/examples/voxel_to_image.yml
```

For Windows:

```sh
image_to_voxel .venv\Lib\site-packages\recon3d\examples\image_to_voxel.yml
voxel_to_image .venv\Lib\site-packages\recon3d\examples\voxel_to_image.yml
```

and shown below, provides the recipe for the ``image_to_voxel`` and ``voxel_to_image`` process.

```yml
<!-- cmdrun cat image_to_voxel.yml -->
```

```yml
<!-- cmdrun cat voxel_to_image.yml -->
```

Images will be written out to a folder with the same "stem" name as the internal `voxel_data_location` key, in this case "machined_cylinder".

The program output follows on Windows:

```sh
Success: database created from file: .venv\Lib\site-packages\recon3d\examples\voxel_to_image.yml
key, value, type
---, -----, ----
cli_entry_points, ['voxel_to_image'], <class 'list'>
voxel_data_path, recon3d/data/machined_cylinder.h5, <class 'str'>    
voxel_data_location, VoxelData/machined_cylinder, <class 'str'>      
image_parent_dir, recon3d/data/output/, <class 'str'>
image_slice_normal, Z, <class 'str'>
image_output_type, .tif, <class 'str'>

Voxel data extracted to the following relative directory:'recon3d\data\output\machined_cylinder'
```