# `instance_analysis`

---

`instance_analysis` takes a semantic image stack and generates an instance image stack for one semantic label at a time, as well as statistics on those instaces, such as size, shape from a best-fit ellispoid, as well as neighborhood metrics. 

---

As opposed to semantic segmentation, which groups objects in an image based on defined categories, instance segmentation can be considered a refined version of semantic segmentation wherein individual instances are independently labelled. For example, an image can be semantically segmented to distinguish between different animals (person, sheep, or dog) in an image, while an instance segmentation will label the various individual animals in the image:
Instance analysis provides the user with various metrics for each *instance* of a **semantic label**, such as the size, shape, or orientation of each individual *instance*, enabling a statistical assessment of instance populations from a dataset. This is typically done on individual populations of semantic labels, so an image with multiple sheep and multiple dogs would return separate data for both the sheep and dog populations. Instance analysis on combined semantic labels would require the generation of a new semantic label (e.g. four-legged animals).

![semantic_vs_instance.png](semantic_vs_instance.png)

With ``recon3d`` installed in a virtual environment called ``.venv``, the ``instance_analysis`` functionality is provided as a command line interface.  Following is an example of the instance_analysis workflow.

Provided a semantic image stack with labeled classes, the instances of each class and associated instance properties can be generated using: 

For macOS/Linux:

```sh
instance_analysis .venv/lib/python3.11/site-packages/recon3d/examples/instance_analysis.yml
```

For Windows:

```sh
instance_analysis .venv\Lib\site-packages\recon3d\examples\instance_analysis.yml
```

and shown below, provides the recipe for the ``instance_analysis`` process.

```yml
<!-- cmdrun cat instance_analysis.yml -->

```

The program output follows on Linux:

```sh
<!-- cmdrun pwd -->
<!-- cmdrun instance_analysis ./docs/userguide/src/instance_analysis/instance_analysis.yml -->
```

```sh
Processing specification file: .venv\Lib\site-packages\recon3d\examples\instance_analysis.yml
Success: database created from file: .venv\Lib\site-packages\recon3d\examples\instance_analysis.yml
key, value, type
---, -----, ----
cli_entry_points, ['instance_analysis', 'image_to_voxel'], <class 'list'>
semantic_images_dir, recon3d/data/cylinder_machined_semantic, <class 'str'>
semantic_images_type, .tif, <class 'str'>
semantic_imagestack_name, machined_cylinder, <class 'str'>
class_labels, {'air': {'value': 0, 'instance_analysis': {'include': False, 'min_feature_size': 0}}, 'metal': {'value': 1, 'instance_analysis': {'include': False, 'min_feature_size': 0}}, 'pores': {'value': 2, 'instance_analysis': {'include': True, 'min_feature_size': 12}}}, <class 'dict'>
voxel_size, {'dx': 1.0, 'dy': 1.0, 'dz': 1.0}, <class 'dict'>        
origin, {'x0': 0.0, 'y0': 0.0, 'z0': 0.0}, <class 'dict'>
pixel_units, micron, <class 'str'>
out_dir, recon3d/data/output, <class 'str'>
h5_filename, cylinder_machined, <class 'str'>
Input path: recon3d\data\cylinder_machined_semantic
Output path: recon3d\data\output
Images read, image array size: (21, 150, 150)
image stack, machined_cylinder, has dimensions (num_images, row, col): (21, 150, 150)
Success: database created from file: .venv\Lib\site-packages\recon3d\examples\instance_analysis.yml
key, value, type
---, -----, ----
cli_entry_points, ['instance_analysis', 'image_to_voxel'], <class 'list'>
semantic_images_dir, recon3d/data/cylinder_machined_semantic, <class 'str'>
semantic_images_type, .tif, <class 'str'>
semantic_imagestack_name, machined_cylinder, <class 'str'>
class_labels, {'air': {'value': 0, 'instance_analysis': {'include': False, 'min_feature_size': 0}}, 'metal': {'value': 1, 'instance_analysis': {'include': False, 'min_feature_size': 0}}, 'pores': {'value': 2, 'instance_analysis': {'include': True, 'min_feature_size': 12}}}, <class 'dict'>
voxel_size, {'dx': 1.0, 'dy': 1.0, 'dz': 1.0}, <class 'dict'>        
origin, {'x0': 0.0, 'y0': 0.0, 'z0': 0.0}, <class 'dict'>
pixel_units, micron, <class 'str'>
out_dir, recon3d/data/output, <class 'str'>
h5_filename, cylinder_machined, <class 'str'>
Processing specification file: .venv\Lib\site-packages\recon3d\examples\instance_analysis.yml
Success: database created from file: .venv\Lib\site-packages\recon3d\examples\instance_analysis.yml
key, value, type
---, -----, ----
cli_entry_points, ['instance_analysis', 'image_to_voxel'], <class 'list'>
semantic_images_dir, recon3d/data/cylinder_machined_semantic, <class 'str'>
semantic_images_type, .tif, <class 'str'>
semantic_imagestack_name, machined_cylinder, <class 'str'>
class_labels, {'air': {'value': 0, 'instance_analysis': {'include': False, 'min_feature_size': 0}}, 'metal': {'value': 1, 'instance_analysis': {'include': False, 'min_feature_size': 0}}, 'pores': {'value': 2, 'instance_analysis': {'include': True, 'min_feature_size': 12}}}, <class 'dict'>
voxel_size, {'dx': 1.0, 'dy': 1.0, 'dz': 1.0}, <class 'dict'>        
origin, {'x0': 0.0, 'y0': 0.0, 'z0': 0.0}, <class 'dict'>
pixel_units, micron, <class 'str'>
out_dir, recon3d/data/output, <class 'str'>
h5_filename, cylinder_machined, <class 'str'>
Input path: recon3d\data\cylinder_machined_semantic
Output path: recon3d\data\output
Images read, image array size: (21, 150, 150)
image stack, machined_cylinder, has dimensions (num_images, row, col): (21, 150, 150)
Wrote output file: recon3d\data\output\cylinder_machined.h5
Success: database created from file: .venv\Lib\site-packages\recon3d\examples\instance_analysis.yml
key, value, type
---, -----, ----
cli_entry_points, ['instance_analysis', 'image_to_voxel'], <class 'list'>
semantic_images_dir, recon3d/data/cylinder_machined_semantic, <class 'str'>
semantic_images_type, .tif, <class 'str'>
semantic_imagestack_name, machined_cylinder, <class 'str'>
class_labels, {'air': {'value': 0, 'instance_analysis': {'include': False, 'min_feature_size': 0}}, 'metal': {'value': 1, 'instance_analysis': {'include': False, 'min_feature_size': 0}}, 'pores': {'value': 2, 'instance_analysis': {'include': True, 'min_feature_size': 12}}}, <class 'dict'>
voxel_size, {'dx': 1.0, 'dy': 1.0, 'dz': 1.0}, <class 'dict'>        
origin, {'x0': 0.0, 'y0': 0.0, 'z0': 0.0}, <class 'dict'>
pixel_units, micron, <class 'str'>
out_dir, recon3d/data/output, <class 'str'>
h5_filename, cylinder_machined, <class 'str'>
Analyzing 'pores' semantic label
        labelling connected components in 'pores'
                with cc3d pacakage, found 103 connected components in 'pores'
                5 connected components remaining after minimum size filter of 12 voxels
        labelled 'pores' voxel data added
        calculating 'pores' instance properties...
                Analyzing feature 1 of 5
        'pores' instance property data added
        'pores' neighborhood property data added

.venv\Lib\site-packages\recon3d\examples\instance_analysis.yml processed!
```

This always outputs an `hdf` file.  The current version of `hdf` file is `hdf5`, so file extensions will terminate in `.h5`. The `hdf` file output can be opened in the HDFView application.
