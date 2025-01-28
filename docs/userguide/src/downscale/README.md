# `downscale`

---

``downscale`` pads an image stack to make its dimensions evenly divisible by a target resolution, then downscales it to that resolution, effectively reducing the number of voxels for computations like Finite Element Analysis, resulting in a significant reduction in image resolution while adding padding to maintain dimension compatibility.

---

With ``recon3d`` installed in a virtual environment called ``.venv``, the ``downscale`` functionality is provided as a command line interface.  Following is an example of the downscale workflow.  

The subject ``thunder_gray.tif`` file, installed to

* developer: ``recon3d/recon3d/data/thunder_gray/thunder_gray.tif``,
* client (macOS/Linux): ``.venv/lib/python3.11/site-packages/recon3d/data/thunder_gray/thunder_gray.tif``,
* client (Windows): ``.venv\Lib\site-packages\recon3d\data\thunder_gray\thunder_gray.tif``

is shown below: 

![thunder_gray_input_illustrated.png](thunder_gray_input_illustrated.png)

The subject file, a grayscale image with ``734 × 348`` of pixel resolution, is shown in the Graphic app to demonstrate pixel height and width.

The ``downscale_thunder.yml`` file, installed to

* developer: ``recon3d/recon3d/examples/downscale_thunder.yml``,
* client (macOS/Linux): ``.venv/lib/python3.11/site-packages/recon3d/examples/downscale_thunder.yml``,
* client (Windows): ``.venv\Lib\site-packages\recon3d\examples\downscale_thunder.yml``,

and shown below, provides the recipe for the ``downscale`` process.

```yml
cli_entry_points:
- downscale  
downscale_tolerance: 0.0001
image_dir: recon3d/data/thunder_gray
image_limit_factor: 2.0
image_type: .tif
out_dir: recon3d/data/output # output path for the processed files
output_stack_type: padded
padding:
    nx: 10
    ny: 10
    nz: 0
resolution_input:
    dx: 10.0
    dy: 10.0
    dz: 10.0
resolution_output:
    dx: 50.0
    dy: 50.0
    dz: 10.0
save_npy: true
writeVTR: true
```

The image is processed with the following command line interface:

For macOS/Linux:

```
downscale .venv/lib/python3.11/site-packages/recon3d/examples/downscale_thunder.yml
```

For Windows:

```
downscale .venv\Lib\site-packages\recon3d\examples\downscale_thunder.yml
```

The program output follows:

```sh
Processing file: recon3d/examples/downscale_thunder.yml
Success: database created from file: recon3d/examples/downscale_thunder.yml
key, value, type
---, -----, ----
cli_entry_points, ['downscale'], <class 'list'>
downscale_tolerance, 0.0001, <class 'float'>
image_dir, recon3d/data/thunder_gray, <class 'str'>
image_limit_factor, 2.0, <class 'float'>
image_type, .tif, <class 'str'>
out_dir, recon3d/data/output, <class 'str'>
output_stack_type, padded, <class 'str'>
padding, {'nx': 10, 'ny': 10, 'nz': 0}, <class 'dict'>
resolution_input, {'dx': 10.0, 'dy': 10.0, 'dz': 10.0}, <class 'dict'>
resolution_output, {'dx': 50.0, 'dy': 50.0, 'dz': 10.0}, <class 'dict'>
save_npy, True, <class 'bool'>
writeVTR, True, <class 'bool'>
Images read, image array size: (1, 348, 734)
Original array size: (1, 348, 734)
New array size: (1, 350, 735)
Downscaled array size: (1, 70, 147)
Cropping image stack to bounding box...
Saving cropped image stack in recon3d/data/output > images_at_resolution_50_dx
.npy file saved.
.vtr file saved.
Finished processing file: recon3d/examples/downscale_thunder.yml
```

The output file(s) image stack is saved to ``recon3d/data/output/images_at_resolution_50_dx``.
The output file ``recon3d/data/output/images_at_resolution_50_dx/0000.tif`` is shown below: 

![thunder_gray_output_illustrated.png](thunder_gray_output_illustrated.png)
  
The output file ``recon3d/data/output/images_at_resolution_50_dx/0000.tif``, a grayscale image with ``734 × 348`` to ``147x70`` (``5x`` reduction) plus ``10+10`` padding to ``167x90`` of pixel resolution, is shown atop the subject image in the Graphic app to demonstrate pixel height and width.  Note the black border padding, with of ``10`` pixels.

The ``50_dx.npy`` and ``50_dx.vtr`` files are also saved to ``recon3d/data/output``.
