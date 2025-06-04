# Utilities

## `grayscale_image_stack_to_segmentation`

This utility takes

* a stack of images in a directory, encoded as 8-bit grayscale (integers `0...255`), and
* a threshold value (defaults to [Otsu](https://scikit-image.org/docs/stable/api/skimage.filters.html#skimage.filters.threshold_otsu) if no threshold value is provided),

to create

* a NumPy (`.npy`) segmentation, ready for additional image processing.

Contents of `grayscale_image_stack_to_segmentation`:

```yml
<!-- cmdrun cat grayscale_image_stack_to_segmentation.yml -->
```
