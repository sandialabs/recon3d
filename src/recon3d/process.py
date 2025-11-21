"""
Basic image processing functions for 3D volumes.
"""

from typing import Tuple, Sequence

import numpy as np


def histogram_stretch(
    vol: np.ndarray,
    clip_percentiles: Tuple[float, float] = (1.0, 99.0),
    ignore_values: Sequence[int] = (0,),
) -> np.ndarray:
    """
    Simple percentile‐based contrast stretch of a volume, optionally
    ignoring certain values when computing the percentiles.

    Parameters
    ----------
    vol : np.ndarray
      Input image volume (3D or 4D).
    clip_percentiles : (float,float)
      e.g. (1,99) to compute the 1st and 99th percentile.
    ignore_values : sequence of ints
      Voxels equal to any of these are omitted from the percentile calc.
      By default only [0] is ignored; 255 (or any other value) will be
      included unless you also put it in this list.

    Returns
    -------
    stretched : np.ndarray
      Same shape/dtype as `vol`, linearly remapped so that `lo-->min`
      and `hi-->max`, clipped into the full dynamic range of the dtype.
    """

    flat = vol.ravel()

    # build a mask: True for voxels we *do* want in our percentile calc
    if ignore_values:
        mask = ~np.isin(flat, ignore_values)
    else:
        mask = np.ones_like(flat, dtype=bool)

    # if after masking there's nothing left, just return a copy
    if not mask.any():
        return vol.copy()

    # compute lo/hi percentiles on the *masked* data
    lo, hi = np.percentile(flat[mask], clip_percentiles)

    # figure out full output range for this dtype
    dtype = vol.dtype
    if np.issubdtype(dtype, np.integer):
        info = np.iinfo(dtype)
    else:
        info = np.finfo(dtype)
    out_min, out_max = info.min, info.max

    # do the linear stretch
    volf = vol.astype(np.float32)
    volf = (volf - lo) / (hi - lo)
    volf = np.clip(volf, 0.0, 1.0)
    volf = volf * (out_max - out_min) + out_min

    return volf.astype(dtype)
