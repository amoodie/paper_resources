import os

import numpy as np


def define_time_idx_set(cube, start=0, num=50, irregular=False):
    if not irregular:
        _ti = np.linspace(start, cube.shape[0]-1,
                          num=num, endpoint=True, dtype=int)
    else:
        _ti0 = np.logspace(0, np.log10(cube.shape[0]-1), num=num, dtype=np.int)
        _ti = np.unique(np.concatenate((np.array([0]), _ti0)))
    return _ti


def _hacky_make_masks(cubes):
    _masks = []
    for c, cube in enumerate(cubes):
        _sigma = cube.meta.sigma
        _mask = _hacky_make_mask(_sigma)
        # _masks[c, :, :] = _mask
        if np.any(_mask):
            _masks.append(_mask)

    _masks = np.array(_masks)

    np.save(
        os.path.join('_analysis', 'selenga_subsd_mask_set.npy'),
        _masks)


def _hacky_make_mask(_sigma, save=False):
    _mask = (_sigma > 0).astype(bool)
    if save:
        np.save(
            os.path.join('_analysis', 'fixed_subsd_mask.npy'),
            _mask)
    else:
        return _mask


def _hacky_get_mask():
    _mask_set = np.load(os.path.join('_analysis', 'subsd_mask_set.npy'))
    _N = _mask_set.shape[0]
    _I = np.random.randint(0, _N)
    _mask = _mask_set[_I, :, :].squeeze()
    return _mask
