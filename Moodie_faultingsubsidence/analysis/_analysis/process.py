import numpy as np
import pandas as pd

import skimage
from sklearn.linear_model import LinearRegression


import deltametrics as dm
from pyDeltaRCM.shared_tools import scale_model_time

from . import utilities

eta_array = np.empty((100, 100, 100))


# use deltametrics to make masks
def _make_OAP(_eta):
    _OAP = dm.plan.OpeningAnglePlanform.from_elevation_data(
        _eta,
        elevation_threshold=0)
    return _OAP


def _make_OAP_eta_t(t):
    global eta_array
    _tth_eta = eta_array[t, :, :]
    return _make_OAP(_tth_eta)


def _make_LandMask(_OAP):
    return dm.mask.LandMask.from_OAP(_OAP)


def _make_FlowMask(_vel):
    return dm.mask.FlowMask(
        _vel,
        flow_threshold=0.3)


def _make_ChannelMask(_OAP, _FlowMask):
    return dm.mask.ChannelMask.from_OAP_and_FlowMask(
        _OAP, _FlowMask)


def _make_ShorelineMask(_OAP):
    return dm.mask.ShorelineMask.from_OAP(_OAP)


def _make_SubsidenceMask(cube):
    if 'sigma' in cube.meta.keys():
        _mask = np.array(cube.meta['sigma']) != 0
        if np.sum(_mask) == 0:
            _mask = utilities._hacky_get_mask()
    else:
        _mask = utilities._hacky_get_mask()
    return _mask


def make_OAP_list(cube, time_idx, pool):
    """Make a list of the OAPs for a single cube.

    Takes a pool as an argument, which is used for processing.
    """
    eta_as_array = cube.export_frozen_variable('eta')
    _eta_list = []
    for t in time_idx:
        _eta_list.append(eta_as_array[t, :, :])
    _OAP_list = pool.starmap(
        _make_OAP,
        zip(_eta_list))

    return _OAP_list


def make_mask_arrays(cube, OAP_list, time_idx):
    """Make the mask arrays.

    This processess again loops through the time slices of a single cube, and
    uses the matching OAP list to create various masks.
    """
    _landmask_list = []
    _flowmask_list = []
    _channelmask_list = []
    _shorelinemask_list = []

    for i, t in enumerate(time_idx):
        _LM = _make_LandMask(OAP_list[i])
        _FM = _make_FlowMask(cube['velocity'][t, :, :])
        _CM = _make_ChannelMask(OAP_list[i], _FM)
        _SM = _make_ShorelineMask(OAP_list[i])

        _landmask_list.append(_LM.mask)
        _flowmask_list.append(_FM.mask)
        _channelmask_list.append(_CM.mask)
        _shorelinemask_list.append(_SM.mask)

    _landmask_array = np.stack(_landmask_list)
    _channelmask_array = np.stack(_channelmask_list)
    _shorelinemask_array = np.stack(_shorelinemask_list)

    _subsidencemask = _make_SubsidenceMask(cube)

    return _landmask_array, _channelmask_array, _shorelinemask_array, _subsidencemask  # noqa: E501


def extract_single_metadata(_path, custom):
    """Make data for a metadata table."""
    cube = dm.cube.DataCube(_path)

    custom_list = []
    for c in custom:
        custom_list.append(str(np.round(np.max(cube.meta[c].data), 2)))

    dx = cube.meta['dx'].data
    B = (cube.meta['N0'].data * dx)
    Hbf = cube.meta['h0'].data
    qw = (cube.meta['u0'][0].data * Hbf)
    Qw = qw * B
    Qs = Qw * (cube.meta['C0_percent'][0].data / 100)

    return [*custom_list, Qw, Qs, Hbf, B, qw]


def make_scalar_datatable(cube, i, If, time_idx, custom):
    """Make the scalar data in the same format as timeseries data.
    """
    _time = scale_model_time(
        cube.t[time_idx], units='years', If=If)
    _ID = np.ones(shape=_time.shape, dtype=int) * i

    custom_arr = np.zeros((_time.shape[0], len(custom)))
    for c, cvar in enumerate(custom):
        if cvar in cube.meta.keys():
            cval = cube.meta[cvar].data
            if cval.ndim == 2:
                cval = np.max(cval)
        else:
            cval = np.nan
        cvallong = np.ones(shape=_time.shape) * cval
        custom_arr[:, c] = cvallong

    scalar_data = pd.DataFrame(
        data=np.column_stack((_ID, time_idx, cube.t[time_idx],
                              _time, custom_arr)),
        columns=['ID', 'time_idx', 'time_sec', 'time', *custom])
    return scalar_data


def fit_line(x, y):
    finite = np.logical_and(np.isfinite(x), np.isfinite(y))
    x = x[finite]
    y = y[finite]
    x = np.atleast_2d(x).T
    y = np.atleast_2d(y).T
    mdl = LinearRegression().fit(x, y)
    return mdl


def _get_fit_coefficient(_xs, _vals):
    _any_valid_vals = np.any(~np.isnan(_vals))
    if _any_valid_vals:
        _vals_fit = fit_line(_xs, _vals)
        _coef = _vals_fit.coef_
    else:
        _coef = np.nan
    return _coef


def _find_xy_shoreline(_edge):
    # find shoreline xys
    _whr_y, _whr_x = np.where(_edge)
    _whr_land = _whr_y <= 4
    _whr_x = _whr_x[~_whr_land]
    _whr_y = _whr_y[~_whr_land]
    return _whr_x, _whr_y


def compute_channel_visitation_array(cube, channelmask_array):
    visitation_array = (np.sum(channelmask_array, axis=0) /
                        channelmask_array.shape[0])
    return visitation_array


def extract_data_from_mask_arrays(cube, mask_array_set,
                                  time_idx):

    # unpack the mask data
    landmask_array, channelmask_array, shorelinemask_array, _ = mask_array_set  # noqa: E501
    __dx = float(cube.meta.dx)
    __CTR = cube.meta['CTR'].data
    __L0 = cube.meta['L0'].data

    delta_area = np.full(time_idx.shape, np.nan)
    avg_radius = np.full(time_idx.shape, np.nan)
    std_radius = np.full(time_idx.shape, np.nan)

    shore_len = np.full(time_idx.shape, np.nan)
    shore_rough = np.full(time_idx.shape, np.nan)

    # loop through each and do computations and appends
    for i, t in enumerate(time_idx):
        # compute delta area
        _da = np.sum(landmask_array[i, :, :].squeeze()) * __dx * __dx
        delta_area[i] = _da

        # compute edge (quicker than edge mask?)
        _edge = skimage.filters.sobel(landmask_array[i, :, :].squeeze())
        _edge[:5, :] = 0
        _edge = skimage.morphology.thin(_edge)

        # find shoreline location across entire topset
        _whr_x, _whr_y = _find_xy_shoreline(_edge)

        # compute radii and stats
        if _whr_x.size > 0:
            _rad = np.sqrt((_whr_x - __CTR)**2 + (_whr_y - __L0)**2)
            _meanrad = np.mean(_rad * __dx)
            _stdrad = np.std(_rad * __dx)
        else:
            _meanrad = np.nan
            _stdrad = np.nan
        avg_radius[i] = _meanrad
        std_radius[i] = _stdrad

        # compute shoreline length
        _smi = shorelinemask_array[i, :, :]
        if (np.sum(_smi) > 5):  # use 5 for safety
            _sl = dm.plan.compute_shoreline_length(_smi) * __dx
        else:
            _sl = np.nan
        shore_len[i] = _sl

        # compute shoreline roughness
        if _da > 0:
            _sr = _sl / np.sqrt(_da)
        else:
            _sr = np.nan
        shore_rough[i] = _sr

    mask_data = pd.DataFrame(
        data=np.column_stack((delta_area, avg_radius, std_radius,
                              shore_len, shore_rough)),
        columns=['delta_area', 'delta_avg_radius', 'delta_std_radius',
                 'delta_shorelength', 'delta_shorerough'])

    return mask_data


def extract_data_from_mask_arrays_subsided_selenga(cube, mask_array_set,
                                                   time_idx):

    # unpack the mask data
    landmask_array, channelmask_array, shorelinemask_array, subsidencemask = mask_array_set  # noqa: E501

    # unpack other data
    __dx = float(cube.meta.dx)
    __CTR = cube.meta['CTR'].data
    __L0 = cube.meta['L0'].data

    avg_subsd_radius = np.full(time_idx.shape, np.nan)
    avg_subsd_prograte = np.full(time_idx.shape, np.nan)
    channelized_area_in_mask = np.full(time_idx.shape, np.nan)
    avg_subsd_radius_alongpath = np.full(time_idx.shape, np.nan)

    # compute mask area
    mask_area = np.sum(subsidencemask)

    # loop through each and do computations and appends
    for i, t in enumerate(time_idx):

        # compute edge (quicker than edge mask?)
        _edge = skimage.filters.sobel(landmask_array[i, :, :].squeeze())
        _edge[:5, :] = 0
        _edge = skimage.morphology.thin(_edge)

        # find shoreline location in subsd side
        _edge = np.logical_and(_edge, subsidencemask)
        _whr_x, _whr_y = _find_xy_shoreline(_edge)
        if _whr_x.size > 0:
            _rad = np.sqrt((_whr_x - __CTR)**2 + (_whr_y - __L0)**2)
            _meanradsubsd = np.mean(_rad * __dx)
        else:
            _meanradsubsd = np.nan
        avg_subsd_radius[i] = _meanradsubsd

        # compute progradation rate in the subsd side
        if i > 0:
            _xdiff = avg_subsd_radius[i] - avg_subsd_radius[i-1]
            _tdiff = cube.t[time_idx[i]] - cube.t[time_idx[i-1]]
            _rate = _xdiff / _tdiff
        else:
            _rate = np.nan
        avg_subsd_prograte[i] = _rate

        # find the area of the mask that is channel
        _channel_in_mask = (subsidencemask.astype(int) *
                            channelmask_array[i, :, :])
        if mask_area > 0:
            _channelized = (np.sum(_channel_in_mask) / mask_area)
        else:
            _channelized = np.nan
        channelized_area_in_mask[i] = _channelized

    mask_data = pd.DataFrame(
        data=np.column_stack((avg_subsd_radius, avg_subsd_prograte,
                              channelized_area_in_mask)),
        columns=['subsd_avg_radius', 'subsd_avg_prograte',
                 'subsd_area_channelized'])

    return mask_data


def _get_section_eta(_section, t):
    _eta = _section['eta'][t, :]
    _eta[_eta < 0] = np.nan
    return _eta


def _get_section_stage(_section, t):
    if 'stage' in _section.variables:
        _stage = _section['stage'][t, :]
    else:
        _stage = _section['eta'][t, :] + _section['eta'][t, :]
    _stage[_stage < 0.05] = np.nan
    return _stage


def extract_topset_measurements(cube, time_idx):

    __dx = float(cube.meta.dx)

    # make the required section
    # register the topset slope line
    slope_subsd = dm.section.RadialSection(cube, azimuth=30)

    _eta_slopes = np.zeros((time_idx.shape[0]))
    _stage_slopes = np.zeros((time_idx.shape[0]))
    _num_dry = np.zeros((time_idx.shape[0]))
    _subaerial_area = np.zeros((time_idx.shape[0]))

    for i, t in enumerate(time_idx):
        _eta = _get_section_eta(slope_subsd, t)
        _stage = _get_section_stage(slope_subsd, t)

        _xs = slope_subsd.s * cube.meta['dx'].data

        _eta_slopes[i] = _get_fit_coefficient(_xs, _eta)
        _stage_slopes[i] = _get_fit_coefficient(_xs, _stage)

        # compute the number of elevated banks along the cross section
        _dry_bool = (cube['depth'][t, :, :] <= 0.1)
        _num_dry[i] = float(np.sum(_dry_bool))

        # compute the number of subaerial pixels
        _subaerial_bool = (cube['eta'][t, :, :] > cube.meta['H_SL'][t])
        _subaerial_area[i] = float(np.sum(_subaerial_bool)) * __dx * __dx

    topset_data = pd.DataFrame(
        data=np.column_stack((_eta_slopes, _stage_slopes,
                              _num_dry, _subaerial_area)),
        columns=['eta_slope', 'stage_slope',
                 'dry_pixels', 'subaerial_area'])

    return topset_data


def extract_shore_flux_record(cube, mask_array_set, time_idx):

    # unpack the mask data
    _, _, _, subsidencemask = mask_array_set  # noqa: E501

    __dx = float(cube.meta.dx)

    # make the required sections
    flux_total_base = dm.section.CircularSection(
        cube, radius=80)

    x, = np.where(subsidencemask.flatten() != 0)
    flat_trace = np.ravel_multi_index(
        np.flipud(flux_total_base.trace.T), cube.shape[1:])
    _match = np.isin(flat_trace, x)

    # register split sections
    flux_subsd = dm.section.PathSection(
        cube, path=flux_total_base.trace[_match, :])
    flux_total = dm.section.PathSection(
        cube, path=flux_total_base.trace)

    # preallocate lists
    flux_water_total = np.full(time_idx.shape, np.nan)
    flux_water_subsd = np.full(time_idx.shape, np.nan)

    flux_sed_total = np.full(time_idx.shape, np.nan)
    flux_sed_subsd = np.full(time_idx.shape, np.nan)

    # loop through each time and compute
    for i, t in enumerate(time_idx):
        # compute water fluxes
        disch_total = flux_total['discharge'][t, :]
        disch_subsd = flux_subsd['discharge'][t, :]

        flux_water_total[i] = np.trapz(
            disch_total, flux_total._s) * __dx
        flux_water_subsd[i] = np.trapz(
            disch_subsd, flux_subsd._s) * __dx

        # compute sediment fluxes
        sedflux_total = flux_total['sedflux'][t, :]
        sedflux_subsd = flux_subsd['sedflux'][t, :]

        flux_sed_total[i] = np.trapz(
            sedflux_total, flux_total._s) * __dx
        flux_sed_subsd[i] = np.trapz(
            sedflux_subsd, flux_subsd._s) * __dx

    # compute the ratio (with nans for safety)
    _fwt = np.copy(flux_water_total)
    _fst = np.copy(flux_sed_total)
    _fwt[_fwt == 0] = np.nan
    _fst[_fst == 0] = np.nan

    flux_water_ratio = (flux_water_subsd / _fwt)
    flux_sed_ratio = (flux_sed_subsd / _fst)

    flux_data = pd.DataFrame(
        data=np.column_stack((flux_water_total, flux_water_ratio,
                              flux_sed_total, flux_sed_ratio,)),
        columns=['shore_flux_water_total', 'shore_flux_water_ratio',
                 'shore_flux_sed_total', 'shore_flux_sed_ratio'])

    return flux_data


def extract_flux_into_subsided_block(spinup_cube, cube,
                                     mask_array_set, time_idx):

    __dx = float(cube.meta.dx)

    # unpack the mask data
    _, _, _, subsidencemask = mask_array_set
    _mask = np.array(subsidencemask)

    # determine the location to place the section
    if 'subsd_radius1' in cube.meta.keys():
        # if mississippi case, then use this R
        _R = cube.meta['subsd_radius1']
    else:
        _R = 40  # selenga default
    flux_into_base = dm.section.CircularSection(
        cube, radius=_R+2)
    x, = np.where(_mask.astype(int).flatten() != 0)
    flat_trace = np.ravel_multi_index(
        np.flipud(flux_into_base.trace.T), cube.shape[1:])
    _match = np.isin(flat_trace, x)

    # register split sections
    flux_subsd = dm.section.PathSection(
        cube, path=flux_into_base.trace[_match, :])
    flux_total = dm.section.PathSection(
        cube, path=flux_into_base.trace)
    flux_subsd_0 = dm.section.PathSection(
        spinup_cube, path=flux_into_base.trace[_match, :])
    flux_total_0 = dm.section.PathSection(
        spinup_cube, path=flux_into_base.trace)

    # preallocate lists
    sediment_accum = np.full(time_idx.shape, np.nan)
    sediment_diff = np.full(time_idx.shape, np.nan)
    flux_water_subsd = np.full(time_idx.shape, np.nan)
    flux_sed_subsd = np.full(time_idx.shape, np.nan)
    flux_water_frac_subsd = np.full(time_idx.shape, np.nan)
    flux_sed_frac_subsd = np.full(time_idx.shape, np.nan)
    flux_water_ratio_c = np.full(time_idx.shape, np.nan)
    flux_sed_ratio_c = np.full(time_idx.shape, np.nan)
    flux_sed_diff_c = np.full(time_idx.shape, np.nan)
    flux_sed_diff_c_norm = np.full(time_idx.shape, np.nan)

    def _get_flux_section_time(section, var, _t):
        amount = section[var][_t, :]
        return np.trapz(amount, section._s) * __dx

    # compute the value for the spinup (fixed)
    water_spinup_total = _get_flux_section_time(
            flux_total_0, 'discharge', -1)
    water_spinup_subsd = _get_flux_section_time(
            flux_subsd_0, 'discharge', -1)
    water_spinup_frac = water_spinup_subsd / water_spinup_total
    sediment_spinup_total = _get_flux_section_time(
            flux_total_0, 'sedflux', -1)
    sediment_spinup_subsd = _get_flux_section_time(
            flux_subsd_0, 'sedflux', -1)
    sediment_spinup_frac = sediment_spinup_subsd / sediment_spinup_total

    # loop through each time and compute
    for i, t in enumerate(time_idx):

        diff = np.array(cube['eta'][t, :, :] - cube['eta'][0, :, :])
        # diff_block = np.sum(diff[_mask]) * (np.sum(_mask) * __dx * __dx)
        diff_block = np.sum(diff[_mask] * __dx * __dx)

        sediment_accum[i] = diff_block
        if i > 0:
            sediment_diff[i] = sediment_accum[i] - sediment_accum[i-1]

        # do computations for the flux ratio before and after
        water_cube_total = _get_flux_section_time(
            flux_total, 'discharge', t)
        water_cube_subsd = _get_flux_section_time(
            flux_subsd, 'discharge', t)
        water_cube_frac = water_cube_subsd / water_cube_total
        if water_spinup_frac > 0:
            water_ratio_c = water_cube_frac / water_spinup_frac
        else:
            water_ratio_c = np.nan

        sediment_cube_total = _get_flux_section_time(
            flux_total, 'sedflux', t)
        sediment_cube_subsd = _get_flux_section_time(
            flux_subsd, 'sedflux', t)
        if sediment_cube_total > 0:
            sediment_cube_frac = sediment_cube_subsd / sediment_cube_total
        else:
            sediment_cube_frac = np.nan

        # make deviation metrics
        if sediment_spinup_frac > 0:
            sediment_ratio_c = sediment_cube_frac / sediment_spinup_frac
        else:
            sediment_ratio_c = np.nan
        sediment_diff_c = sediment_cube_frac - sediment_spinup_frac
        sediment_diff_c_norm = (sediment_diff_c) / (1 - sediment_spinup_frac)

        # store all values into some preallocated arrays
        flux_water_subsd[i] = water_cube_subsd
        flux_water_frac_subsd[i] = water_cube_frac
        flux_water_ratio_c[i] = water_ratio_c

        flux_sed_subsd[i] = sediment_cube_subsd
        flux_sed_frac_subsd[i] = sediment_cube_frac
        flux_sed_ratio_c[i] = sediment_ratio_c
        flux_sed_diff_c[i] = sediment_diff_c
        flux_sed_diff_c_norm[i] = sediment_diff_c_norm

    spacetime_data = pd.DataFrame(
        data=np.column_stack((sediment_accum, sediment_diff,
                              flux_water_subsd, flux_water_frac_subsd,
                              flux_water_ratio_c, flux_sed_subsd,
                              flux_sed_frac_subsd, flux_sed_ratio_c,
                              flux_sed_diff_c, flux_sed_diff_c_norm)),
        columns=['sediment_accum', 'sediment_diff',
                 'flux_water_subsd', 'flux_water_frac_subsd',
                 'flux_water_ratio_c', 'flux_sed_subsd',
                 'flux_sed_frac_subsd', 'flux_sed_ratio_c',
                 'flux_sed_diff_c', 'flux_sed_diff_c_norm'])

    return spacetime_data


def extract_correlation_with_inlet(cube, **kwargs):
    """Find a correlation with the inlet, over space.

    This could use a different method, but here we simply use the normalized
    disharge as a proxy for the connectedness.
    """

    dx = cube.meta['dx'].data
    B = (cube.meta['N0'].data * dx)
    Hbf = cube.meta['h0'].data
    qw0 = (cube.meta['u0'][0].data * Hbf)
    Qw0 = qw0 * B
    _arr = np.array(cube['discharge'][-1:, :, :].squeeze()) / qw0

    return _arr


def extract_connectivity_subsided_block(cube, time_idx, mask_array_set,
                                        connectivity_array):
    """Find the connectivity with the inlet.
    """
    # correct the nans in the connectivity array to zeros
    connectivity_array[np.isnan(connectivity_array)] = 0

    # unpack the mask data
    _, _, _, subsidencemask = mask_array_set
    subsidencemask = np.array(subsidencemask)

    __dx = float(cube.meta.dx)
    Hbf = cube.meta['h0'].data
    qw0 = (cube.meta['u0'][0].data * Hbf)
    _subsd_area = float(np.sum(subsidencemask)) * __dx * __dx

    # find the connectivity values in the masked area
    _subsd_conn0 = connectivity_array[subsidencemask]
    # grab the max for the table
    _conn_max0 = np.max(_subsd_conn0)

    conn_max = np.full(time_idx.shape, np.nan)
    # loop through each time and compute
    for i, t in enumerate(time_idx):

        # make a connectivity estimate for the current cube at t0
        _arr = np.array(cube['discharge'][t, :, :].squeeze()) / qw0
        _subsd_conn = _arr[subsidencemask]

        conn_max[i] = np.max(_subsd_conn)

    # proces into a table
    subsd_connectivity_max0 = np.ones(shape=time_idx.shape) * _conn_max0
    subsd_connectivity_max = conn_max
    subsd_connectivity_ratio = conn_max / subsd_connectivity_max0
    subsd_area = np.ones(shape=time_idx.shape) * _subsd_area

    conn_data = pd.DataFrame(
        data=np.column_stack((subsd_connectivity_max0, subsd_connectivity_max,
                              subsd_connectivity_ratio, subsd_area)),
        columns=['subsd_connectivity_max0', 'subsd_connectivity_max',
                 'subsd_connectivity_ratio', 'subsd_area'])

    return conn_data
