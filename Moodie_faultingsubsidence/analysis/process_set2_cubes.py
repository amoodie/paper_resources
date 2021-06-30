import os
import sys
import tempfile
import gc

import pandas as pd

from multiprocessing import Pool, freeze_support

import deltametrics as dm

from _analysis import process, utilities
from tqdm import tqdm


if __name__ == '__main__':

    freeze_support()

    _tmpdir = tempfile.TemporaryDirectory().name
    # dask.config.set({'temporary-directory': _tmpdir})

    # process the input argument
    if len(sys.argv) == 2:
        _input_flag = sys.argv[1].strip('-')
    else:
        raise ValueError('No arguments supplied.')

    # determine the run configuration
    MODEL_SET = os.path.join(os.path.expandvars('$SCRATCH'),
                             'faults_set2')

    # determine scaling parameters and time slices
    If = (7 / 365.25)

    # set up the spinup cubes
    spinups_path = os.path.join(MODEL_SET, 'spin_up')
    spinups_fldrs = [f.path for f in os.scandir(spinups_path) if f.is_dir()]
    spinups_fldrs.sort()

    # set up the matrix cube list
    matrix_path = os.path.join(MODEL_SET, 'matrix')
    matrix_fldrs = [f.path for f in os.scandir(matrix_path) if f.is_dir()]
    matrix_fldrs.sort()

    SPINUPS_SAMPLE = dm.cube.DataCube(
        os.path.join(spinups_fldrs[0], 'pyDeltaRCM_output.nc'))
    CUBES_SAMPLE = dm.cube.DataCube(
        os.path.join(matrix_fldrs[-1], 'pyDeltaRCM_output.nc'))

    # determine the length of objcts to run
    if _input_flag == 'debug':
        # here, we limit the operations to run more quickly
        _nspinups = min(2, len(spinups_fldrs))
        spinups_fldrs = spinups_fldrs[:_nspinups]
        _ntsspinups = SPINUPS_SAMPLE.shape[0] // 30

        _nmatrix = 20
        _intmatrix = len(matrix_fldrs) // _nmatrix
        _whchmatrix = range(0, len(matrix_fldrs), _intmatrix)
        matrix_fldrs = matrix_fldrs[::_intmatrix]
        _ntsmatrix = CUBES_SAMPLE.shape[0] // 30

        _nsets = ((_nspinups * _ntsspinups) +
                  (_nmatrix * _ntsmatrix))

    else:
        _nspinups = len(spinups_fldrs)
        _nmatrix = len(matrix_fldrs)

        _ntsspinups = SPINUPS_SAMPLE.shape[0]
        _ntsmatrix = CUBES_SAMPLE.shape[0] // 3

        _nsets = ((_nspinups * _ntsspinups) +
                  (_nmatrix * _ntsmatrix))

    print(
        'Running in {mode} mode:\n'
        '  n spinups:    {ns}\n'
        '  n matrix:     {nm}\n'
        '  n sets:       {nsets}\n'.format(
            mode=_input_flag, ns=_nspinups,
            nm=_nmatrix, nsets=_nsets)
    )

    # determine the number of workers and threads to use
    if _input_flag == 'stampede':
        n_workers = 32
    else:
        n_workers = 32

    # determine the time indexes to compute data for
    spinup_time_idx = utilities.define_time_idx_set(
        SPINUPS_SAMPLE, start=0, num=_ntsspinups)

    matrix_time_idx = utilities.define_time_idx_set(
        CUBES_SAMPLE, start=0, num=_ntsmatrix)

    # determine the metadata variables to pull
    custom = ['sigma', 'subsd_width', 'subsd_length',
              'subsd_depth', 'subsd_radius1']

    def process_spinup_job(_path, i, time_idx, pool):

        # create the cube as a delayed operation
        cube = dm.cube.DataCube(_path)

        # make the computation of the all of the OAPs on the cluster
        _OAP_list = process.make_OAP_list(cube, time_idx, pool)

        mask_array_set = process.make_mask_arrays(
            cube, _OAP_list, time_idx)

        delta_mask_data = process.extract_data_from_mask_arrays(
            cube, mask_array_set, time_idx)
        subsd_mask_data = process.extract_data_from_mask_arrays_subsided_selenga(
            cube, mask_array_set, time_idx)

        topset_data = process.extract_topset_measurements(
            cube, time_idx)

        correlation_arr = process.extract_correlation_with_inlet(
            cube, method='discharge')

        scalar_data = process.make_scalar_datatable(
            cube, i, If, time_idx, custom)

        results = pd.concat(
            [scalar_data, delta_mask_data,
             subsd_mask_data, topset_data], axis=1)

        return results, correlation_arr

    def process_matrix_job(_path, i, time_idx, pool, connectivity_array):

        # create the cube as a delayed operation
        cube = dm.cube.DataCube(_path)

        # make the computation of the all of the OAPs on the cluster
        _OAP_list = process.make_OAP_list(cube, time_idx, pool)

        mask_array_set = process.make_mask_arrays(
            cube, _OAP_list, time_idx)

        delta_mask_data = process.extract_data_from_mask_arrays(
            cube, mask_array_set, time_idx)
        subsd_mask_data = process.extract_data_from_mask_arrays_subsided_selenga(
            cube, mask_array_set, time_idx)

        spacetime_data = process.extract_flux_into_subsided_block(
            SPINUPS_SAMPLE, cube, mask_array_set, time_idx)
        connectivity_data = process.extract_connectivity_subsided_block(
                cube, time_idx, mask_array_set, connectivity_array)

        topset_data = process.extract_topset_measurements(
            cube, time_idx)

        scalar_data = process.make_scalar_datatable(
            cube, i, If, time_idx, custom)

        results = pd.concat(
            [scalar_data, delta_mask_data, subsd_mask_data,
             topset_data, spacetime_data, connectivity_data],
            axis=1)

        return results

    # set up the method for processing
    with Pool(n_workers) as pool:

        # process the spinup cubes
        print('\nProcessing spin ups:')
        spinup_table_results = []
        spinup_array_results = []
        for i, spinups_fldr in tqdm(enumerate(spinups_fldrs),
                                    total=len(spinups_fldrs)):
            _path = os.path.join(spinups_fldr, 'pyDeltaRCM_output.nc')
            data_table, data_arrays = process_spinup_job(
                _path, i, spinup_time_idx, pool=pool)
            spinup_table_results.append(data_table)
            spinup_array_results.append(data_arrays)

        # spinup_results = dask.compute(*spinup_table_results)
        spinup_results = spinup_table_results

        # accumulate the tables into one long results table
        spinups_long = pd.concat(spinup_results)
        spinups_long.to_csv(
            os.path.join('output_tables', 'set2_spinups_table.csv'),
            index=True)

        spinup_results = None
        gc.collect()

        # process the matrix cubes
        print('\nProcessing matrix:')
        matrix_table_results = []
        for i, matrix_fldr in tqdm(enumerate(matrix_fldrs),
                                   total=len(matrix_fldrs)):
            _path = os.path.join(matrix_fldr, 'pyDeltaRCM_output.nc')
            data_table = process_matrix_job(
                _path, i, matrix_time_idx, pool=pool,
                connectivity_array=spinup_array_results[0])
            matrix_table_results.append(data_table)

        matrix_results = matrix_table_results

        # accumulate the tables into one long results table
        matrix = pd.concat(matrix_results)
        matrix.to_csv(
            os.path.join('output_tables', 'set2_matrix_table.csv'),
            index=True)

        matrix_results = None
        gc.collect()

    # process the matrix cubes into a single metadata table
    #    this could be folded into the return from the above process?
    summ_tbl = pd.DataFrame(columns=[*custom, 'Qw', 'Qs', 'Hbf', 'B', 'qw'])
    for i, matrix_fldr in enumerate(matrix_fldrs):
        _path = os.path.join(matrix_fldr, 'pyDeltaRCM_output.nc')
        _meta = process.extract_single_metadata(_path, custom)
        summ_tbl.loc[matrix_fldr] = _meta

    # accumulated result written to file
    summ_tbl.to_csv(
        os.path.join('output_tables', 'set2_metasummary_table.csv'),
        index=True)
