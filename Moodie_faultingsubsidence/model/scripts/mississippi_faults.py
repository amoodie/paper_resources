import numpy as np

import os
import shutil
import sys
import yaml

import pyDeltaRCM
from pyDeltaRCM.shared_tools import sec_in_day, day_in_yr

from smt.sampling_methods import LHS


class MississippiFaultsModel(pyDeltaRCM.DeltaModel):
    """MississippiModel.

    This model subclass is used to generate the fault-slip configuration in
    Set 2 of Moodie and Passalacqua.
    """

    def __init__(self, input_file, **kwargs):

        # before loading, need to pop relevant keys from the input_file
        self._spin_up = kwargs.pop('spin_up', False)

        # Non-default params are not set to the model, so we have to open
        #   the file and load the dictionary before it is cleared.
        #   NOTE: as of newer pyDeltaRCM versions, there is a much
        #         cleaner way to do this, see docs.
        _file = open(input_file, mode='r')
        _loader = pyDeltaRCM.shared_tools.custom_yaml_loader()
        param_dict = yaml.load(_file, Loader=_loader)
        _file.close()

        # set params to model for later use
        self._width = param_dict.pop('width', None)
        self._length = param_dict.pop('length', None)
        self._depth = param_dict.pop('depth', None)
        self._R_lims = param_dict.pop('R_lims', None)
        self.subsidence_rate = self._depth
        self._width_meters = None
        self._length_meters = None
        self._R = None

        # before loading random state from file, need to get a new random seed
        self._new_seed = np.random.randint((2**32) - 1, dtype='u8')

        # inherit from base model
        super().__init__(input_file, **kwargs)

        # set up custom initial conditions
        self._slipped = False
        self._seed_changed = False

        if not self._spin_up:
            # manually trigger the subsidence field to be generated
            self._toggle_subsidence = True
            self.init_subsidence()

    def hook_init_output_file(self):
        # save dimension variables to output file
        self._save_var_list['meta']['subsd_width_radian'] = ['_width',
                                                             'radians',
                                                             'f4', ()]
        self._save_var_list['meta']['subsd_width'] = ['_width_meters',
                                                      'meters',
                                                      'f4', ()]
        self._save_var_list['meta']['subsd_length'] = ['_length_meters',
                                                       'meters',
                                                       'f4', ()]
        self._save_var_list['meta']['subsd_depth'] = ['_depth',
                                                      'meters',
                                                      'f4', ()]
        self._save_var_list['meta']['subsd_radius1'] = ['_R',
                                                        'pixels',
                                                        'f4', ()]

    def init_subsidence(self):
        """Initialize subsidence pattern constrained to a tighter region.

        In the Set 2 experiments, we vary the size and location of the
        subsided region. Use the parameters grabbed from the YAML file during
        init to set up the sigma array. Note that calculations are
        implemented in other methods.
        """
        _msg = 'Initializing custom subsidence field'
        self.log_info(_msg, verbosity=1)

        if self._toggle_subsidence:

            # pick the radian angle of the block
            self._R, self._theta = self.pick_center_radian_random(
                self._width, self._R_lims)

            # set up the sigma mask array
            self.subsidence_mask, width_m = self.make_subsidence_block(
                self._R, self._theta, self._width, self._length, self._depth)
            self._width_meters = width_m
            self._length_meters = self._length * self.dx

            # set the sigma field with no time dependence
            self.sigma = self.subsidence_mask * self._depth

    def hook_run_one_timestep(self):
        """hook one timestep.

        Hook the `run_one_timestep` method to change the seed when resuming
        from model runs. This method only has an effect on the first timestep
        after resuming from the checkpoint.
        """
        # only do anything on the first timestep
        if not self._seed_changed:

            # get the new seed
            _ns = int(self._new_seed)

            # change the state of the random number generator
            pyDeltaRCM.shared_tools.set_random_seed(_ns)
            self._seed_changed = True

            # always write the seed to file for record and reproducability
            _msg = 'Random seed was changed to: %s ' % str(_ns)
            self.log_info(_msg, verbosity=0)

    def apply_subsidence(self):
        """Overwrite this method to handle our spin-up and slip cases.

        Here, we only take an action if the run is not a spin-up case, and has
        not already "slipped". That is, this method only has an effect on the
        first timestep after resuming from the checkpoint.
        """
        if not self._spin_up:

            if not self._slipped:

                self.apply_earthquake_slip_event()
                self._slipped = True

                # save the initial state
                self._save_time_since_data = float("inf")
                self.log_model_time()
                self.output_data()

    def apply_earthquake_slip_event(self):
        """Apply the earthquake slip.

        Custom method to apply the earthquake slip. Here, we simply adjust eta
        by `sigma` (calcualted in `init_subsidence`). Then we apply light
        smoothing of the bed topography by 10 iterations of `topo_diffusion`
        and smoothing of the water surface by 5 iterations of the water
        routing routines.
        """

        _msg = 'Applying fault slip event: {slp} m'.format(
            slp=self.subsidence_rate)
        self.log_info(_msg, verbosity=1)

        # adjust eta by sigma
        self.eta[:] = self.eta - self.sigma

        # apply iterations of topo diffusion
        for _ in range(10):
            self.topo_diffusion()

        # rebalance the stage
        self.stage[:] = np.maximum(self.stage, self.H_SL)
        Hsmth = pyDeltaRCM.water_tools._smooth_free_surface(
            np.copy(self.stage), self.cell_type, self._Nsmooth, self._Csmooth)
        self.stage = (((1 - self._omega_sfc) * self.stage) +
                      (self._omega_sfc * Hsmth))

        # apply a flooding correction
        self.flooding_correction()

        # update depth
        self.depth[:] = np.maximum(0, self.stage - self.eta)

        # adjust for the abrupt change in stage.
        for iteration in range(5):
            self.init_water_iteration()
            self.run_water_iteration()
            self.compute_free_surface()
            self.finalize_water_iteration(iteration)

    def pick_center_radian_random(self, _W, _R_lims):
        # is constrained random pick, based on the width of the block to place

        # pick the radius
        _R = np.random.randint(_R_lims[0], _R_lims[1], size=1)

        # then solve for width at that distance in radians
        _W_astheta = _W / _R
        # _W_astheta = (_W)  # / _R
        _half_W_astheta = _W_astheta / 2

        _abs_theta_min = -1.5707963267948966
        _abs_theta_max = 1.5707963267948966
        _theta_min = _abs_theta_min + _half_W_astheta
        _theta_max = _abs_theta_max - _half_W_astheta

        _theta = np.random.uniform(_theta_min, _theta_max)

        return _R, _theta

    def make_subsidence_block(self, _R, _theta, _W, _L, _D):

        _W_astheta = (_W / _R)
        _half_W_astheta = _W_astheta / 2

        # theta limits
        theta1 = _theta - _half_W_astheta
        theta2 = _theta + _half_W_astheta

        # radial limits
        R1 = _R
        R2 = _R + _L

        # compute the masked area
        Rloc = np.sqrt((self.y - self.L0)**2 + (self.x - self.W / 2.)**2)

        thetaloc = np.zeros((self.L, self.W))
        thetaloc[self.y > self.L0 - 1] = np.arctan(
            (self.x[self.y > self.L0 - 1] - self.W / 2.)
            / (self.y[self.y > self.L0 - 1] - self.L0 + 1))
        subsidence_mask = ((R1 <= Rloc) & (Rloc <= R2) &
                           (theta1 <= thetaloc) &
                           (thetaloc <= theta2))
        subsidence_mask[:self.L0, :] = False

        # convert the close corners to points to compute the width in m
        left_x, left_y = R1*np.cos(theta1), R1*np.sin(theta1)
        right_x, right_y = R1*np.cos(theta2), R1*np.sin(theta2)

        width_m = np.sqrt((left_x-right_x)**2+(left_y-right_y)**2) * self._dx

        return subsidence_mask, width_m


# Set up the script to run the model when called by name
if __name__ == '__main__':

    # get script input argument
    #   input should be `--spinup` or `--matrix`.
    _arg = sys.argv
    if len(_arg) == 3:
        _input_flag = sys.argv[1].strip('-')
        _force = sys.argv[2].strip('-')
    elif len(_arg) == 2:
        _input_flag = sys.argv[1].strip('-')
        _force = False
    else:
        raise ValueError('No arguments supplied.')

    # parameter choices for scaling
    If = 7 / 365.25  # intermittency factor for year-scaling
    SLR_mmyr = 2
    SLR = pyDeltaRCM.preprocessor.scale_relative_sea_level_rise_rate(
        SLR_mmyr, If=If)
    print("SLR is:", SLR)

    # base yaml configuration
    base_output = '/scratch/faults_mississippi'
    base_yaml = './yaml/mississippi.yaml'
    spin_up_dir = 'spin_up/set_00'
    matrix_dir = 'matrix_03'

    # if running the spinup
    if _input_flag == 'spinup':

        _mdl = MississippiFaultsModel(
            input_file=base_yaml,
            out_dir=os.path.join(base_output, spin_up_dir),
            subsidence_rate=0,
            SLR=SLR,
            save_checkpoint=True,
            save_dt=8640000,
            spin_up=True)

        # solve for how many timesteps
        targ_dur = 2000  # target spinup duration
        targ_dur_mdl = (targ_dur*sec_in_day*day_in_yr) * If
        tsteps = int((targ_dur_mdl // _mdl.time_step) + 1)
        print("expected tsteps:", tsteps)

        try:
            for i in range(tsteps):
                _mdl.update()
        except Exception as e:
            print("ERROR!")
            _mdl.logger.exception(e)

        # finalize
        _mdl.finalize()

    elif _input_flag == 'tempresume':

        _mdl = MississippiFaultsModel(
            input_file=base_yaml,
            out_dir=os.path.join(base_output, spin_up_dir),
            subsidence_rate=0,
            SLR=SLR,
            save_checkpoint=True,
            resume_checkpoint=True,
            save_dt=1,
            spin_up=True)

        # run two timesteps
        _mdl.update()
        _mdl.update()

        # finalize
        _mdl.finalize()

    # if running the matrix
    elif _input_flag == 'matrix':

        # open the base yaml file and ammend it as a dict input to pp
        _file = open(base_yaml, mode='r')
        _loader = pyDeltaRCM.shared_tools.custom_yaml_loader()
        param_dict = yaml.load(_file, Loader=_loader)
        _file.close()

        # make the set list for simulations
        set_list = []
        D_set = [0.1, 1, 5]
        W_set = np.array([3600, 12000, 24000]) / 1200
        L_set = 9600 / 1200
        for d in D_set:
            for w in W_set:
                _d = {'width': w,
                      'length': L_set,
                      'depth': d}
                set_list.append(_d)
        param_dict['set'] = set_list
        param_dict['ensemble'] = 11
        param_dict['R_lims'] = [40, 75]  # limiting distance from apex (cells)

        # add other configurations
        param_dict.update(
            {'out_dir': os.path.join(base_output, matrix_dir),
             'SLR': SLR,
             'save_checkpoint': False,
             'resume_checkpoint': True,
             'save_dt': 1728000,  # 1728000
             'parallel': 14})  # how many to run in parallel

        # let the preprocessor write the initial matrix and
        #    create a new output netcdf file
        pp = pyDeltaRCM.Preprocessor(
            param_dict,
            timesteps=5000)

        # prepare the additonal information for each folder of the matrix
        for j, j_file in enumerate(pp.file_list):
            # get folder
            j_folder = j_file.parent

            # copy the spinup checkpoint to each of the folders
            shutil.copy(
                src=os.path.join(base_output, spin_up_dir, 'checkpoint.npz'),
                dst=os.path.join(base_output, matrix_dir, j_folder))

        # run the jobs
        pp.run_jobs(DeltaModel=MississippiFaultsModel)

    # otherwise
    else:
        raise ValueError('Invalid argument supplied.')
