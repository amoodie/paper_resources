import numpy as np

import os
import shutil
import sys
import yaml

import pyDeltaRCM
from pyDeltaRCM.shared_tools import sec_in_day, day_in_yr


class SelengaModel(pyDeltaRCM.DeltaModel):
    """SelengaModel.

    This model subclass is used to generate the fault-slip configuration in
    Set 1 of Moodie and Passalacqua.
    """

    def __init__(self, input_file, **kwargs):

        # before loading random state from file, need to get a new random seed
        self._new_seed = np.random.randint((2**32) - 1, dtype='u8')

        # inherit from base model
        super().__init__(input_file, **kwargs)

        # set up custom initial conditions
        self._spin_up = kwargs.pop('spin_up', False)
        self._slipped = False
        self._seed_changed = False

        if not self._spin_up:
            # manually trigger the subsidence field to be generated
            self._toggle_subsidence = True
            self.init_subsidence()

    def init_subsidence(self):
        """Initialize subsidence pattern constrained to a tighter region.

        Uses theta1 and theta2 to set the angular bounds for the
        subsiding region. theta1 and theta2 are set in relation to the
        inlet orientation. The inlet channel is at an angle of 0, if
        theta1 is -pi/3 radians, this means that the angle to the left of
        the inlet that will be included in the subsiding region is 30
        degrees. theta2 defines the right angular bounds for the subsiding
        region in a similar fashion.
        """
        _msg = 'Initializing custom subsidence field'
        self.log_info(_msg, verbosity=1)

        if self._toggle_subsidence:

            _width = (np.pi / 3)

            # right
            _theta1_limits = (-1.5707963267948966, 1.5707963267948966 - _width)
            theta1 = np.random.uniform(_theta1_limits[0], _theta1_limits[1])
            theta2 = theta1 + _width

            R1 = 0.3 * self.L  # radial limits (fractions of L)
            R2 = 0.9 * self.L

            Rloc = np.sqrt((self.y - self.L0)**2 + (self.x - self.W / 2.)**2)

            thetaloc = np.zeros((self.L, self.W))
            thetaloc[self.y > self.L0 - 1] = np.arctan(
                (self.x[self.y > self.L0 - 1] - self.W / 2.)
                / (self.y[self.y > self.L0 - 1] - self.L0 + 1))
            self.subsidence_mask = ((R1 <= Rloc) & (Rloc <= R2) &
                                    (theta1 <= thetaloc) &
                                    (thetaloc <= theta2))
            self.subsidence_mask[:self.L0, :] = False

            # set the sigma field with no time dependence
            self.sigma = self.subsidence_mask * self.subsidence_rate

            _msg1 = 'Theta1 set to: {0}'.format(theta1)
            _msg2 = 'Theta2 set to: {0}'.format(theta2)
            self.log_info(_msg1, verbosity=0)
            self.log_info(_msg2, verbosity=0)

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

                # apply the slip event
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

    # base yaml configuration
    base_output = '/scratch/faults_selenga'
    base_yaml = './yaml/selenga.yaml'
    spin_up_dir = 'spin_up'
    matrix_dir = 'matrix'

    # if running the spinup
    if _input_flag == 'spinup':

        _mdl = SelengaModel(
            input_file=base_yaml,
            out_dir=os.path.join(base_output, spin_up_dir),
            subsidence_rate=0,
            save_checkpoint=True,
            save_dt=864000,
            spin_up=True)

        # solve for how many timesteps
        targ_dur = 600  # target spinup duration
        targ_dur_mdl = (targ_dur*sec_in_day*day_in_yr) * If
        tsteps = int((targ_dur_mdl // _mdl.time_step) + 1)

        try:
            for i in range(tsteps):
                _mdl.update()
        except Exception as e:
            print("ERROR!")
            _mdl.logger.exception(e)

        # finalize
        _mdl.finalize()

    # if running the matrix
    elif _input_flag == 'matrix':

        # open the base yaml file and ammend it as a dict input to pp
        _file = open(base_yaml, mode='r')
        _loader = pyDeltaRCM.shared_tools.custom_yaml_loader()
        param_dict = yaml.load(_file, Loader=_loader)
        _file.close()

        # add a matrix with sigmas to the dict
        _matrix = {'subsidence_rate': [0, 0.01, 0.02, 0.05,
                                       0.1, 0.2, 0.5, 1, 2, 5]}
        param_dict['matrix'] = _matrix
        param_dict['ensemble'] = 7

        # add other configurations
        param_dict.update(
            {'out_dir': os.path.join(base_output, matrix_dir),
             'save_checkpoint': False,
             'resume_checkpoint': True,
             'save_dt': 864000,  # 864000, 10 days
             'parallel': 12})  # 8

        # let the preprocessor write the initial matrix and
        #    create a new output netcdf file
        pp = pyDeltaRCM.Preprocessor(
            param_dict,
            timesteps=2100)  # 2100

        # prepare the additonal information for each folder of the matrix
        for j, j_file in enumerate(pp.file_list):
            # get folder
            j_folder = j_file.parent

            # copy the spinup checkpoint to each of the folders
            shutil.copy(
                src=os.path.join(base_output, spin_up_dir, 'checkpoint.npz'),
                dst=os.path.join(base_output, matrix_dir, j_folder))

        # run the jobs
        pp.run_jobs(DeltaModel=SelengaModel)

    # otherwise
    else:
        raise ValueError('Invalid argument supplied.')
