import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import mpl_toolkits.axes_grid1 as axtk

import os
import shutil
import sys
import yaml
import warnings
from packaging import version

import pyDeltaRCM
from pyDeltaRCM.shared_tools import sec_in_day, day_in_yr

import craterstats as cst

# load and format paths for craterstats functions
functions_path = os.path.join(cst.PATH, "config/functions.txt")
craterstats_functions = cst.gm.read_textstructure(functions_path)


def crater_production_function_inverse(pf, irange, y, tol=1e-6):
    """Inverse of the crater production function.

    Parameters
    ----------
    pf
        Production function

    irange
        Range of interest
    """
    log_y = np.log10(y)
    # xrange = np.log10(pf.range)
    xrange = np.log10(irange)
    divisions = 9
    count = 0

    while True:
        x0 = np.linspace(
            xrange[1], xrange[0], divisions
        )  # reverse order because np.interp() requires increasing 'x-values'
        y0 = np.log10(pf.evaluate("cumulative", 10.0**x0))
        x = np.interp(log_y, y0, x0)
        y1 = np.log10(pf.evaluate("cumulative", 10.0**x))
        q = np.searchsorted(y0, log_y)

        xrange = x0[[q, q - 1]]
        count += 1
        if abs(y1 - log_y) < tol or count > 99:
            break

    return 10**x


def generate_CSFD_from_production_function(
    time_interval,
    size_interval,
    domain_area,
    poisson_intervals=True,
    rng_seed=None,
):
    """CSFD using Poisson space events.

    Parameters
    ----------
    time_interval
        2-element list of start time and end time in Ga.

    size_interval
        2-element list of smallest and largest diameter craters to generate in km.

    domain_area
        domain area in km2.

    poisson_intervals
        Use Poisson spaced events? (otherwise just expectation interval).

    rng_seed
        Seed for a random number generator used only for generating the CSFD.
        This enables simulations to use the same CSFD, but place craters in
        different locations (i.e., use the same random number generator from pyDeltaRCM).
    """

    # production and chronology functions
    production_function = cst.Productionfn(craterstats_functions, "Mars, Ivanov (2001)")
    chronology_function = cst.Chronologyfn(
        craterstats_functions, "Mars, Hartmann & Neukum (2001)"
    )

    t_interval = time_interval  # time start, finish in Ga
    Area = domain_area  # modelled area in km2
    poisson_intervals = True  #

    list_d = []  # list of diameters
    # list_dt = []  # list of interarrival times

    diameter_range = size_interval  # in km
    if (diameter_range[0] < production_function.range[0]) or (
        diameter_range[1] > production_function.range[1]
    ):
        warnings.warn(
            "Crater range of interest was outside of production functions's defined range. Extrapolating."
        )

    N = production_function.evaluate(
        "cumulative", [1.0, diameter_range[0], diameter_range[1]]
    )  # default a0
    N1_ratio = N[1] / N[0]

    t = t_interval[0]
    i = 0

    # set up the random number generator for crater sizes and interarrival
    if not (rng_seed is None):
        rng = np.random.default_rng(rng_seed)
        _get_random_uniform = lambda: rng.uniform(0, 1)
    else:
        _get_random_uniform = lambda: pyDeltaRCM.shared_tools.get_random_uniform(1)

    while t <= t_interval[1]:
        # generate crater diameter
        u = _get_random_uniform()
        y = u * (N[1] - N[2]) + N[2]
        d = crater_production_function_inverse(production_function, diameter_range, y)

        # time interval
        phi = chronology_function.phi(t)
        lam = phi * N1_ratio * Area
        if poisson_intervals:
            u = _get_random_uniform()
            dt = -np.log(u) / lam  # proper poisson interval
        else:
            dt = 1.0 / lam  # mean interval

        list_d.append(d)
        # list_dt.append(dt)

        t += dt
        i += 1

    return list_d


def generate_CSFD_from_random_uniform(size_interval, number, rng_seed=None):
    if not (rng_seed is None):
        rng = np.random.default_rng(rng_seed)
        _get_random_uniform = lambda: rng.uniform(0, 1)
    else:
        _get_random_uniform = lambda: pyDeltaRCM.shared_tools.get_random_uniform(1)

    low, high = size_interval

    list_d = np.zeros(number)
    for i in range(number):
        u = _get_random_uniform()
        u_d = u * (high - low) + (low)  # uniform resized to diameter range
        list_d[i] = u_d

    return list_d


class CraterModel(pyDeltaRCM.DeltaModel):
    """A model that includes craters.

    This model samples from a distribution to determine possible impact
    events.

    These events change the elevation of the bed instantaneously, and this is
    encoded (how?) into the strata, so we can determine whethher these events
    are preserved or destroyed.
    """

    def __init__(self, input_file, **kwargs):
        # inherit from base model
        super().__init__(input_file, **kwargs)

    def hook_import_files(self):
        """Define the custom YAML parameters."""
        # whether to run vegetation
        self.subclass_parameters["add_craters"] = {"type": "bool", "default": False}

        # sets the generator for crater diameters.
        #    Given as a float to use a single size, or as a string 'sfd' [default]
        #    to use the Howard et al., size-frequency distribution.
        self.subclass_parameters["p_crater_size"] = {
            "type": ["int", "float", "str"],
            "default": "sfd",  # size-frequency dist (i.e., production function)
        }

        # accumulation time, for crater sfd simulations
        #   this controls the total number of crater generated over
        #   the duration of the run.
        self.subclass_parameters["p_crater_accumulation_duration"] = {
            "type": ["int", "float"],
            "default": 1,  # million years of duration
        }

        # how many timesteps the model will be run for.
        #   this is a requirement for this custom model, because we need to space
        #   the craters uniformly over the duration of the `p_crater_accumulation_duration`.
        #   be careful to set this to the same value sent to the preprocessor (if using).
        self.subclass_parameters["p_timesteps"] = {
            "type": ["int"],
            "default": 1,  # million years of duration
        }

        # accumulation number, for crater fixed size simulations
        #   this will allow some scaling between the flood and
        #   year-based crater metrics
        self.subclass_parameters["p_crater_accumulation_number"] = {
            "type": ["int", "float"],
            "default": 200,  # 200 craters of fixed size
        }

        # crater size range generated by the production function
        #   (only applies to sfd simulations)
        self.subclass_parameters["p_sfd_range"] = {
            "type": ["list"],
            "default": [10, 300],  # 5 to 300 m
        }

        # sets the minimum diameter for craters made with the Howard model.
        #    Given as a multiple of dx, to yield the minimum diameter (recommend <=3)
        self.subclass_parameters["p_howard_min_multiple"] = {
            "type": ["int", "float"],
            "default": 1,  # 1 cell
        }

    def hook_init_output_file(self):
        """Add non-standard grids, figures and metadata to be saved."""

        # save the locations of any craters
        self._save_var_list["crater_rim"] = [
            "cratered_rim_sincesave",
            "ID",
            "i4",
            ("time", "x", "y"),
        ]
        self._save_var_list["crater_ejecta"] = [
            "cratered_ejecta_sincesave",
            "ID",
            "i4",
            ("time", "x", "y"),
        ]

    def hook_after_create_domain(self):
        """Add parameters and fields to the model for all cratering operations."""

        self.domain_area = (self.L * self._dx) * (self.W * self.dx)
        self.domain_area_km = self.domain_area / 1e6

        # -- crater size settings --
        # set up the parameters for crater diameters
        if isinstance(self.p_crater_size, str):
            # size frequency distribution params
            if self.p_crater_size == "sfd":  # size frequency distribution (bounded)
                # set time inteval in absolute time billion years
                time_interval = [
                    -3.5,
                    -3.5 + self.p_crater_accumulation_duration / 1000,
                ]

                # set size interval in km
                size_interval = [bound / 1000 for bound in self.p_sfd_range]

                # call the generator (returns a list of crater sizes)
                crater_list = generate_CSFD_from_production_function(
                    time_interval, size_interval, self.domain_area_km
                )
                # crater_list = [c * 1000 for c in crater_list]
                crater_list = (
                    np.array(crater_list) * 1000
                )  # convert to nparray and meters

                # randomly (uniformly) distribute over the duration of the run
                crater_when = [
                    self._get_random_integer_rescaled(0, self.p_timesteps)
                    for _ in range(len(crater_list))
                ]
                crater_when = np.array(crater_when, dtype=int)
            elif (
                self.p_crater_size == "uniform"
            ):  # size frequency distribution (bounded)
                # set size interval in km
                size_interval = [bound / 1000 for bound in self.p_sfd_range]

                crater_list = generate_CSFD_from_random_uniform(
                    size_interval, self.p_crater_accumulation_number
                )
                crater_list = (
                    np.array(crater_list) * 1000
                )  # convert to nparray and meters

                crater_when = [
                    self._get_random_integer_rescaled(0, self.p_timesteps)
                    for _ in range(len(crater_list))
                ]
                crater_when = np.array(crater_when, dtype=int)
            else:
                raise ValueError()

        elif isinstance(self.p_crater_size, (int, float)):
            # raise NotImplementedError()
            # fixed size every crater
            crater_list = [
                self.p_crater_size for _ in range(self.p_crater_accumulation_number)
            ]

            # randomly (uniformly) distribute over the duration of the run
            crater_when = [
                self._get_random_integer_rescaled(0, self.p_timesteps)
                for _ in range(self.p_crater_accumulation_number)
            ]
            crater_when = np.array(crater_when, dtype=int)

        else:
            raise ValueError("Bad value for p_crater_size.")

        # create array with crater size and timing
        self.crater_schedule = dict(diameter=crater_list, time_iteration=crater_when)

        self._howard_Dref = 7 * 1000  # a reference diameter, 7 km, from Howard, 2007
        self._rad_mult_max_dist = 6  # how many times the radius is affected by crater
        self._rad_rim_ll, self._rad_rim_ul = (
            0.9,
            1.41,
        )  # how far (fraction of radius) to apply the rim tagging

        # -- addtl model fields --
        # field to track where craters have occured (inclusive of all sub-parts)
        #    this is used to force the sand_frac var to -1 on data output
        self.cratered = np.zeros_like(self.depth).astype(int)
        self.cratered_sincesave = np.zeros_like(self.depth).astype(int)

        # field to track where the rim cells are located
        #    this stores ID integers on data output
        self.cratered_rim_sincesave = np.zeros_like(self.depth).astype(int)

        # field to track where the ejecta cells are located
        #    this stores ID integers on data output
        self.cratered_ejecta_sincesave = np.zeros_like(self.depth).astype(int)
        self.cratered_count = int(0)  # start with zero, but first fill will be 0+1

        # what value to fill in sandfrac field
        self.crater_sandfrac_value = -1

        # track bed elevation change
        self.eta_change = np.zeros_like(self.depth)

    def hook_solve_water_and_sediment_timestep(self):
        # called at the top of each time iteration
        #   reset the crater tracker
        self.cratered[:] = 0

    def hook_after_route_sediment(self):
        """Apply cratering according to crater_schedule."""
        # determine change in bed elevation on this time iteration
        self.eta_change = self.eta - self.eta0

        # if cratering is on, run the routine
        if self.add_craters:
            # find craters with crater timestep set to current time iteration
            scheduled_crater_idxs = np.where(
                self.crater_schedule["time_iteration"] == self._time_iter
            )[0]
            # crater the surface with each crater
            for i, crater_idx in enumerate(scheduled_crater_idxs):
                self._crater_surface(self.crater_schedule["diameter"][crater_idx])

    def _get_random_integer_rescaled(self, low=0, high=1):
        """
        Return random integer from low to high (exclusive).
        """
        uniform = pyDeltaRCM.shared_tools.get_random_uniform(1)
        integer = int(low + (uniform * (high - low)))
        return integer

    def _crater_surface(self, diameter):
        """Adds a single crater to the surface.

        Parameters
        ---------
        diameter
            Crater diameter to be added, in meters.
        """

        # grab a prob, this maintains the seeds from legacy CraterModel code,
        #   it could be deleted, but a few sensitivity tests would need
        #   to be rewritten, and this seemed easier.
        prob = pyDeltaRCM.shared_tools.get_random_uniform(1)

        # random choice of location, in cell coordinates
        #   NOTE: excludes the inlet channel plus one
        center = (
            self._get_random_integer_rescaled(self.L0 + 1, self.L),
            self._get_random_integer_rescaled(0, self.W),
        )

        # apply the Howard model for craters, including simple and complex
        radius = diameter / 2
        dists = self._dx * np.sqrt(
            (self.y - center[0]) ** 2 + (self.x - center[1]) ** 2
        )
        dists_flat = dists.flatten()

        # determine the parameters to modify the bed elevations
        H1, H2, H2H1, m, n = self._get_crater_parameters(diameter)

        # apply the landscape change inside the crater
        incrater = dists_flat[dists_flat <= radius]
        in_idx = np.where(dists_flat <= radius)[0]
        # equation for inside the crater
        inDeltaH = H2H1 + H1 * (incrater / radius) ** m
        Ii = 0.9  # rim follows slope
        dists[dists == 0] = self._dx
        weights = (dists / radius) ** (-n)
        Er = np.average(self.eta, weights=weights)  # reference elev
        Ei = np.mean(self.eta.flat[in_idx])
        inDeltaE = inDeltaH + (Er - Ei) * (1 - (Ii * (incrater / radius) ** 2))
        self.eta.flat[in_idx] += inDeltaE

        out_logical = np.logical_and(
            dists_flat > radius,
            dists_flat
            < np.maximum(np.sqrt(2) * self._dx, self._rad_mult_max_dist * radius),
        )  # mod so always grabs one cell out!
        outcrater = dists_flat[out_logical]
        out_idx = np.where(out_logical)[0]
        # equation for outside the crater (ejecta!)
        outDeltaH = H2 * ((2 * (outcrater)) / (diameter)) ** -n
        G = np.minimum((1 - Ii), (outDeltaH / H2))
        if diameter <= (self.dx):
            outDeltaH = np.maximum(outDeltaH, 0.15)
        outDeltaE = outDeltaH + G * (Er - Ei)
        self.eta.flat[out_idx] += outDeltaE

        # solve for rim and ejecta locations
        _rad_rim_ul_adj = np.maximum((2) * self._dx, self._rad_rim_ul * radius)
        _rad_rim_ll_off = np.maximum(
            self._dx, _rad_rim_ul_adj - (self._rad_rim_ll * radius)
        )
        _rad_rim_ll_adj = _rad_rim_ul_adj - _rad_rim_ll_off
        rim_logical = np.logical_and(
            dists_flat >= _rad_rim_ll_adj, dists_flat < (_rad_rim_ul_adj)
        )
        ejecta_logical = np.logical_and(
            dists_flat >= (_rad_rim_ul_adj),
            dists_flat <= (self._rad_mult_max_dist * radius),
        )
        rim_elev = np.max(self.eta.flat[in_idx])

        # update the crater info tracker
        #  set the sand_frac updater to true everywhere inside and outside the crater
        self.cratered_count = int(self.cratered_count + 1)
        self.cratered.flat[in_idx] = self.cratered_count
        self.cratered_sincesave.flat[in_idx] = self.cratered_count
        if out_idx.size > 1:
            self.cratered.flat[out_idx] = self.cratered_count
            self.cratered_sincesave.flat[out_idx] = self.cratered_count

        rim_idx = np.where(rim_logical)[0]
        self.cratered_rim_sincesave.flat[rim_idx] = self.cratered_count
        # solve for the ejecta location and update that tracker
        # where not inside or rim, but is out

        ejecta_idx = np.where(ejecta_logical)[0]
        self.cratered_ejecta_sincesave.flat[ejecta_idx] = self.cratered_count

        _msg = (
            f"Applied crater at "
            f"({center[0]*self.dx:}, {center[1]*self.dx:}), "
            f"size {diameter:0.1f} m diameter by "
            f"{np.max(outDeltaE) + np.abs(np.min(inDeltaE)):0.1f} m tall, "
            f"with {np.abs(np.min(inDeltaE)):0.1f} m below ref elev, "
            f"and ref elevation {Er:0.1f}."
        )
        self.log_info(_msg, verbosity=0)

    def hook_after_finalize_timestep(self):
        # hard change to the surface composition to adjust
        #   to where crater has been formed: -1
        self.sand_frac[(self.cratered > 0)] = self.crater_sandfrac_value
        self.active_layer[(self.cratered > 0)] = self.crater_sandfrac_value

    def hook_output_data(self):
        # called just before output_data, change the surface to crater values.
        #   this is imperfect, but should be *close* assuming the
        #   save interval is frequent enough
        if self._save_time_since_data >= self.save_dt:
            # if the conditions to save on this step will be met
            #   (i.e., met in `ouput_data`), replace the sandfrac data with the
            #   output tracking field.
            self.sand_frac[(self.cratered_sincesave > 0)] = self.crater_sandfrac_value

            # replace current timestep data with the info for since saved
            self.cratered[:] = self.cratered_sincesave[:]

    def hook_output_checkpoint(self):
        # called just after output_data, need to reset the crater tracker
        if self._save_time_since_data == 0:
            # if the counter was just reset, then clear the tracker
            #   WARNING: does doing this here break checkpointing??
            self.cratered_sincesave[:] = 0
            self.cratered_rim_sincesave[:] = 0
            self.cratered_ejecta_sincesave[:] = 0

    def _get_crater_parameters(self, diameter):
        if diameter <= self._howard_Dref:
            H1 = 2.54 * diameter**0.67
            H2 = 1.93 * diameter**0.52
            m = 0.73 * diameter**0.113  # value: 2 to 3
        elif diameter > self._howard_Dref:
            H1 = 12.20 * diameter**0.49
            H2 = 0.79 * diameter**0.6
            m = 0.64 * diameter**0.13  # value: 2 to 3
        else:
            raise ValueError("How did we get here?")
        H2H1 = H2 - H1

        # Howard, 2007:
        # "The exponent n is constrained such that volume deposited on the rim
        # equals the volume excavated from the bowl and ranges from a value of
        # about 3 for a 7 km crater to 3.5 for a 250 km crater."
        n = 3

        return H1, H2, H2H1, m, n


if __name__ == "__main__":
    # set up models and run them in parallel
    param_dict = {
        "Length": 6000.0,
        "Width": 12000.0,
        "dx": 20,
        "verbose": 1,
        "N0_meters": 150,
        "L0_meters": 75,
        "Np_water": 1250,  # 1000
        "Np_sed": 250,  # 125
        "u0": 1.5,
        "h0": 6,
        "SLR": 1e-8,
        "f_bedload": 0.5,
        "itermax": 1,
        "save_eta_figs": True,
        "save_sandfrac_figs": False,
        "save_eta_grids": True,
        "save_stage_grids": True,
        "save_velocity_grids": True,
        "save_sandfrac_grids": True,
        "clobber_netcdf": True,
        "save_checkpoint": True,
        "out_dir": "/scratch/crater_csfd_matrix_3",
        "save_dt": 500000,
        "ensemble": 3,
        "add_craters": True,
    }

    ## for size-frequency distribution runs
    param_dict["p_crater_size"] = "sfd"  # "uniform", "sfd"
    param_dict["p_sfd_range"] = [10, 300]
    param_dict["p_timesteps"] = 10000
    param_dict["matrix"] = {
        "p_crater_accumulation_duration": [1, 10, 100]  # 1, 10, 100 Ma
    }

    ## for uniform frequency runs
    # param_dict["p_crater_size"] = "uniform"  # "uniform", "sfd"
    # param_dict["p_sfd_range"] = [10, 300]
    # param_dict["p_timesteps"] = 6000
    # param_dict["p_crater_accumulation_number"] = 250  # number for uniform
    # param_dict["matrix"] = {
    #     "f_bedload": [0.8, 0.5, 0.2]
    # }

    pp = pyDeltaRCM.Preprocessor(
        param_dict, defer_output=True, timesteps=param_dict["p_timesteps"], parallel=6
    )
    pp.run_jobs(DeltaModel=CraterModel)
