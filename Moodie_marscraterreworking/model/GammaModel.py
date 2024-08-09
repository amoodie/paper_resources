import pyDeltaRCM


class GammaModel(pyDeltaRCM.DeltaModel):
    """
    Model class implemented to evaluate sensitivity to the gamma parameter.
    
    Gamma is calculated using the gravitational acceleration constant. Note, 
    simulations using this model class are not a good representation of delta
    building on Mars, as there are additional processes that would need to be
    modulated by change in gravity and are not parameterized as such in the
    pyDeltaRCM model as of date of publication. 
    
    Simulations using this model class are described in the paper's 
    Supplementary Information, Section S1.
    """
    def __init__(self, input_file=None, **kwargs):
        super().__init__(input_file, **kwargs)

    def hook_import_files(self):
        """Define custom gravity parameter."""
        # custom numeric parameter
        self.subclass_parameters["gravity"] = {"type": ["float"], "default": 9.81}

    def hook_create_domain(self):
        """Called before gamma is set"""
        self.g = self.gravity
        self.gamma = (
            self.g * self.S0 * self._dx / (self.u0**2)
        )  # water weighting coeff

        _msg = f"Modifying parameters for sensitivity test. g = {self.g} and gamma = {self.gamma}."
        self.log_info(_msg, verbosity=0)


if __name__ == "__main__":
    # get script input argument

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
        "itermax": 1,
        "save_eta_figs": True,
        "save_eta_grids": True,
        "clobber_netcdf": True,
        "save_checkpoint": False,
        "out_dir": "/scratch/gamma_sensitivity",
        "save_dt": 500000,
    }

    param_dict["matrix"] = {
        "gravity": [9.81, 3.71],
        "f_bedload": [0.2, 0.5, 0.8],
    }
    param_dict["timesteps"] = 10000
    param_dict["parallel"] = 6

    pp = pyDeltaRCM.Preprocessor(param_dict)

    pp.run_jobs(DeltaModel=GammaModel)
