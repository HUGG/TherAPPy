"""
TherAPPy
"""

import itertools
import numpy as np
import astropy.units as u

import therappy.thermochron_models as tm
import therappy.default_thermochron_model_parameters as dtp


class thermochron_object():

    """
    Base class for thermochron objects. These objects can store thermochronological data (track lengths, densities, He conentrations etc..), 
    and derived properties such as thermochonometer ages and temperature histories. Each thermochron object represents a single mineral grain. 

    The aim is to construct this library in such a way that several thermochron objects can be combined to a single sample, 
    and several samples can be combined to a section that for instance represents a borehole or a vertical section of outcrops.

    """

    def __init__(self, mineral, thermochron_model=None, thermochron_data=None, thermochron_model_parameters=None, thermochron_age=None, temperature_history=None):

        # check if thermochron system is available
        available_minerals = ["apatite", "zircon", "muscovite", "biotite", "hornblende"]
        try:
            assert mineral in available_minerals
        except AssertionError:
            raise AssertionError(f"error, the mineral {mineral} that you tried to assign to this object is not implemented in the code"
                                 f", choose one of these isntead: {available_minerals}")

        self.mineral = mineral
        self.model = thermochron_model
        self.thermochron_data =  thermochron_data


    def model_thermochron(self, temperature_history, thermochronometer, model=None, thermochron_parameters=None):
        
        """
        Forward model of thermochronological data (He concentration, fission track density, length distributions, etc..)
        using temperature history and thermochronology parameters as an input

        """

        available_models = ["Dodson", "Ketcham2007", "Meesters2002"]
        available_thermochronometers = ["FT", "He", "Ar"]
        if model is None:
            # take the model specified in the class
            model = self.model
        
        try:
            assert thermochronometer in available_thermochronometers
        except AssertionError:
            raise AssertionError(f"error, thermochronometer {thermochronometer} is not in the list of supported models, "
                                 f"choose one of these instead: {available_thermochronometers}")

        try:
            assert model in available_models
        except AssertionError:
            raise AssertionError(f"error, model {model} is not in the list of supported models, choose one of these instead: {available_models}")

        # take the default parameters
        default_thermochron_parameters = dtp.default_thermochron_params(self.mineral, thermochronometer, model)

        final_thermochron_parameters = default_thermochron_parameters.copy()

        # overwrite default parameters if specified
        if thermochron_parameters is not None:
            for k in thermochron_parameters.keys():
                try:
                    assert k in final_thermochron_parameters.keys()
                except:
                    msg = f"parameter {k} is not in the default thermochron parameter set"
                # update parameter:
                final_thermochron_parameters[k] = thermochron_parameters[k]

        assert final_thermochron_parameters is not None

        # store model and thermochron parameters in class
        # not sure if this is the best way to structure this..
        self.model = model
        self.thermochronometer = thermochronometer
        self.thermochron_parameters = thermochron_parameters

        # get time and temperature
        time_bp, temp = temperature_history["time"], temperature_history["temperature"]
        time = time_bp.max() - time_bp

        self.time_bp = time_bp
        self.time = time

        # model thermochron
        if model == "Dodson":
            modelled_thermochron_age = tm.calculate_closure_age(time, temp, final_thermochron_parameters)
            
            modelled_thermochron_age_bp = time.max() - modelled_thermochron_age 
            model_results = {"modelled_thermochron_age_bp": modelled_thermochron_age_bp}

        elif self.mineral == "apatite" and thermochronometer == "FT":
            model_results = tm.model_AFT_age_and_lengths(time, temp, final_thermochron_parameters)

        elif self.mineral == "apatite" and thermochronometer == "He":
            model_results = tm.model_AHe_age(self, time, temp, final_thermochron_parameters)
        
        else:
            raise ValueError("warning, no thermochron model executed, check input")

        return model_results

    def calculate_thermochron_age(self, thermochron_data):

        """
        Calculate thermochronometer age from measured thermochron data (i.e., He concentration, fission track density, etc..)

        todo...
        
        """
        
        pass
            

def calculate_misfit(predicted_ages, measured_ages, measured_stdev, misfit_type=1, num_params=0):
    """
    Calculates misfit value between measured and predicted thermochronometer ages

    type 1 = Braun et al. (2012) equation 8 (Default)
    type 2 = Braun et al. (2012) equation 9
    type 3 = Braun et al. (2012) equation 10

    Braun, J., Van Der Beek, P., Valla, P., Robert, X., Herman, F., Glotzbach, C., Pedersen, V., Perry, C.,
    Simon-Labric, T. and Prigent, C., 2012. Quantifying rates of landscape evolution and tectonic processes
    by thermochronology and numerical modeling of crustal heat transport using PECUBE. Tectonophysics, 524,
    pp.1-28.

    Parameters
    ----------
    predicted_ages : numpy array
        Array of predicted thermochronometer ages.
    measured_ages : numpy arrray
        Array of measured thermochronometer ages.
    measured_stdev : numpy arrray
        Array of standard deviations for the measured thermochronometer ages.
    misfit_type : int, default=1
        Misfit type to calculate. See equations 8-10 of Braun et al. (2012).
    num_params : int, default=0
        Number of model parameters if using misfit type 2.

    Returns
    -------
    misfit : float
        Calculated misfit value.
    """
    # Caclulate general misfit, modify for types 1 and 2 as needed
    misfit = ((predicted_ages - measured_ages) ** 2 / measured_stdev**2).sum()

    if misfit_type == 1:
        misfit = np.sqrt(misfit) / len(predicted_ages)

    if misfit_type == 2:
        misfit = misfit / (len(predicted_ages) - num_params - 1)

    return misfit