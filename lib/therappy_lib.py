"""
TherAPPy
"""

import itertools
import numpy as np
import astropy.units as u

import lib.thermochron_models as tm
import lib.default_thermochron_model_parameters as dtp


class thermochron_object():

    """
    Base class for thermochron objects. These objects can store thermochronological data (track lengths, densities, He conentrations etc..), 
    and derived properties such as thermochonometer ages and temperature histories. Each thermochron object represents a single mineral grain. 

    The aim is to construct this library in such a way that several thermochron objects can be combined to a single sample, 
    and several samples can be combined to a section that for instance represents a borehole or a vertical section of outcrops.

    """

    def __init__(self, thermochron_system, thermochron_model=None, thermochron_data=None, thermochron_model_parameters=None, thermochron_age=None, temperature_history=None):

        # check if thermochron system is available
        try:
            assert thermochron_system in ["AHe", "AFT", "ZHe", "ZFT"]
        except AssertionError:
            raise AssertionError(f"error, the thermochron system {thermochron_system} that you tried to assign to this object is not implemented in the code")

        self.system = thermochron_system
        self.model = thermochron_model
        self.thermochron_data =  thermochron_data


    def model_thermochron(self, temperature_history, model=None, thermochron_parameters=None):
        
        """
        Forward model of thermochronological data (He concentration, fission track density, length distributions, etc..)
        using temperature history and thermochronology parameters as an input

        """

        available_models = ["Dodson"]

        if model is None:
            # take the model specified in the class
            model = self.model

        assert model in available_models

        if thermochron_parameters is None:
            # take the default parameters
            thermochron_parameters = dtp.default_thermochron_params(self.system, model)

        assert thermochron_parameters is not None

        # store model and thermochron parameters in class
        # not sure if this is the best way to structure this..
        self.model = model
        self.thermochron_parameters = thermochron_parameters

        # get time and temperature
        time, temp = temperature_history["time"], temperature_history["temperature"]

        # model thermochron
        if model == "Dodson":
            self.modelled_thermochron_age = tm.calculate_closure_age(time, temp, thermochron_parameters)
            
    def calculate_thermochron_age(self, thermochron_data):

        """
        Calculate thermochronometer age from measured thermochron data (i.e., He concentration, fission track density, etc..)

        todo...
        
        """
        
        pass
            

