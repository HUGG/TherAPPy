import numpy as np
import astropy.units as u
from shapely.geometry import LineString


def calculate_closure_temp(time, temp_K, ea, omega, geom,
                           min_gradient=1e-7 * u.K / u.year):
    """
    calculate closure temperatures using Dodson's (1973) equations
    
    based on glide's fortran code:
    https://github.com/cirederf13/glide/blob/master/transient_geotherm/src/closure_temps.f90
    https://github.com/cirederf13/glide/blob/master/transient_geotherm/src/dodson.f90

    """

    R = 8.314 * u.J / (u.K * u.mol)

    
    cooling_K = np.gradient(temp_K, time)
    #cooling_K = cooling.to(u.K / u.year, )

    ind = cooling_K < min_gradient
    cooling_K[ind] = min_gradient

    tau = R * temp_K ** 2 / (ea * cooling_K)
    Tc = ea / (R * np.log(geom * tau * omega))

    return Tc


def calculate_closure_age(time, temp, thermochron_parameters, verbose=False,
                          closure_temperature_error=0.0, subsample_T_history=5):
    """
    first calculate closure temperatures and then find the intersection between the closure temperature vs time curve 
    and the temperature vs time curve
    """

    ea, geom, omega = thermochron_parameters["Ea"], thermochron_parameters["geom"], thermochron_parameters["omega"]

    temp_K = temp.to(u.K, equivalencies=u.temperature())

    if subsample_T_history is not None:
        # take every x temperature history points to avoid spurious changes in dT/dx with small timesteps
        # todo: replace this with a more elegant moving average algorithm
        time = time[subsample_T_history:-subsample_T_history:subsample_T_history]
        temp_K = temp_K[subsample_T_history:-subsample_T_history:subsample_T_history]

    Tc = calculate_closure_temp(time, temp_K, ea, omega, geom)

    if closure_temperature_error is not None:
        Tc += closure_temperature_error

    if verbose is True:
        print('closure temp = ', Tc)

    if Tc.max() < temp_K.min() or Tc.min() > temp_K.max():
        print('warning, cooling temp outside of range of thermochron temps')
        print('range: ', Tc.min(), Tc.max())
        raise ValueError

        #return np.nan, np.nan

    else:
        xy1 = np.vstack([time.to(u.year).value, temp_K.value]).T
        xy2 = np.vstack([time.to(u.year).value, Tc.value]).T
        line1 = LineString(xy1)
        line2 = LineString(xy2)

        int_pt = line1.intersection(line2)
        if int_pt.geom_type == 'MultiPoint':
            xi, yi = int_pt[0].x, int_pt[0].y
        elif int_pt.geom_type == 'LineString':
            if len(int_pt.coords) > 0:
                # xi, yi = int_pt.coords.xy[:, 0], int_pt.coords.xy[:, 1]
                return np.nan, np.nan
            else:
                return np.nan, np.nan
        else:
            xi, yi = int_pt.x, int_pt.y

        age = xi * u.year
        Tc_int = yi * u.K

        return age
    

def model_AFT_age_and_lengths(time, temperature, AFT_parameters, AFT_model_parameters=None, model="Ketcham2007"):

    available_models = ["Ketcham2007"]
    available_kinetic_parameters = ["Clwt", "Dpar"]

    try:
        assert model in available_models
    except AssertionError:
        msg = f"error, {model} model for AFT is not supported, choose one of {available_models}"
        raise AssertionError(msg)
    
    if model == "Ketcham2007":
        import lib.AFT_model_lib as AFT_model_lib

        time_Myr = time.value / 1e6
        kinetic_parameter = AFT_parameters["kinetic_parameter"]

        try:
            assert kinetic_parameter in ["Clwt", "Dpar"]
        except AssertionError:
            msg = f"error, kinetic_parameter {kinetic_parameter} for AFT is not supported, choose one of {available_kinetic_parameters}"
            raise AssertionError(msg)
        
        if kinetic_parameter == "Clwt":
            kinetic_value = AFT_parameters["Clwt"]
        elif kinetic_parameter == "Dpar":
            kinetic_value = AFT_parameters["Dpar"]
        
        model_results = AFT_model_lib.simulate_AFT_annealing(time_Myr, temperature.value, kinetic_value)

    return model_results