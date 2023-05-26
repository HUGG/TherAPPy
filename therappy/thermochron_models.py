import pdb
import numpy as np
import astropy.units as u
from astropy.units import cds
#from shapely.geometry import LineString
from astropy.units import imperial


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


def find_intersections(A, B):
    """
    Find intersection between two lines

    source for this function: https://stackoverflow.com/questions/3252194/numpy-and-line-intersections

    Parameters
    ----------
    A : numpy array
        two column array containing the x and y coordinates of the first line
    B : _type_
        two-column array containing the x and y coordinates of the second line

    Returns
    -------
    xi : numpy array
        array with x-coordinates of intersection points
    yi : numpy array
        array with y-coorindates of intersection points
    """
    
    # min, max and all for arrays
    amin = lambda x1, x2: np.where(x1<x2, x1, x2)
    amax = lambda x1, x2: np.where(x1>x2, x1, x2)
    aall = lambda abools: np.dstack(abools).all(axis=2)
    slope = lambda line: (lambda d: d[:,1]/d[:,0])(np.diff(line, axis=0))

    x11, x21 = np.meshgrid(A[:-1, 0], B[:-1, 0])
    x12, x22 = np.meshgrid(A[1:, 0], B[1:, 0])
    y11, y21 = np.meshgrid(A[:-1, 1], B[:-1, 1])
    y12, y22 = np.meshgrid(A[1:, 1], B[1:, 1])

    m1, m2 = np.meshgrid(slope(A), slope(B))
    m1inv, m2inv = 1/m1, 1/m2

    yi = (m1*(x21-x11-m2inv*y21) + y11)/(1 - m1*m2inv)
    xi = (yi - y21)*m2inv + x21

    xconds = (amin(x11, x12) < xi, xi <= amax(x11, x12), 
              amin(x21, x22) < xi, xi <= amax(x21, x22) )
    yconds = (amin(y11, y12) < yi, yi <= amax(y11, y12),
              amin(y21, y22) < yi, yi <= amax(y21, y22) )

    return xi[aall(xconds)], yi[aall(yconds)]


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
        # line1 = LineString(xy1)
        # line2 = LineString(xy2)

        # int_pt = line1.intersection(line2)
        # if int_pt.geom_type == 'MultiPoint':
        #     xi_check, yi_check = int_pt[0].x, int_pt[0].y
        # elif int_pt.geom_type == 'LineString':
        #     if len(int_pt.coords) > 0:
        #         # xi, yi = int_pt.coords.xy[:, 0], int_pt.coords.xy[:, 1]
        #         return np.nan, np.nan
        #     else:
        #         return np.nan, np.nan
        # else:
        #     xi_check, yi_check = int_pt.x, int_pt.y

        xi_, yi_ = find_intersections(xy1, xy2)
        xi, yi = xi_[0], yi_[0]

        # try:
        #     assert np.allclose(xi, xi_check)
        # except AssertionError:
        #     pdb.set_trace()

        age = xi * u.year
        Tc_int = yi * u.K

        return age
    

def model_AFT_age_and_lengths(time, temperature, model_parameters, model="Ketcham2007"):

    available_models = ["Ketcham2007"]
    available_kinetic_parameters = ["Clwt", "Dpar"]

    try:
        assert model in available_models
    except AssertionError:
        msg = f"error, {model} model for AFT is not supported, choose one of {available_models}"
        raise AssertionError(msg)
    
    if model == "Ketcham2007":
        import therappy.AFT_model_lib as AFT_model_lib

        # take kinetic value (Clwt, Dpar or rmr0) out of dictionary to make it a separate argument for the AFT model
        kinetic_par = model_parameters["kinetic_parameter"]
        kinetic_value = model_parameters[kinetic_par]
        final_model_parameters = model_parameters.copy()
        del final_model_parameters[kinetic_par]
        
        model_results = AFT_model_lib.simulate_AFT_annealing(time.value, temperature.value, kinetic_value,  **final_model_parameters)

    return model_results


def model_AHe_age(mineral, t, T, thermochron_parameters, use_fortran_algorithm=True, n_eigenmodes=15):
                #   D0_div_a2=np.exp(13.4),
                #   Ea=32.9 * 4184,
                #   R=8.3144621,
                #   decay_constant_238U=4.916e-18,
                #   decay_constant_232Th=1.57e-18,
                #   decay_constant_235U=3.12e-17,
                #   alpha_ejection=True,
                #   stopping_distance=20e-6,
                #   method='RDAAM',
                #   alpha=0.04672,
                #   C0=0.39528,
                #   C1=0.01073,
                #   C2=-65.12969,
                #   C3=-7.91715,
                  
                #   ):

    """
    """

    import therappy.AHe_models as AHe_models

    method = thermochron_parameters["diffusivity_model"]
    available_methods = ["Farley2000", "RDAAM", "Wolf1996"]
    try:
        assert method in available_methods
    except AssertionError:
        msg = f"the diffusivity model {method} is not implemented. Choose one of the following methods instead: {available_methods}"
        raise AssertionError(msg)
    
    radius_, U_, Th_ = mineral.radius, mineral.U, mineral.Th

    temperature_K = T.to(u.K, equivalencies=u.temperature())

    t_sec = t.to(u.s)

    radius = radius_.to(u.m).value #* 1e-6
    print(f"radius = {radius}")
    
    # convert U from ppm to kg/kg
    U = U_.to(cds.ppm).value * 1e-6
    Th = Th_.to(cds.ppm).value * 1e-6

    print(f"U = {U}")
    print(f"U = {Th}")

    # calculate He production:
    decay_constant_238U = thermochron_parameters["decay_constant_238U"]
    decay_constant_232Th = thermochron_parameters["decay_constant_232Th"]   
    decay_constant_235U = thermochron_parameters["decay_constant_235U"]

    U238 = (137.88 / 138.88) * U
    U235 = (1.0 / 138.88) * U
    Th232 = Th
    Ur0 = 8 * U238 * decay_constant_238U + 7 * U235 * decay_constant_235U \
          + 6 * Th232 * decay_constant_232Th
    decay_constant = Ur0 / (8*U238 + 7*U235 + 6*Th232)

    if method is 'Farley2000':

        #D0 = D0_div_a2 * radius ** 2
        # values in HeFTy 1.8.3:
        D0 = 50.0 / 1e4     # m2/sec
        Ea = 32.9 * 4184.0  # J/mol
        D_div_a2 = D0 / (radius**2) * np.exp(-Ea / (R*temperature_K.value))
        D = D_div_a2 * radius**2
        #print 'using Farley (2000) diffusion parameters'

    elif method is 'RDAAM':
        #print 'using RDAAM model to calculate helium diffusivity'
        #print 'with U238=%0.3e, U235=%0.3e, Th232=%0.3e, radius=%0.3e' % \
        #      (U238, U235, Th232, radius)

        alpha = thermochron_parameters["alpha"]
        C0, C1, C2, C3 = thermochron_parameters["C0"], thermochron_parameters["C1"], thermochron_parameters["C2"], thermochron_parameters["C3"]

        D = AHe_models.calculate_He_diffusivity_RDAAM(temperature_K.value, t_sec.value, U238, U235, Th232, radius,
                                                      alpha=alpha, C0=C0, C1=C1,
                                                      C2=C2, C3=C3,
                                                      use_fortran_algorithm=use_fortran_algorithm)

    elif method is 'Wolf1996':
        # diffusivity params Wolf et al (1996), table 7, Durango
        # tested, values are really given in log10 instead of ln
        # big difference with D0 values in Flowers (2009), not sure why
        R=8.3144621
        #log_D0_div_a2 = 7.7 #(1/sec)
        log_D0_div_a2 = 7.82 #(1/sec) , value given in HeFTy 1.8.3

        D0_div_a2 = 10**log_D0_div_a2
        #D0_div_a2 = np.exp(log_D0_div_a2)
        Ea = 36.3 * 4184
        D0 = D0_div_a2 * radius ** 2
        D = (D0 / radius**2 * np.exp(-Ea / (R*temperature_K.value))) * radius**2

    else:
        msg = 'error, cannot determine method for calculating helium ' \
              'diffusivity, choose "Wolf1996", "Farley2000", ' \
              'or "RDAAM", current method = %s' % method
        raise ValueError(msg)

    alpha_ejection = thermochron_parameters["model_alpha_ejection"]
    stopping_distance = thermochron_parameters["stopping_distance"]

    ahe_age = AHe_models.He_diffusion_Meesters_and_Dunai_2002(
        t_sec.value, D, radius, Ur0,
        decay_constant=decay_constant,
        U_function='exponential',
        n_eigenmodes=n_eigenmodes,
        alpha_ejection=alpha_ejection,
        stopping_distance=stopping_distance)
    
    ahe_age_yr = (ahe_age * u.s).to(u.year)
    model_results = {"modelled_thermochron_age": ahe_age_yr[-1], "modelled_thermochron_ages": ahe_age_yr,
                     "He_diffusivity": D, "He_diffusion_model": "Meesters2002", "He_diffusivity_model": method}

    return model_results


def model_vitrinite_reflectance(times_forward, temperatures_in, method='easyRo', default_vr_method="easyRo"):
    
    """
    calculate vitrinite reflectance using the easy%Ro algorithm of 
    Burnham & Sweeney (1989) or the basin%Ro model by Nielsen et al. (2015)
    
    this function has been verified by reproducing the %Ro values 
    provided in Figure 8 of Sweeney & Burnham (1990)
    
    References:
    Burnham, A. K., & Sweeney, J. J. (1989). A chemical kinetic model of vitrinite maturation and reflectance. Geochimica et Cosmochimica Acta, 53(10), 2649–2657.
    Nielsen, S. B., Clausen, O. R., & McGregor, E. (2015). basin%Ro: A vitrinite reflectance model derived from basin and laboratory data. Basin Research, 29, 515–536. https://doi.org/10.1111/bre.12160

    Parameters:
    -----------
    times_bp : array
        time in years bp. The last item should be 0 years bp
    temperatures : array
       temperature at each timestep, in degrees C
    vr_method : string      
        VR method, choose either 'easyRo' or 'basinRo'. default is 'easyRo'
        
    Returns:
    --------
    Ro : array           
        calculated %Ro values at each timestep
    """

    vr_methods = ["easyRo", "basinRo"]


    if method == None:
        method = default_vr_method

    try:
        assert method in vr_methods
    except AssertionError:
        msg = f"error, vr_method {method} not recognized. Choose one of these instead: {vr_methods}"
        raise AssertionError(msg)

    #Myr = 1e6 * 365.25 * 24 * 60 * 60
    Myr = (1e6 * u.year).to(u.s)
    
    timesteps = len(times_forward)

    # convert T to Kelvin
    temperatures_K = temperatures_in.to(u.K, equivalencies=u.temperature())

    # fixed parameters:
    if method is 'easyRo':
        weights = np.array([0.03, 0.03, 0.04, 0.04, 0.05, 0.05, 0.06,
                            0.04, 0.04, 0.07, 0.06, 0.06, 0.06, 0.05,
                            0.05, 0.04, 0.03, 0.02, 0.02, 0.01 ])
    elif method is 'basinRo':
        weights = np.array([0.0185,	0.0143,	0.0569,	0.0478,	0.0497,	0.0344,	0.0344,
                            0.0322,	0.0282,	0.0062,	0.1155,	0.1041,	0.1023,	0.076,
                            0.0593,	0.0512,	0.0477,	0.0086,	0.0246,	0.0096])

    # activation energy (cal/mol)
    activationEnergy = np.array([34., 36., 38., 40., 42., 44., 46., 48., 50.,
                                 52., 54., 56., 58., 60., 62., 64., 66., 68.,
                                 70., 72. ]) * 1000.0 * imperial.cal / u.mol
    
    components = len(activationEnergy)
    
    # preexponential factor (s^-1)
    if method is 'easyRo':
        preExp = 1.0e13 / u.s
    elif method is 'basinRo':
        preExp = np.exp(60.9856) / Myr
    
    # universal gas constant (cal K-1 mol-1)  
    R = 1.987 * imperial.cal / (u.K * u.mol)
    
    # rearrange time & activation energy arrays
    T_all_components = np.resize(temperatures_K.value,(components, timesteps)) * u.K
    activationEnergies_all_timesteps = np.resize(activationEnergy,(timesteps, components)).T
    
    # ERT
    EdivRT = activationEnergies_all_timesteps / (R*temperatures_K)
    
    # I
    I = preExp * T_all_components * np.exp(-EdivRT) * (1-(EdivRT**2+2.334733*EdivRT+0.250621) / 
                                        (EdivRT**2+3.330657*EdivRT+1.681534) )
    
    # calculate heating rates (K / sec)
    heatingRates = np.diff(temperatures_K) / np.diff(times_forward.to(u.s))
    
    # zero heating rate at first timestep
    heatingRates = np.insert(heatingRates, [0], 0)
    
    # remove zero heating rates
    heatingRates[heatingRates.value==0] = 1.0e-30 * u.K / u.s
    
    # delta I
    deltaI = np.zeros(I.shape)
    for j in range(1, timesteps, 1):
        deltaI[:, j] = deltaI[:, j-1] + (I[:, j]-I[:, j-1]) / heatingRates[j]
    
    #cumulativeReacted
    cumulativeReacted = weights * (1.0 - np.exp(-deltaI.T))
    for k in range(components):
        cumulativeReacted[deltaI[k, :] > 220.0, k] = weights[k]
    cumulativeReacted[deltaI.T < 1.0e-20] = 0
    
    #sumReacted
    sumReacted = cumulativeReacted.sum(axis=1)
    
    # calculate Ro at each timestep
    if method is 'easyRo':
        Ro = np.exp(-1.6 + 3.7 * sumReacted)
    elif method is 'basinRo':
        Rzero = 0.2104
        Ro = Rzero * np.exp(3.7 * sumReacted)
    
    model_results = {"modelled_Ro": Ro[-1], "modelled_Ro_all_timesteps": Ro, "vr_method": method}

    return model_results