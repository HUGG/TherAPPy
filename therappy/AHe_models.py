__author__ = 'elco'

import numpy as np

from . import AFT_model_lib as AFT

# import fortran module
try:
    #from . import AFT_model_core
    from . import AFT_model_core
except ImportError:
    print('-' * 30)
    print('warning: failed to import fortran annealing module')
    print('use slower python implementation of AFT annealing module instead')
    print('compile the fortran module by running the following command ' \
          'in the source directory of this module:')
    print('f2py -c AFT_model_core.f90 ' \
          '-m AFT_model_core')
    print('-' * 30)


def He_diffusion_Meesters_and_Dunai_2002(t, D, radius, Ur0,
                                         U_function='constant',
                                         shape='sphere',
                                         decay_constant=1e-8,
                                         n_eigenmodes=50,
                                         x=0,
                                         all_timesteps=True,
                                         alpha_ejection=True,
                                         stopping_distance=20e-6):

    """
    Helium diffusion model by Meesters and Dunai (2002)

    t : numpy array
        time (s)
    D : numpy array
        diffusivity for each timestep (m2 s-1)
    radius : float
        radius of diffusion domain (m)
    Ur0 : ?

    U_function : string, optional
        'constant' or 'exponential'
    shape : string, optional
        shape of the modeled diffusion domain, default is 'sphere'
    decay_constant: float
        decay constant
    n_eigenmodes : int
        number of eigenmodes to evaluate
    x : float
        initial concentration of daughter isotope
    all_timesteps : bool
        report AHe age for all timesteps or just the last
    alpha_ejection : bool
        model alpha ejection or not
    stopping_distance : float
        alpha ejection stopping distance (m)

    """

    # find shape params
    n = np.arange(1, n_eigenmodes+1, 1)
    if shape == 'sphere':
        # eq 13:
        mu = (n * np.pi / radius) ** 2
        gamma = 6.0 / ((n * np.pi) ** 2)
    elif shape == 'finite cylinder':
        print('need to figure out bessel functions in numpy')
        # http://docs.scipy.org/doc/scipy-0.13.0/reference/special.html
        raise AssertionError(f"error, shape {shape} has not been implemented yet in the AHe module")
    elif shape == 'infinite cylinder':
        print('need to figure out bessel functions in numpy')
        raise AssertionError(f"error, shape {shape} has not been implemented yet in the AHe module")

    # experimental, implement alpha ejection algorithm
    # this is basically adjusting the gamma term, see eq. 24 in Meesters
    # and Dunai (2002) pt. 2
    # sigma is alpha ejection distance. acc. Farley et al. (1996) GCA, 
    # varies between 10-30 um for apatites
    if alpha_ejection is True:
        a = radius
        k = n * np.pi / a
        sigma = stopping_distance
        gamma_old = gamma
        gamma = 3.0 / (n * np.pi)**2 * (1.0 - sigma / (2 * a) 
                       + 1 / (n * np.pi) * 1 / (k * sigma) 
                       * (1.0 - np.cos(k * sigma)) 
                       + 1.0 / (k * sigma) * np.sin(k * sigma))

    nt = len(t)

    # calculate decay time
    decay_time = 1.0 / decay_constant

    # calculate F function (eq ..)
    if U_function == 'constant':
        F = t
    elif U_function == 'exponential':
        F = decay_time * (1.0 - np.exp(-t / decay_time))
    else:
        msg = 'please supply value for U_function, '
        msg += 'choose either "constant" or "exponential"'
        raise ValueError(msg)

    # eq. 5, time integration of diffusivity:
    xi = np.zeros(nt)
    for j in range(0, nt-1):
        xi[j + 1] = xi[j] + (D[j] + D[j+1]) / 2.0 * (t[j+1] - t[j])

    # eq. 6
    Fa = (F[1:] - F[:-1]) / (xi[1:] - xi[:-1])

    # iterate over all eigenmodes:
    beta = np.zeros((nt, nt))

    if all_timesteps is True:

        cn = np.zeros((nt, n_eigenmodes))

        for N in range(nt):

            for n in range(n_eigenmodes):

                # eq. 8
                #beta[:J] = np.exp(-mu[n] * (xi[-1] - xi[:J]))
                beta[N, :] = np.exp(-mu[n] * (xi[N] - xi))

                # right hand side of eq. 7:
                #beta_sum = 0
                #for j in range(0, N):
                #    beta_sum += (beta[N, j+1] - beta[N, j]) * Fa[j]
                beta_sum = np.sum((beta[N, 1:N+1] - beta[N, :N]) * Fa[:N])

                # eq. 7
                cn[N, n] = x + Ur0 * gamma[n] / mu[n] * beta_sum

        # eq. 1
        Cav = cn.sum(axis=1)

    else:
        cn = np.zeros(n_eigenmodes)

        for n in range(n_eigenmodes):

            # eq. 8
            beta[n, :] = np.exp(-mu[n] * (xi[-1] - xi[:]))

            # right hand side of eq. 7:
            beta_sum = np.sum((beta[n, 1:] - beta[n, :-1]) * Fa)

            # eq. 7
            cn[n] = x + Ur0 * gamma[n] / mu[n] * beta_sum

        # eq. 1
        Cav = cn.sum()

    # find age, not sure if Ur0 is the correct production term here...
    t_c = Cav / Ur0

    return t_c


def calculate_He_diffusivity_RDAAM(temperature_input, time, U238, U235, Th232, radius,
                                   use_fortran_algorithm=False,
                                   kinetic_parameter='Clwt',
                                   kinetic_value=0.0,
                                   rmr0_min=0,
                                   rmr0_max=0.85,
                                   alpha=0.04672,
                                   C0=0.39528,
                                   C1=0.01073,
                                   C2=-65.12969,
                                   C3=-7.91715):

    """
    calculate He diffusivity as a function of radiation damage
    acc. to RDAAM model, FLowers et al. (2009) GCA

    """

    # convert temperature from degr. C to K
    #temperature = temperature_input + 273.15
    temperature = temperature_input

    try:
        a = type(AFT_model_core)
    except NameError:
        use_fortran_algorithm = False

    log_omega_p = -22.0
    log_phi_p = -13.0
    Etrap = 34.0 * 1000.0 # J/mol
    ln_D0_L_div_a2 = 9.733
    E_L = 122.3 * 1000.0 # J/mol
    L = 8.1e-4  # in cm (!)
    R = 8.3144621

    # calculate reduced track density
        # get annealing kinetics:
    if kinetic_parameter != 'rmr0':
        rmr0, kappa = \
            AFT.calculate_kinetic_parameters(kinetic_parameter, kinetic_value)
    else:
        print('using rmr0 as kinetic parameter')
        rmr0 = kinetic_value
        kappa = 1.04 - rmr0

    if np.isnan(rmr0) is True or rmr0 <= rmr0_min:
        print('!! warning, rmr0 lower than minimum')
        print('!! %s = %0.3f' % (kinetic_parameter, kinetic_value))
        print('!! setting rmr0 to %0.3f' % rmr0_min)
        rmr0 = rmr0_min
        kappa = 1.04 - rmr0
    elif rmr0 > rmr0_max:
        print('!! warning, rmr0 value exceeds most resistant apatite in ' \
              'Carlson (1999) dataset')
        print('!! adjusting rmr0 from %0.3f to %0.3f' % (rmr0, rmr0_max))
        rmr0 = rmr0_max
        kappa = 1.04 - rmr0

    # calculate length of timesteps in sec
    dts = (time[1:] - time[:-1])
    temperature_midpoint = (temperature[1:] + temperature[:-1]) / 2.0
    nsteps = len(dts)

    if use_fortran_algorithm is True:

        # fortran module for reduced track lengths:
        # call fortran module to calculate reduced fission track lengths
        #rmf, rcf = AFT_model_core.reduced_ln(dts, temperature_midpoint, rmr0, kappa,
        #                                                    alpha, C0, C1, C2, C3, nsteps)
        rcf = AFT_model_core.reduced_ln(
            dts, temperature_midpoint, rmr0, kappa, alpha, C0, C1, C2, C3, nsteps)
        rmf = AFT.caxis_project_reduced_lengths(rcf)
        #rmf, rcf = AFT_model_core.reduced_ln(
        #    dts, temperature, rmr0, kappa, alpha, C0, C1, C2, C3, nsteps)

        # correct 0 length tracks:
        rmf[rmf < 0] = 0.0

        rm = rmf
        rc = rcf

    else:
        print('use python reduced track ln function:')
        # warning, reduced length is not correct
        # check against fortran function
        # python reduced track length function:
        r_cmod = AFT.calculate_reduced_track_lengths(dts, temperature, 
                                                     C0=C0, C1=C1, C2=C2, C3=C3,
                                                     alpha=alpha)
        rcp = AFT.kinetic_modifier_reduced_lengths(r_cmod, rmr0, kappa)
        rmp = AFT.caxis_project_reduced_lengths(rcp)
        rm = rmp
        rc = rcp

    rho_r = AFT.calculate_normalized_density(rc)

    decay_const_U238 = 4.916e-18
    decay_const_Th232 = 1.57e-18
    decay_const_U235 = 3.12e-17

    Na = 6.022e23 # avagadro number

    # density of apatite
    density = 3190.0

    atomic_mass_U238 = 238.05078826
    atomic_mass_U235 = 235.0439299
    atomic_mass_Th232 = 232.0377

    U238_g = U238 * 1000
    U238_mol = U238_g / atomic_mass_U238
    U238_atoms_per_kg = U238_mol * Na
    U238_atoms_per_cm3 = U238_atoms_per_kg * density / (100**3)

    U235_g = U235 * 1000
    U235_mol = U235_g / atomic_mass_U235
    U235_atoms_per_kg = U235_mol * Na
    U235_atoms_per_cm3 = U235_atoms_per_kg * density / (100**3)

    Th232_g = Th232 * 1000
    Th232_mol = Th232_g / atomic_mass_Th232
    Th232_atoms_per_kg = Th232_mol * Na
    Th232_atoms_per_cm3 = Th232_atoms_per_kg * density / (100**3)

    # calculate volume density of tracks
    # unit rho_v = tracks/cm3
    t1 = time[:-1]
    t2 = time[1:]
    rho_v = (8.0/8.0 * U238_atoms_per_cm3 * (np.exp(decay_const_U238 * t2) - np.exp(decay_const_U238 * t1))
             + 7.0/8.0 * U235_atoms_per_cm3 * (np.exp(decay_const_U235 * t2) - np.exp(decay_const_U235 * t1))
             + 6.0/8.0 * Th232_atoms_per_cm3 * (np.exp(decay_const_Th232 * t2) - np.exp(decay_const_Th232 * t1)))

    # calculate rho_s
    eta_q = 0.91
    lambda_f_yr = 8.46e-17  # fission track decay constant (yr-1)
    #year = 365.25 * 24 * 60.0 * 60.0
    year = 365 * 24 * 60.0 * 60.0
    lambda_f = lambda_f_yr / year  # fission track decay constant (s-1)
    decay_const = 1.551e-4
    lambda_D = decay_const / year / 1e6

    e_rho_s = lambda_f / lambda_D * rho_v * eta_q * L * rho_r

    # cumulative sum of radiation damage, e_rho_s:
    e_rho_s_sum = np.cumsum(e_rho_s)

    # calculate diffusivity
    C = 10**log_phi_p * e_rho_s_sum + 10**log_omega_p * e_rho_s_sum**3
    D_div_a2 = np.exp(ln_D0_L_div_a2) * np.exp(-E_L / (R * temperature_midpoint)) / (C * np.exp(Etrap / (R*temperature_midpoint)) + 1)
    D = D_div_a2 * radius**2

    # convert diffusivity from midpoint to full array
    D_final = np.zeros_like(temperature)
    D_final[0] = D[0]
    D_final[-1] = D[-1]
    D_final[1:-1] = (D[1:] + D[:-1]) / 2.0

    return D_final



