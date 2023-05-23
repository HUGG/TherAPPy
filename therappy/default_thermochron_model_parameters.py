

import astropy.units as u


def default_thermochron_params(mineral, thermochronometer, thermochron_model):
    
    """
    Default parameters for the Dodson cooling age model are based on 
    Reiners and Brandon (2006, https://doi.org/10.1146/annurev.earth.34.031405.125202)
    
    """
    
    assert thermochron_model in ["Dodson"]

    if thermochron_model == "Dodson":
        # default grain size for calculating omega for the MAr thermochronometer
        a_Mar = 100.0 * 1e-6 * u.m

        if thermochronometer == "FT":
        # these values are taken from reiners' review paper and correspond to Ea (x1000.) and omega (x secinyr)
            if mineral == 'apatite':
                # here are Ea and D0/a2 values for AFT from ketcham. 1999
                # taken from reiners 2004
                ea = 147 * 1e3 * u.J / u.mol
                geom = 1
                omega = 2.05e6 / u.s
            # if thermochron_system == 'AFT_min':
            #     # here are Ea and D0/a2 values for AFT from ketcham. 1999
            #     # taken from reiners 2004
            #     ea = 138 * 1e3 * u.J / u.mol
            #     geom = 1
            #     omega = 5.08e5 / u.s
            # if thermochron_system == 'AFT_max':
            #     # here are Ea and D0/a2 values for AFT from ketcham. 1999
            #     # taken from reiners 2004
            #     ea = 187 * 1e3 * u.J / u.mol
            #     geom = 1
            #     omega = 1.57e8 / u.s
            elif mineral == 'zircon':
                # here are Ea and D0/a2 values for ZFT from reiners2004
                # taken from reiners 2004
                ea = 208 * 1e3 * u.J / u.mol  # !/(4.184*1.e3)
                geom = 1
                omega = 4.0e8 / u.s
                # energy=224.d3
                # geom=1.d0
                # diff=1.24e8*3600.*24.*365.25e6
                # diff = 1.24e8*3600.*24.*365.25e6
        elif thermochronometer == 'He':
            if mineral == 'apatite':
                # here are Ea and D0/a2 values for AHe from Farley et al. 2000
                # taken from reiners 2004
                ea = 138 * 1e3 * u.J / u.mol
                geom = 1
                omega = 7.64e7 / u.s
            elif mineral == 'zircon':
                # here are Ea and D0/a2 values for ZHe from reiners2004
                # taken from reiners 2004
                ea = 169 * 1e3 * u.J / u.mol  # !/(4.184*1.e3)
                geom = 1
                omega = 7.03e5 / u.s
                # energy=178.d3
                # geom=1.d0
                # diff=7.03d5*3600.d0*24.d0*365.25d6

        # the following are for argon argon, might be a bit much for most, 51,52,53 in glide
        elif thermochronometer == 'Ar':
            if mineral == "hornblende":
                # here are Ea and D0/a2 values for hbl from harrison81
                # taken from reiners 2004
                ea = 268 * 1e3 * u.J / u.mol
                geom = 1
                omega = 1320 / u.s
            # elif thermochron_system == 'MAr_old':
            #     # here are Ea and D0/a2 values for mus from hames&bowring1994,robbins72
            #     # taken from reiners 2004
            #     ea = 180 * 1e3 * u.J / u.mol
            #     geom = 1
            #     omega = 3.91 / u.s

            elif mineral == "muscovite":
                # Values from Harrison et al. (2009, https://doi.org/10.1016/j.gca.2008.09.038)
                ea_cal = 63 * 1e3 * u.imperial.cal / u.mol
                ea = ea_cal.to(u.J / u.mol)
                geom = 1
                D0 = 2.3e-4 * u.m**2 / u.s

                omega = 55 * D0 / a_Mar**2

            # elif thermochron_system == 'MAr_min':
            #     # Values from Harrison et al. (2009, https://doi.org/10.1016/j.gca.2008.09.038)
            #     ea_cal = 56 * 1e3 * u.imperial.cal / u.mol
            #     ea = ea_cal.to(u.J / u.mol)
            #     geom = 1
            #     D0 = 72.3 * 1e-4 * u.m**2 / u.s
            #     omega = 55 * D0 / a_Mar**2

            # elif thermochron_system == 'MAr_max':
            #     # Values from Harrison et al. (2009, https://doi.org/10.1016/j.gca.2008.09.038)
            #     ea_cal = 70 * 1e3 * u.imperial.cal / u.mol
            #     ea = ea_cal.to(u.J / u.mol)
            #     geom = 1
            #     D0 = 0.1 * 1e-4 * u.m**2 / u.s
            #     omega = 55 * D0 / a_Mar**2

            elif mineral == 'biotite':
                # here are Ea and D0/a2 values for bio from grove&harrison1996
                # taken from reiners 2004
                ea = 197 * 1e3 * u.J / u.mol
                geom = 1
                omega = 733. / u.s

        thermochron_parameters = {"Ea": ea, "geom": geom, "omega": omega}

    return thermochron_parameters

    