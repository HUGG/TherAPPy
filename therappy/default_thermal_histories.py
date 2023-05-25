"""collection of thermal histories 
"""

import numpy as np
import astropy.units as u


def thermal_history(thermal_history_name, dt=1e4*u.year):

    thermal_history_collection = ["simple", "Wolf1998_1", "Wolf1998_2", 
                                  "Wolf1998_3", "Wolf1998_4", "Wolf1998_5"]

    assert thermal_history_name in thermal_history_collection

    dtf = dt.to(u.year).value

    if thermal_history_name == "simple":
        time = np.arange(0, 1e8, dtf)[::-1] * u.year
        temperature =  np.linspace(300, 10, len(time))  * u.deg_C

    elif thermal_history_name == "Wolf1998_1":
        time = np.arange(0, 1e8, dtf)[::-1] * u.year
        temperature = np.zeros(len(time)) * u.deg_C

        sect1 = time.value > 4e7
        sect2 = time.value <= 4e7
        temperature[sect1] = 135.0 * u.deg_C 
        temperature[sect2] = 15.0 * u.deg_C
        
    elif thermal_history_name == "Wolf1998_2":
        time = np.arange(0, 1e8, dtf)[::-1] * u.year
        temperature =  np.linspace(135, 15, len(time))  * u.deg_C

    elif thermal_history_name == "Wolf1998_3":
        time = np.arange(0, 1e8, dtf)[::-1] * u.year
        temperature = np.zeros(len(time)) * u.deg_C

        temperature[time>20e6 * u.year] = 60.0 * u.deg_C 
        temperature[time<=20e6 * u.year] = 15.0 * u.deg_C

    elif thermal_history_name == "Wolf1998_4":
        time = np.arange(0, 1e8, dtf)[::-1] * u.year
        temperature = np.zeros(len(time)) * u.deg_C

        sect1 = time>75e6 * u.year
        start1 = np.where(sect1)[0][0]
        end1 = np.where(sect1)[0][-1] + 1
        nn1 = end1 - start1

        sect3 = time<=25e6 * u.year
        start3 = np.where(sect3)[0][0]
        end3 = np.where(sect3)[0][-1] + 1
        nn3 = end3 - start3
        
        temperature[:] = 60.0 * u.deg_C
        temperature[sect1] = np.linspace(100, 60.0, nn1) * u.deg_C
        temperature[sect3] = np.linspace(60, 15.0, nn3) * u.deg_C

    elif thermal_history_name == "Wolf1998_5":
        time = np.arange(0, 1e8, dtf)[::-1] * u.year
        temperature = np.zeros(len(time)) * u.deg_C

        sect1 = time > 5e6 * u.year
        start1 = np.where(sect1)[0][0]
        end1 = np.where(sect1)[0][-1] + 1
        nn1 = end1 - start1

        sect2 = time <= 5e6 * u.year
        start2 = np.where(sect2)[0][0]
        end2 = np.where(sect2)[0][-1] + 1
        nn2 = end2 - start2
        
        temperature[:] = 60.0 * u.deg_C
        temperature[sect1] = np.linspace(15, 65.0, nn1) * u.deg_C
        temperature[sect2] = np.linspace(65.0, 15.0, nn2) * u.deg_C

    thermal_history = {"time": time, "temperature": temperature, 
                       "thermal_history_name": thermal_history_name}

    return thermal_history 
