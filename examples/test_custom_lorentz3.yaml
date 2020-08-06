import numpy as np
from numpy import inf, power

name = "lorentzian3"
title = "This model calculates an empirical functional form for SAS data \
characterized by three Lorentzian-type functions."
description = """I(q) = scale_1/(1.0 + pow((q*length_1),exponent_1))
             + scale_2/(1.0 + pow((q*length_2),exponent_2) )
             + scale_3/(1.0 + pow((q*length_3),exponent_3) )+ background
             scale_1    = Lorentzian term scaling #1
             length_1   = Lorentzian screening length #1 [A]
             exponent_1 = Lorentzian exponent #1
             scale_2    = Lorentzian term scaling #2
             length_2   = Lorentzian screening length #2 [A]
             exponent_2 = Lorentzian exponent #2
             scale_3    = Lorentzian term scaling #3
             length_3   = Lorentzian screening length #3 [A]
             exponent_3 = Lorentzian exponent #3
             background = Incoherent background
        """
category = "shape-independent"

# pylint: disable=bad-whitespace, line-too-long
#            ["name", "units", default, [lower, upper], "type", "description"],
parameters = [["lorentz_scale_1",  "",     100.0, [-inf, inf], "", "First power law scale factor"],
              ["lorentz_length_1", "Ang", 1000.0, [-inf, inf], "", "First Lorentzian screening length"],
              ["lorentz_exp_1",    "",      3.0, [-inf, inf], "", "First exponent of power law"],
              ["lorentz_scale_2",  "",      10.0, [-inf, inf], "", "Second scale factor for broad Lorentzian peak"],
              ["lorentz_length_2", "Ang",  100.0, [-inf, inf], "", "Second Lorentzian screening length"],
              ["lorentz_exp_2",    "",      2.0, [-inf, inf], "", "Second exponent of power law"],
              ["lorentz_scale_3",  "",      1.0, [-inf, inf], "", "Third scale factor for broad Lorentzian peak"],
              ["lorentz_length_3", "Ang",  10.0, [-inf, inf], "", "Third Lorentzian screening length"],
              ["lorentz_exp_3",    "",      1.0, [-inf, inf], "", "Third exponent of power law"],
             ]
# pylint: enable=bad-whitespace, line-too-long


def Iq(q,
       lorentz_scale_1=100.0,
       lorentz_length_1=1000.0,
       lorentz_exp_1=3.0,
       lorentz_scale_2=10.0,
       lorentz_length_2=100.0,
       lorentz_exp_2=2.0)
       lorentz_scale_2=1.0,
       lorentz_length_2=10.0,
       lorentz_exp_2=1.0):

    """
    :param q:                   Input q-value (float or [float, float])
    :param lorentz_scale_1:     Second scale factor for broad Lorentzian peak
    :param lorentz_length_1:    First Lorentzian screening length
    :param lorentz_exp_1:       Exponent of the first Lorentz function
    :param lorentz_scale_2:     Second scale factor for broad Lorentzian peak
    :param lorentz_length_2:    Second Lorentzian screening length
    :param lorentz_exp_2:       Exponent of the second Lorentz function
    :param lorentz_scale_3:     Third scale factor for broad Lorentzian peak
    :param lorentz_length_3:    Third Lorentzian screening length
    :param lorentz_exp_3:       Exponent of the Third Lorentz function
    :return:                    Calculated intensity
    """
# pylint: disable=bad-whitespace
    intensity  = lorentz_scale_1/(1.0 +
                                  power(q*lorentz_length_1, lorentz_exp_1))
    intensity += lorentz_scale_2/(1.0 +
                                  power(q*lorentz_length_2, lorentz_exp_2))
    intensity += lorentz_scale_3/(1.0 +
                                  power(q*lorentz_length_3, lorentz_exp_3))
# pylint: enable=bad-whitespace
    return intensity

Iq.vectorized = True  # Iq accepts an array of q values

def random():
    """Return a random parameter set for the model."""
    scale = 10**np.random.uniform(0, 4, 3)
    length = 10**np.random.uniform(1, 4, 3)
    expon = np.random.uniform(1, 6, 3)

    pars = dict(
        #background=0,
        scale=1, # scale provided in model
        lorentz_scale_1=scale[0],
        lorentz_length_1=length[0],
        lorentz_exp_1=expon[0],
        lorentz_scale_2=scale[1],
        lorentz_length_2=length[1],
        lorentz_exp_2=expon[1],
        lorentz_scale_3=scale[2],
        lorentz_length_3=length[2],
        lorentz_exp_3=expon[2],
    )
    return pars

