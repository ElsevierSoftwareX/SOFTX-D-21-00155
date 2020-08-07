import numpy as np
from numpy import inf, power

name = "loren3"
title = "This model calculates an empirical functional form for SAS data \
characterized by three Lorentzian-type functions."
description = """I(q) = scale_1/(1.0 + pow((q*length_1),exp_1))
             + scale_2/(1.0 + pow((q*length_2),exp_2) )
             + scale_3/(1.0 + pow((q*length_3),exp_3) )+ background
             scale_1    = Lorentzian term scaling #1
             length_1   = Lorentzian screening length #1 [A]
             exp_1 = Lorentzian exponent #1
             scale_2    = Lorentzian term scaling #2
             length_2   = Lorentzian screening length #2 [A]
             exp_2 = Lorentzian exponent #2
             scale_3    = Lorentzian term scaling #3
             length_3   = Lorentzian screening length #3 [A]
             exp_3 = Lorentzian exponent #3
             background = Incoherent background
        """
category = "shape-independent"

# pylint: disable=bad-whitespace, line-too-long
#            ["name", "units", default, [lower, upper], "type", "description"],
parameters = [["scale_1",  "",     100.0, [-inf, inf], "", "First power law scale factor"],
              ["length_1", "Ang", 1000.0, [-inf, inf], "", "First Lorentzian screening length"],
              ["exp_1",    "",      3.0, [-inf, inf], "", "First exponent of power law"],
              ["scale_2",  "",      10.0, [-inf, inf], "", "Second scale factor for broad Lorentzian peak"],
              ["length_2", "Ang",  100.0, [-inf, inf], "", "Second Lorentzian screening length"],
              ["exp_2",    "",      2.0, [-inf, inf], "", "Second exponent of power law"],
              ["scale_3",  "",      1.0, [-inf, inf], "", "Third scale factor for broad Lorentzian peak"],
              ["length_3", "Ang",  10.0, [-inf, inf], "", "Third Lorentzian screening length"],
              ["exp_3",    "",      1.0, [-inf, inf], "", "Third exponent of power law"],
             ]
# pylint: enable=bad-whitespace, line-too-long


def Iq(q,
       scale_1=100.0,
       length_1=1000.0,
       exp_1=3.0,
       scale_2=10.0,
       length_2=100.0,
       exp_2=2.0,
       scale_3=1.0,
       length_3=10.0,
       exp_3=1.0):

    """
    :param q:                   Input q-value (float or [float, float])
    :param scale_1:     Second scale factor for broad Lorentzian peak
    :param length_1:    First Lorentzian screening length
    :param exp_1:       Exponent of the first Lorentz function
    :param scale_2:     Second scale factor for broad Lorentzian peak
    :param length_2:    Second Lorentzian screening length
    :param exp_2:       Exponent of the second Lorentz function
    :param scale_3:     Third scale factor for broad Lorentzian peak
    :param length_3:    Third Lorentzian screening length
    :param exp_3:       Exponent of the Third Lorentz function
    :return:                    Calculated intensity
    """
# pylint: disable=bad-whitespace
    intensity  = scale_1/(1.0 + power(q*length_1, exp_1))
    intensity += scale_2/(1.0 + power(q*length_2, exp_2))
    intensity += scale_3/(1.0 + power(q*length_3, exp_3))
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
        scale_1=scale[0],
        length_1=length[0],
        exp_1=expon[0],
        scale_2=scale[1],
        length_2=length[1],
        exp_2=expon[1],
        scale_3=scale[2],
        length_3=length[2],
        exp_3=expon[2],
    )
    return pars

