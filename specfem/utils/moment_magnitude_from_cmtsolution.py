import sys
import numpy as np


def seismic_moment(mt):
    """
    Return the seismic moment based on a moment tensor.
    Can take a list of tensor components, or a Tensor object from ObsPy.

    :type mt: list of floats or obspy.core.event.source.Tensor
    :param mt: the components of the moment tensor M_ij
    :rtype: float
    :return: the seismic moment, in units of N*m
    """
    return 1 / np.sqrt(2) * np.sqrt(sum([_ ** 2 for _ in mt]))


def moment_magnitude(moment, c=10.7):
    """
    Return the moment magitude, M_w, based on a seismic moment. Equation from
    Hanks & Kanamori (1979)

    :type c: float
    :param c: correction factor for conversion, 10.7 for units of N*m,
        16.1 for units of dyne*cm
    :type moment: float
    :param moment: the seismic moment, in units of N*m
    :rtype: float
    :return: moment magnitude, M_w
    """
    return 2 / 3 * np.log10(moment) - c

if __name__ == "__main__":
    with open(sys.argv[1], "r") as f:
        lines = f.readlines()

    fm = []
    for line in lines[7:]:
        print(line)
        key, val = line.strip().split(":")
        fm.append(float(val))

    print(moment_magnitude(seismic_moment(fm))

