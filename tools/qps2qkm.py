"""
Convert Qp/Qs to Qkappa/Qmu or vice versa
Ported from Carl Tape's GEOTOOLS/matlab_util/util_vel/qps2qkm.m
"""


def negligible_qkappa(qp, qs, vp, vs):
    """
    Specfem has an equation that defines whether a converted bulk attenuation
    is negligible, if it is, then the bulk attenuation is replaced by a default
    maximum value
    
    :type qp: float or np.array
    :param qp: P-wave attenuation
    :type qs: float or np.array
    :param qs: S-wave attenuation
    :type vp: float or np.array
    :param vp: P-wave velocity for the given point, m/s or km/s
    :type vs: float or np.array
    :param vs: S-wave velocity for the given point, same units as vp
    :rtype: tuple (bool, float) or np.array of tuples
    :return: (negligible qkappa?, returned value of equation)
    """
    f = 4/3 * (vs/vp) ** 2
    neg = abs(qs - f * qp)
    return neg <= 1E-5, neg


def qps2qkm(qp_or_qk, qs_or_qm, vp, vs, method):
    """
    Convert Qp, Qs to Qkappa, Qmu, or vice versa
    Ported from Carl Tape's qps2qkm.m
    Based on Dahlen and Tromp, Eq. 9.59
    
    :type qp_or_qk: float or np.array
    :param qp_or_qk: Qp or Qkappa
    :type qs_or_qm: float or np.array
    :param qs_or_qm: Qs or Qmu
    :type vp: float or np.array
    :param vp: P-wave velocity for the given point, m/s or km/s
    :type vs: float or np.array
    :param vs: S-wave velocity for the given point, same units as vp
    :type method: str
    :param method: 'qps2qkm' for Qp, Qs to Qkappa, Qmu  OR
                   'qkm2qps' for Qkappa, Qmu to Qp, Qs
    :rtype: float, float or np.array
    :return: Qp, Qs or Qkappa, Qmu, depending on method
    """
    methods = ["qps2qkm", "qkm2qps"]
    assert(method in methods), f"method must be in {methods}"
    
    f = 4/3 * (vs/vp) ** 2
    
    if method == "qps2qkm":
        qk = (1 - f) / (1 / qp_or_qk - f / qs_or_qm)
        return qk, qs_or_qm
    elif method == "qkm2qps":
        qp = 1 / ((1 - f) / qp_or_qk + f / qs_or_qm)
        return qp, qs_or_qm


