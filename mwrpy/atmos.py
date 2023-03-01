"""Module for atmsopheric functions."""
import numpy as np
from numpy import ma
import metpy.calc as mpcalc
from metpy.units import masked_array
import scipy.constants
import mwrpy.constants as con

HPA_TO_P = 100


def spec_heat(T: np.ndarray) -> np.ndarray:
    "Specific heat for evaporation (J/kg)"
    return con.LATENT_HEAT - 2420.0 * (T - con.T0)


def vap_pres(q: np.ndarray, T: np.ndarray) -> np.ndarray:
    "Water vapor pressure (Pa)"
    return q * con.RW * T


def t_dew_rh(T: np.ndarray, rh: np.ndarray) -> np.ndarray:
    "Dew point temperature (K) from relative humidity ()"
    return (
        mpcalc.dewpoint_from_relative_humidity(
            masked_array(T, data_units="K"), masked_array(rh, data_units="")
        ).magnitude
        + con.T0
    )


def pot_tem(T: np.ndarray, q: np.ndarray, p: np.ndarray, z: np.ndarray) -> np.ndarray:
    "Potential temperature (K)"
    p_baro = calc_p_baro(T, q, p, z)
    return mpcalc.potential_temperature(
        masked_array(p_baro, data_units="Pa"), masked_array(T, data_units="K")
    ).magnitude


def eq_pot_tem(T: np.ndarray, q: np.ndarray, p: np.ndarray, z: np.ndarray) -> np.ndarray:
    "Equivalent potential temperature (K)"
    e = vap_pres(q, T)
    p_baro = calc_p_baro(T, q, p, z)
    Theta = pot_tem(T, q, p, z)
    return Theta + (spec_heat(T) * con.MW_RATIO * e / (p_baro - e) / con.SPECIFIC_HEAT) * Theta / T


def rel_hum(T: np.ndarray, q: np.ndarray) -> np.ndarray:
    "Relative humidity ()"
    return vap_pres(q, T) / calc_saturation_vapor_pressure(T)


def rh_err(T: np.ndarray, q: np.ndarray, dT: np.ndarray, dq: np.ndarray) -> np.ndarray:
    "Calculates relative humidity error propagation from absolute humidity and temperature ()"
    es = calc_saturation_vapor_pressure(T)
    drh_dq = con.RW * T / es
    des_dT = es * 17.67 * 243.5 / ((T - con.T0) + 243.5) ** 2
    drh_dT = q * con.RW / es**2 * (es - T * des_dT)
    drh = np.sqrt((drh_dq * dq) ** 2 + (drh_dT * dT) ** 2)

    return drh


def abs_hum(T: np.ndarray, rh: np.ndarray) -> np.ndarray:
    "Absolute humidity (kg/m^3)"
    es = calc_saturation_vapor_pressure(T)
    return (rh * es) / (con.RW * T)


def calc_p_baro(T: np.ndarray, q: np.ndarray, p: np.ndarray, z: np.ndarray) -> np.ndarray:
    "Calculate pressure (Pa) in each level using barometric height formula"
    Tv = mpcalc.virtual_temperature(
        masked_array(T, data_units="K"), masked_array(q, data_units="")
    ).magnitude
    p_baro = ma.masked_all(T.shape)
    p_baro[(~ma.getmaskarray(q).any(axis=1)) & (~ma.getmaskarray(T).any(axis=1)), 0] = (
        p[(~ma.getmaskarray(q).any(axis=1)) & (~ma.getmaskarray(T).any(axis=1))] * HPA_TO_P
    )
    for ialt in np.arange(len(z) - 1) + 1:
        p_baro[:, ialt] = p_baro[:, ialt - 1] * np.exp(
            -scipy.constants.g
            * (z[ialt] - z[ialt - 1])
            / (con.RS * np.mean([Tv[:, ialt], Tv[:, ialt - 1]]))
        )

    return p_baro


def calc_saturation_vapor_pressure(temperature: np.ndarray) -> np.ndarray:
    """Goff-Gratch formula for saturation vapor pressure (Pa) over water adopted by WMO.
    Args:
        temperature: Temperature (K).
    Returns:
        Saturation vapor pressure (Pa).
    """
    ratio = con.T0 / temperature
    inv_ratio = ratio**-1
    return (
        10
        ** (
            10.79574 * (1 - ratio)
            - 5.028 * np.log10(inv_ratio)
            + 1.50475e-4 * (1 - (10 ** (-8.2969 * (inv_ratio - 1))))
            + 0.42873e-3 * (10 ** (4.76955 * (1 - ratio)) - 1)
            + 0.78614
        )
    ) * HPA_TO_P


def c2k(T: np.ndarray) -> np.ndarray:
    """Converts Celsius to Kelvins."""
    return ma.array(T) + con.T0


def dir_avg(time: np.ndarray, spd: np.ndarray, drc: np.ndarray, win: float = 0.5):
    """Computes average wind direction (DEG) for a certain window length"""
    width = len(time[time <= time[0] + win])
    if (width % 2) != 0:
        width = width + 1
    seq = range(len(time))
    avg_dir = []
    for i in range(len(seq) - width + 1):
        avg_dir.append(winddir(spd[seq[i : i + width]], drc[seq[i : i + width]]))
    return np.array(avg_dir), width


def winddir(spd: np.ndarray, drc: np.ndarray):
    """Computes mean wind direction (deg)"""
    ve = -np.mean(spd * np.sin(np.deg2rad(drc)))
    vn = -np.mean(spd * np.cos(np.deg2rad(drc)))
    vdir = np.rad2deg(np.arctan2(ve, vn))
    if vdir < 180.0:
        Dv = vdir + 180.0
    elif vdir > 180.0:
        Dv = vdir - 180
    return Dv
