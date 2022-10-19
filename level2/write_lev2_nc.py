from read_specs import get_site_specs
from level2.lev2_meta_nc import get_data_attributes
from level1.quality_control import spectral_consistency
from level2.get_ret_coeff import get_mvr_coeff
from level2.lwp_offset import correct_lwp_offset
from utils import rebin_2d, get_coeff_list, isbit
import rpg_mwr
import os
from datetime import datetime, date
import pytz
from tzwhere import tzwhere
import numpy as np
from numpy import ma
import netCDF4 as nc
import pandas as pd
import metpy.calc as mpcalc
from metpy.units import units, masked_array

Fill_Value_Float = -999.0
Fill_Value_Int = -99
"""specific gas constant for water vapor (J/kg K)
and vapor pressure e0 (Pa) at T0 (K)"""
Rw, e0, T0 = 462.0, 611.0, 273.15


def lev2_to_nc(
    date: str, site: str, data_type: str, lev1_path: str, lev2_path: str
) -> dict:
    """This function reads Level 1 files,
    applies retrieval coefficients for Level 2 products and writes it into a netCDF file.

    Args:
        site: Name of site.
        data_type: Data type of the netCDF file.
        lev1_path: Path of Level 1 file.
        lev2_path: Path of Level 2 output directory.

    Examples:
        >>> from level2.write_lev2_nc import lev2_to_nc
        >>> lev2_to_nc('date', site_name', '2P00', '/path/to/lev1_file/lev1_data.nc', '/path/to/lev2_file/')
    """

    if data_type not in (
        "2P01",
        "2P02",
        "2P03",
        "2P04",
        "2P07",
        "2P08",
        "2I01",
        "2I02",
        "2S02",
    ):
        raise RuntimeError(
            ["Data type " + data_type + " not supported for file writing."]
        )
    with nc.Dataset(lev1_path) as lev1:
        BL_scan = _test_BL_scan(site, lev1)
        if (data_type == "2P02") & (BL_scan == False):
            data_type = "2P01"
        if data_type in ("2P04", "2P07", "2P08"):
            T_prod = "2P02"
            if BL_scan == False:
                T_prod = "2P01"
            for d_type in [T_prod, "2P03"]:
                global_attributes, params = get_site_specs(site, d_type)
                if not os.path.isfile(
                    lev2_path
                    + "MWR_"
                    + d_type
                    + "_"
                    + global_attributes["wigos_station_id"]
                    + "_"
                    + date
                    + ".nc"
                ):
                    rpg_dat, coeff, index, flag = get_products(
                        site, lev1, d_type, params, global_attributes
                    )
                    _combine_lev1(lev1, rpg_dat, index)
                    _mask_flag(rpg_dat, flag)
                    _add_att(global_attributes, coeff, 0)
                    hatpro = rpg_mwr.Rpg(rpg_dat)
                    hatpro.data = get_data_attributes(hatpro.data, d_type)
                    output_file = (
                        lev2_path
                        + "MWR_"
                        + d_type
                        + "_"
                        + global_attributes["wigos_station_id"]
                        + "_"
                        + date
                        + ".nc"
                    )
                    rpg_mwr.save_rpg(
                        hatpro, output_file, global_attributes, d_type, params
                    )

        global_attributes, params = get_site_specs(site, data_type)
        rpg_dat, coeff, index, flag = get_products(
            site, lev1, data_type, params, global_attributes
        )
        _combine_lev1(lev1, rpg_dat, index)
        _mask_flag(rpg_dat, flag)
        _add_att(global_attributes, coeff, 0)
        hatpro = rpg_mwr.Rpg(rpg_dat)
        hatpro.data = get_data_attributes(hatpro.data, data_type)
        output_file = (
            lev2_path
            + "MWR_"
            + data_type
            + "_"
            + global_attributes["wigos_station_id"]
            + "_"
            + date
            + ".nc"
        )
        rpg_mwr.save_rpg(hatpro, output_file, global_attributes, data_type, params)


def get_products(
    site: str, lev1: dict, data_type: str, params: dict, global_attributes: dict
) -> dict:
    "Derive specified Level 2 products"

    rpg_dat = dict()

    if data_type in ("2I01", "2I02"):
        if data_type == "2I01":
            prod = "lwp"
        else:
            prod = "iwv"

        coeff, _, _, _, _, _ = get_mvr_coeff(site, prod, lev1["frequency"][:])
        if coeff["rt"] < 2:
            coeff, offset, lin, quad, ran_err, sys_err = get_mvr_coeff(
                site, prod, lev1["frequency"][:]
            )
        else:
            coeff, ns_ta, ns_sc, ns_os, w1, w2 = get_mvr_coeff(
                site, prod, lev1["frequency"][:]
            )
        ret_in = retrieval_input(lev1, coeff)

        _, freq_ind, _ = np.intersect1d(
            lev1["frequency"][:],
            coeff["freq"][:, 0],
            assume_unique=False,
            return_indices=True,
        )
        index = np.where(
            (lev1["pointing_flag"][:] == 0)
            & np.any(
                np.abs(
                    (np.ones((len(lev1["ele"][:]), len(coeff["ele"]))) * coeff["ele"])
                    - np.transpose(
                        np.ones((len(coeff["ele"]), len(lev1["ele"][:])))
                        * lev1["ele"][:]
                    )
                )
                < 0.6,
                axis=1,
            )
        )[0]
        if params["flag_status"][3] == 1:
            flag = np.where(
                np.sum(isbit(lev1["quality_flag"][index, freq_ind], 3), axis=1) > 0
            )[0]
        else:
            flag = np.where(np.sum(lev1["quality_flag"][index, freq_ind], axis=1) > 0)[
                0
            ]

        if len(index) == 0:
            raise RuntimeError(
                ["No suitable data found for processing for data type: " + data_type]
            )

        if coeff["rt"] < 2:
            coeff_offset = offset(lev1["ele"][index])
            coeff_lin = lin(lev1["ele"][index])
            coeff_quad = quad(lev1["ele"][index])
            if prod == "lwp":
                lwp = (
                    coeff_offset
                    + np.sum(coeff_lin * ret_in[index, :], axis=1)
                    + np.sum(coeff_quad * ret_in[index, :] ** 2, axis=1)
                )
            else:
                rpg_dat["iwv"] = (
                    coeff_offset
                    + np.sum(coeff_lin * ret_in[index, :], axis=1)
                    + np.sum(coeff_quad * ret_in[index, :] ** 2, axis=1)
                )
            rpg_dat[str(prod + "_random_error")] = ran_err(lev1["ele"][index])
            rpg_dat[str(prod + "_systematic_error")] = sys_err(lev1["ele"][index])

        # else:
        #     c_w1 = w1(lev1["ele"][index])
        #     c_w2 = w2(lev1["ele"][index])
        #     hl = np.ones((len(index), c_w1.shape[2]), np.float32) * Fill_Value_Float
        #     for ind in range(c_w1.shape[2]):
        #         hl[:, ind] = np.tanh(
        #             np.sum(
        #                 (
        #                     ret_in[index, 1:] * coeff["in_sc"][0, :]
        #                     + coeff["in_os"][0, :]
        #                 )
        #                 * c_w1[:, 1:, ind]
        #             )
        #             + c_w1[:, 0, ind]
        #         )

        if prod == "lwp":
            freq_31 = np.where(np.round(lev1["frequency"][:], 1) == 31.4)[0]
            if len(freq_31) != 1:
                rpg_dat["lwp"], rpg_dat["lwp_offset"] = (
                    lwp,
                    np.ones(len(index)) * Fill_Value_Float,
                )
            else:
                rpg_dat["lwp"], rpg_dat["lwp_offset"] = correct_lwp_offset(
                    lev1.variables, lwp, index, site
                )

    elif data_type in ("2P01", "2P03"):
        if data_type == "2P01":
            prod, ret = "temperature", "tze"
        else:
            prod, ret = "water_vapor_vmr", "hze"

        coeff, _, _, _, _, _ = get_mvr_coeff(site, ret, lev1["frequency"][:])
        if coeff["rt"] < 2:
            coeff, offset, lin, quad, ran_err, sys_err = get_mvr_coeff(
                site, ret, lev1["frequency"][:]
            )
        ret_in = retrieval_input(lev1, coeff)

        _, freq_ind, _ = np.intersect1d(
            lev1["frequency"][:],
            coeff["freq"][:, 0],
            assume_unique=False,
            return_indices=True,
        )
        index = np.where(
            (lev1["pointing_flag"][:] == 0)
            & np.any(
                np.abs(
                    (np.ones((len(lev1["ele"][:]), len(coeff["ele"]))) * coeff["ele"])
                    - np.transpose(
                        np.ones((len(coeff["ele"]), len(lev1["ele"][:])))
                        * lev1["ele"][:]
                    )
                )
                < 0.6,
                axis=1,
            )
        )[0]
        if params["flag_status"][3] == 1:
            flag = np.where(
                np.sum(isbit(lev1["quality_flag"][index, freq_ind], 3), axis=1) > 0
            )[0]
        else:
            flag = np.where(np.sum(lev1["quality_flag"][index, freq_ind], axis=1) > 0)[
                0
            ]

        if len(index) == 0:
            raise RuntimeError(
                ["No suitable data found for processing for data type: " + data_type]
            )

        rpg_dat["altitude"] = coeff["height_grid"][:] + params["station_altitude"]
        rpg_dat[prod] = ma.masked_all((len(index), coeff["n_height_grid"]))

        if coeff["rt"] < 2:
            coeff_offset = offset(lev1["ele"][index])
            coeff_lin = lin(lev1["ele"][index])
            coeff_quad = quad(lev1["ele"][index])

            for ialt, _ in enumerate(rpg_dat["altitude"]):
                rpg_dat[prod][:, ialt] = (
                    coeff_offset[:, ialt]
                    + np.sum(coeff_lin[:, ialt, :] * ret_in[index, :], axis=1)
                    + np.sum(coeff_quad[:, ialt, :] * ret_in[index, :] ** 2, axis=1)
                )
            rpg_dat[str(prod + "_random_error")] = ran_err(lev1["ele"][index])
            rpg_dat[str(prod + "_systematic_error")] = sys_err(lev1["ele"][index])

    elif data_type == "2P02":

        coeff, _, _, _, _, _ = get_mvr_coeff(site, "tel", lev1["frequency"][:])
        if coeff["rt"] < 2:
            coeff, offset, lin, quad, ran_err, sys_err = get_mvr_coeff(
                site, "tel", lev1["frequency"][:]
            )

        _, freq_ind, _ = np.intersect1d(
            lev1["frequency"][:],
            coeff["freq"],
            assume_unique=False,
            return_indices=True,
        )
        _, freq_bl, _ = np.intersect1d(
            coeff["freq"], coeff["freq_bl"], assume_unique=False, return_indices=True
        )

        ix0 = np.where(
            (lev1["ele"][:] > coeff["ele"][0] - 0.6)
            & (lev1["ele"][:] < coeff["ele"][0] + 0.6)
            & (lev1["pointing_flag"][:] == 1)
        )[0]
        if len(ix0) == 0:
            raise RuntimeError(
                ["No suitable data found for processing for data type: " + data_type]
            )

        ele, flag = [], []
        for i in range(len(coeff["ele"])):
            ele.append(lev1["ele"][ix0 + i])
            flag.append(
                np.sum(isbit(lev1["quality_flag"][ix0 + i, freq_ind], 3), axis=1)
            )
        ele = np.stack(ele) * np.ones((len(coeff["ele"]), len(lev1["ele"][ix0])))
        flag = np.stack(flag) * np.ones((len(coeff["ele"]), len(lev1["ele"][ix0])))

        index = np.where(
            np.any(
                np.abs(
                    (np.ones((len(lev1["ele"][ix0]), len(coeff["ele"]))) * coeff["ele"])
                    - np.transpose(ele)
                )
                < 0.6,
                axis=1,
            )
        )[0]
        index = ix0[index]
        flag = np.where(pd.DataFrame((flag > 0)).all())[0]

        if len(index) == 0:
            raise RuntimeError(
                ["No suitable data found for processing for data type: " + data_type]
            )

        tb = (
            np.ones((len(coeff["freq"]), len(coeff["ele"]), len(index)))
            * Fill_Value_Float
        )
        for ii, iv in enumerate(index):
            tb[:, :, ii] = lev1["tb"][iv : iv + len(coeff["ele"]), freq_ind].T

        tb_alg = []
        if len(freq_ind) - len(freq_bl) > 0:
            tb_alg = np.squeeze(tb[0 : len(freq_ind) - len(freq_bl), 0, :])
        for ifq, _ in enumerate(coeff["freq_bl"]):
            tb_alg = np.append(tb_alg, np.squeeze(tb[freq_bl[ifq], :, :]), axis=0)

        rpg_dat["altitude"] = coeff["height_grid"][:] + params["station_altitude"]
        rpg_dat["temperature"] = ma.masked_all((len(index), coeff["n_height_grid"]))

        if coeff["rt"] < 2:
            for ialt, _ in enumerate(rpg_dat["altitude"]):
                rpg_dat["temperature"][:, ialt] = offset[ialt] + np.sum(
                    lin[:, ialt] * tb_alg[:, :].T, axis=1
                )

            rpg_dat["temperature_random_error"] = (
                np.repeat(ran_err, len(index)).reshape((len(ran_err), len(index))).T
            )
            rpg_dat["temperature_systematic_error"] = (
                np.repeat(sys_err, len(index)).reshape((len(sys_err), len(index))).T
            )

    elif data_type in ("2P04", "2P07", "2P08"):

        tem_dat = load_product(site, datetime.strptime(lev1.date, "%Y-%m-%d"), "2P02")
        hum_dat = load_product(site, datetime.strptime(lev1.date, "%Y-%m-%d"), "2P03")
        coeff, index = dict(), []
        _add_att(coeff, tem_dat, 1)

        hum = rebin_2d(
            hum_dat.variables["time"][:],
            hum_dat.variables["water_vapor_vmr"][:, :],
            tem_dat.variables["time"][:],
        )
        hum_re = rebin_2d(
            hum_dat.variables["time"][:],
            hum_dat.variables["water_vapor_vmr_random_error"][:, :],
            tem_dat.variables["time"][:],
        )
        hum_se = rebin_2d(
            hum_dat.variables["time"][:],
            hum_dat.variables["water_vapor_vmr_systematic_error"][:, :],
            tem_dat.variables["time"][:],
        )
        flag = ma.any(
            np.ma.mask_or(
                ma.getmask(tem_dat.variables["temperature"][:, :]),
                ma.getmask(hum[0]),
                shrink=False,
            ),
            axis=1,
        )

        rpg_dat["altitude"] = tem_dat.variables["altitude"][:]
        pres = np.interp(
            tem_dat.variables["time"][:], lev1["time"][:], lev1["air_pressure"][:]
        )
        p_baro = calc_p_baro(
            tem_dat.variables["temperature"][:, :], hum[0], pres, rpg_dat["altitude"]
        )
        if data_type == "2P04":
            rpg_dat["relative_humidity"] = (
                vap_pres(hum[0], tem_dat.variables["temperature"][:, :])
                / mpcalc.saturation_vapor_pressure(
                    masked_array(tem_dat.variables["temperature"][:, :], data_units="K")
                ).magnitude
            )
            rpg_dat["relative_humidity_random_error"] = rh_err(
                tem_dat.variables["temperature"][:, :],
                hum[0],
                tem_dat.variables["temperature_random_error"],
                hum_re[0],
            )
            rpg_dat["relative_humidity_systematic_error"] = rh_err(
                tem_dat.variables["temperature"][:, :],
                hum[0],
                tem_dat.variables["temperature_systematic_error"],
                hum_se[0],
            )
        if data_type == "2P07":
            rpg_dat["potential_temperature"] = mpcalc.potential_temperature(
                masked_array(p_baro, data_units="Pa"),
                masked_array(tem_dat.variables["temperature"][:, :], data_units="K"),
            ).magnitude
        if data_type == "2P08":
            Theta = mpcalc.potential_temperature(
                masked_array(p_baro, data_units="Pa"),
                masked_array(tem_dat.variables["temperature"][:, :], data_units="K"),
            ).magnitude
            e = vap_pres(hum[0], tem_dat.variables["temperature"][:, :])
            rpg_dat["equivalent_potential_temperature"] = (
                Theta
                + (
                    spec_heat(tem_dat.variables["temperature"][:, :])
                    * 0.622
                    * e
                    / (p_baro - e)
                    / 1004.0
                )
                * Theta
                / tem_dat.variables["temperature"][:, :]
            )
        _combine_lev1(
            tem_dat.variables, rpg_dat, np.arange(len(tem_dat.variables["time"][:]))
        )

    elif data_type == "2S02":

        c_list = get_coeff_list(site, "tbx")
        coeff, _, _, _, _, _ = get_mvr_coeff(site, "tbx", lev1["frequency"][:])
        index = np.where(
            (lev1["ele"][:] > coeff["ele"] - 0.6)
            & (lev1["ele"][:] < coeff["ele"] + 0.6)
            & (lev1["pointing_flag"][:] == 0)
        )[0]
        _, tb_ret = spectral_consistency(lev1, c_list)
        rpg_dat["tb_spectrum"] = tb_ret[index, :]
        rpg_dat["tb"] = lev1["tb"][index, :]
        fields = ["frequency", "receiver_nb", "receiver"]
        for name in fields:
            rpg_dat[name] = lev1[name][:]
        flag = np.array(np.zeros(len(index)), dtype=bool)

    return rpg_dat, coeff, index, flag


def _combine_lev1(lev1: dict, rpg_dat: dict, index: np.ndarray) -> None:
    "add level1 data"
    lev1_vars = [
        "time",
        "time_bnds",
        "station_latitude",
        "station_longitude",
        "azi",
        "ele",
    ]
    if index != []:
        for ivars in lev1_vars:
            if lev1[ivars].ndim > 1:
                rpg_dat[ivars] = lev1[ivars][index, :]
            else:
                rpg_dat[ivars] = lev1[ivars][index]


def _mask_flag(rpg_dat: dict, flag: np.ndarray) -> None:
    "mask flagged data"
    for ivars in rpg_dat:
        if (rpg_dat[ivars].ndim > 1) & (len(rpg_dat[ivars]) == len(rpg_dat["time"])):
            rpg_dat[ivars][flag, :] = ma.masked
        elif (
            (rpg_dat[ivars].ndim == 1)
            & (len(rpg_dat[ivars]) == len(rpg_dat["time"]))
            & (ivars != "time")
        ):
            rpg_dat[ivars][flag] = ma.masked


def _add_att(global_attributes: dict, coeff: dict, lev2: int) -> None:
    "add retrieval attributes"
    fields = [
        "retrieval_type",
        "retrieval_elevation_angles",
        "retrieval_frequencies",
        "retrieval_auxiliary_input",
        "retrieval_description",
    ]
    for name in fields:
        if lev2:
            global_attributes[name] = eval("coeff." + name)
        else:
            global_attributes[name] = coeff[name]


def load_product(site: str, date: str, prod: str):
    "load existing lev2 file for deriving other products"
    file = []
    global_attributes, params = get_site_specs(site, "1C01")
    ID = global_attributes["wigos_station_id"]
    data_out_l2 = params["data_out"] + "level2/" + date.strftime("%Y/%m/%d/")
    file_name = (
        data_out_l2 + "MWR_" + prod + "_" + ID + "_" + date.strftime("%Y%m%d") + ".nc"
    )

    if os.path.isfile(file_name):
        file = nc.Dataset(file_name)
    "load single pointing T if no BL scans are performed"
    if (not os.path.isfile(file_name)) & (prod == "2P02"):
        file_name = (
            data_out_l2 + "MWR_2P01_" + ID + "_" + date.strftime("%Y%m%d") + ".nc"
        )
        if os.path.isfile(file_name):
            file = nc.Dataset(file_name)
    return file


def _test_BL_scan(site: str, lev1: dict) -> bool:
    "Check for existing BL scans in lev1 data"
    BL_scan = True
    coeff, _, _, _, _, _ = get_mvr_coeff(site, "tel", lev1["frequency"][:])
    ix0 = np.where(
        (lev1["ele"][:] > coeff["ele"][0] - 0.6)
        & (lev1["ele"][:] < coeff["ele"][0] + 0.6)
        & (lev1["pointing_flag"][:] == 1)
    )[0]
    if len(ix0) == 0:
        BL_scan = False
    return BL_scan


"specific heat for evaporation (J/kg)"


def spec_heat(T):
    return (2500.0 - 2.42 * (T - T0)) * 1000.0


"water vapor pressure"


def vap_pres(qs, T):
    return qs * Rw * T


def rh_err(T: np.ndarray, qs: np.ndarray, dT: np.ndarray, dq: np.ndarray) -> np.ndarray:
    "Calculates relative humidity error from absolute humidity and temperature"

    L = spec_heat(T)
    e = vap_pres(qs, T)
    es = mpcalc.saturation_vapor_pressure(masked_array(T, data_units="K")).magnitude
    # es = e0*np.exp((L/(Rw*T0))*((T-T0)/T))

    "error propagation"
    drh_dq = Rw * T / es
    # des_dT = es*(T*(T-T0)*(-2420)+T0*L)/(Rw*T0*T**2)
    des_dT = es * 17.67 * 243.5 / ((T - T0) + 243.5) ** 2
    drh_dT = qs * Rw / es**2 * (es - T * des_dT)
    drh = np.sqrt((drh_dq * dq) ** 2 + (drh_dT * dT) ** 2)

    return drh


def calc_p_baro(
    T: np.ndarray, qs: np.ndarray, p: np.ndarray, z: np.ndarray
) -> np.ndarray:
    "Calculate pressure in each level using barometric height formula"

    Tv = T * (1.0 + 0.608 * qs)
    p_baro = ma.masked_all(T.shape)
    p_baro[
        (~ma.getmaskarray(qs).any(axis=1)) & (~ma.getmaskarray(T).any(axis=1)), 0
    ] = (
        p[(~ma.getmaskarray(qs).any(axis=1)) & (~ma.getmaskarray(T).any(axis=1))]
        * 100.0
    )
    for ialt in np.arange(len(z) - 1) + 1:
        p_baro[:, ialt] = p_baro[:, ialt - 1] * np.exp(
            -9.81
            * (z[ialt] - z[ialt - 1])
            / (287.0 * (Tv[:, ialt] + Tv[:, ialt - 1]) / 2.0)
        )

    return p_baro


def retrieval_input(lev1: dict, cf: list) -> np.ndarray:
    "Get retrieval input"

    bias = np.ones((len(lev1["time"][:]), 1), np.float32)
    doy = np.ones((len(lev1["time"][:]), 2), np.float32) * Fill_Value_Float
    tz = tzwhere.tzwhere()
    sun = np.ones((len(lev1["time"][:]), 2), np.float32) * Fill_Value_Float
    for i, time in enumerate(lev1["time"][:].data):
        doy[i, 0] = np.cos(
            datetime.fromtimestamp(time).timetuple().tm_yday / 365 * 2 * np.pi
        )
        doy[i, 1] = np.sin(
            datetime.fromtimestamp(time).timetuple().tm_yday / 365 * 2 * np.pi
        )
        timezone_str = tz.tzNameAt(
            lev1["station_latitude"][i], lev1["station_longitude"][i]
        )
        timezone = pytz.timezone(timezone_str)
        dt = datetime.fromtimestamp(time)
        dt = dt + timezone.utcoffset(dt)
        sun[i, 0] = np.cos(
            (dt.hour + dt.minute / 60 + dt.second / 3600) / 365 * 2 * np.pi
        )
        sun[i, 1] = np.sin(
            (dt.hour + dt.minute / 60 + dt.second / 3600) / 365 * 2 * np.pi
        )

    fields = [
        "air_temperature",
        "relative_humidity",
        "air_pressure",
        "irt",
        "doy",
        "sun",
    ]

    if cf["rt"] < 2:
        ret_in = lev1["tb"][:, :]
        for field in fields:
            if field in lev1.variables:
                if lev1[field].ndim == 1:
                    ret_in = np.concatenate(
                        (
                            ret_in,
                            np.reshape(lev1[field][:].data, (len(lev1["time"][:]), 1)),
                        ),
                        axis=1,
                    )
                else:
                    ret_in = np.concatenate((ret_in, lev1[field][:, :]), axis=1)
            else:
                ret_in = np.concatenate((ret_in, eval(field)), axis=1)
    else:
        ret_in = np.concatenate((bias, lev1["tb"][:, :]), axis=1)
        for i, field in enumerate(cf["aux"]):
            if (field == "TS") & (cf["aux_flg"][i] == 1):
                ret_in = np.concatenate(
                    (
                        ret_in,
                        np.reshape(
                            lev1["air_temperature"][:].data, (len(lev1["time"][:]), 1)
                        ),
                    ),
                    axis=1,
                )
            if (field == "HS") & (cf["aux_flg"][i] == 1):
                ret_in = np.concatenate(
                    (
                        ret_in,
                        np.reshape(
                            lev1["relative_humidity"][:].data, (len(lev1["time"][:]), 1)
                        ),
                    ),
                    axis=1,
                )
            if (field == "PS") & (cf["aux_flg"][i] == 1):
                ret_in = np.concatenate(
                    (
                        ret_in,
                        np.reshape(
                            lev1["air_pressure"][:].data, (len(lev1["time"][:]), 1)
                        ),
                    ),
                    axis=1,
                )
            if (field == "ZS") & (cf["aux_flg"][i] == 1):
                ret_in = np.concatenate((ret_in, lev1["irt"][:, :]), axis=1)
            if (field == "IR") & (cf["aux_flg"][i] == 1):
                ret_in = np.concatenate((ret_in, lev1["irt"][:, :]), axis=1)
            if (field == "I1") & (cf["aux_flg"][i] == 1):
                ret_in = np.concatenate(
                    (
                        ret_in,
                        np.reshape(lev1["irt"][:, 0].data, (len(lev1["time"][:]), 1)),
                    ),
                    axis=1,
                )
            if (field == "I2") & (cf["aux_flg"][i] == 1):
                ret_in = np.concatenate(
                    (
                        ret_in,
                        np.reshape(lev1["irt"][:, 1].data, (len(lev1["time"][:]), 1)),
                    ),
                    axis=1,
                )
            if (field == "DY") & (cf["aux_flg"][i] == 1):
                ret_in = np.concatenate((ret_in, doy), axis=1)
            if (field == "SU") & (cf["aux_flg"][i] == 1):
                ret_in = np.concatenate((ret_in, sun), axis=1)

    return ret_in
