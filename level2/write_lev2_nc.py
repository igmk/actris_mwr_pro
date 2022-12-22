import os
from datetime import datetime

import netCDF4 as nc
import numpy as np
import pytz
from numpy import ma
from tzwhere import tzwhere

import rpg_mwr
from atmos import eq_pot_tem, pot_tem, rel_hum, rh_err
from level1.quality_control import spectral_consistency
from level2.get_ret_coeff import get_mvr_coeff
from level2.lev2_meta_nc import get_data_attributes
from level2.lwp_offset import correct_lwp_offset
from read_specs import get_site_specs
from utils import add_time_bounds, get_coeff_list, get_ret_ang, get_ret_freq, interpol_2d, get_file_list
from level1.rpg_bin import get_rpg_bin

Fill_Value_Float = -999.0
Fill_Value_Int = -99


def lev2_to_nc(date_in: str, site: str, data_type: str, lev1_path: str, lev2_path: str) -> dict:
    """This function reads Level 1 files,
    applies retrieval coefficients for Level 2 products and writes it into a netCDF file.

    Args:
        site: Name of site.
        data_type: Data type of the netCDF file.
        lev1_path: Path of Level 1 file.
        lev2_path: Path of Level 2 output directory.

    Examples:
        >>> from level2.write_lev2_nc import lev2_to_nc
        >>> lev2_to_nc('date', site_name', '2P00', 
        '/path/to/lev1_file/lev1_data.nc', '/path/to/lev2_file/')
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
        raise RuntimeError(["Data type " + data_type + " not supported for file writing."])

    with nc.Dataset(lev1_path) as lev1:
        BL_scan = _test_BL_scan(site, lev1)
        if (data_type == "2P02") & (not BL_scan):
            data_type = "2P01"
        if data_type in ("2P04", "2P07", "2P08"):
            T_prod = "2P02"
            if not BL_scan:
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
                    + date_in
                    + ".nc"
                ):
                    rpg_dat, coeff, index = get_products(
                        site, lev1, d_type, params
                    )
                    _combine_lev1(lev1, rpg_dat, index, d_type, params)
                    hatpro = rpg_mwr.Rpg(rpg_dat)
                    hatpro.data = get_data_attributes(hatpro.data, d_type)
                    output_file = (
                        lev2_path
                        + "MWR_"
                        + d_type
                        + "_"
                        + global_attributes["wigos_station_id"]
                        + "_"
                        + date_in
                        + ".nc"
                    )
                    rpg_mwr.save_rpg(hatpro, output_file, global_attributes, d_type)

        global_attributes, params = get_site_specs(site, data_type)
        rpg_dat, coeff, index = get_products(site, lev1, data_type, params)
        _combine_lev1(lev1, rpg_dat, index, data_type, params)
        _add_att(global_attributes, coeff, lev1)
        hatpro = rpg_mwr.Rpg(rpg_dat)
        hatpro.data = get_data_attributes(hatpro.data, data_type)
        output_file = (
            lev2_path
            + "MWR_"
            + data_type
            + "_"
            + global_attributes["wigos_station_id"]
            + "_"
            + date_in
            + ".nc"
        )
        rpg_mwr.save_rpg(hatpro, output_file, global_attributes, data_type)


def get_products(
    site: str, lev1: dict, data_type: str, params: dict
) -> dict:
    "Derive specified Level 2 products"

    rpg_dat = {}

    if data_type in ("2I01", "2I02"):
        if data_type == "2I01":
            prod = "lwp"
        else:
            prod = "iwv"

        coeff = get_mvr_coeff(site, prod, lev1["frequency"][:])
        if coeff[0]["rt"] < 2:
            coeff, offset, lin, quad, ran_err, sys_err = get_mvr_coeff(
                site, prod, lev1["frequency"][:]
            )
        else:
            coeff, in_sc, in_os, op_sc, op_os, w1, w2, pn = get_mvr_coeff(site, prod, lev1["frequency"][:])
        ret_in = retrieval_input(lev1, coeff)
        _, freq_ind, _ = np.intersect1d(
            lev1["frequency"][:], coeff["freq"][:, 0], assume_unique=False, return_indices=True
        )

        index = np.where(
            (lev1["pointing_flag"][:] == 0)
            & np.any(
                np.abs(
                    (np.ones((len(lev1["ele"][:]), len(coeff["ele"]))) * coeff["ele"])
                    - np.transpose(
                        np.ones((len(coeff["ele"]), len(lev1["ele"][:]))) * lev1["ele"][:]
                    )
                )
                < 0.5,
                axis=1,
            )
        )[0]
        if len(index) == 0:
            raise RuntimeError(
                ["No suitable data found for processing for data type: " + data_type]
            )
        coeff["retrieval_elevation_angles"] = str(
            np.sort(np.unique(ele_ret(lev1["ele"][index], coeff)))
        )

        if coeff["rt"] < 2:
            coeff_offset = offset(lev1["ele"][index])
            coeff_lin = lin(lev1["ele"][index])
            coeff_quad = quad(lev1["ele"][index])
            tmp = (
                coeff_offset
                + np.sum(coeff_lin * ret_in[index, :], axis=1)
                + np.sum(coeff_quad * ret_in[index, :] ** 2, axis=1)
            )
            rpg_dat[str(prod + "_random_error")] = ran_err(lev1["ele"][index])
            rpg_dat[str(prod + "_systematic_error")] = sys_err(lev1["ele"][index])

        else:    
            tmp = np.ones(len(index), np.float32) * Fill_Value_Float
            c_w1, c_w2, sp = w1(lev1["ele"][index]), w2(lev1["ele"][index]), pn(lev1["ele"][index])
            i_sc, i_os = in_sc(lev1["ele"][index]), in_os(lev1["ele"][index])
            o_sc, o_os = op_sc(lev1["ele"][index]), op_os(lev1["ele"][index])
            for ix, iv in enumerate(index):
                ret_in[iv, 1:] = (ret_in[iv, 1:] - i_os[ix, :]) * i_sc[ix, :]
                hl = np.ones(c_w1.shape[2] + 1, np.float32)
                hl[1:] = np.tanh(sp[ix] * ret_in[iv, :].dot(c_w1[ix, :, :]))
                tmp[ix] = np.tanh(sp[ix] * hl.dot(c_w2[ix, :])) * o_sc[ix, :] + o_os[ix, :]

        if prod == "lwp":
            freq_31 = np.where(np.round(lev1["frequency"][:], 1) == 31.4)[0]
            if len(freq_31) != 1:
                rpg_dat["lwp"], rpg_dat["lwp_offset"] = (
                    tmp,
                    np.ones(len(index)) * Fill_Value_Float,
                )
            else:
                rpg_dat["lwp"], rpg_dat["lwp_offset"] = correct_lwp_offset(
                    lev1.variables, tmp, index, site
                )
        else:
            rpg_dat["iwv"] = tmp
            

    elif data_type in ("2P01", "2P03"):
        if data_type == "2P01":
            prod, ret = "temperature", "tze"
        else:
            prod, ret = "water_vapor_vmr", "hze"

        coeff = get_mvr_coeff(site, ret, lev1["frequency"][:])
        if coeff[0]["rt"] < 2:
            coeff, offset, lin, quad, ran_err, sys_err = get_mvr_coeff(
                site, ret, lev1["frequency"][:]
            )
        else:
            coeff, in_sc, in_os, op_sc, op_os, w1, w2, pn = get_mvr_coeff(site, ret, lev1["frequency"][:])
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
                        np.ones((len(coeff["ele"]), len(lev1["ele"][:]))) * lev1["ele"][:]
                    )
                )
                < 0.5,
                axis=1,
            )
        )[0]
        if len(index) == 0:
            raise RuntimeError(
                ["No suitable data found for processing for data type: " + data_type]
            )
        coeff["retrieval_elevation_angles"] = str(
            np.sort(np.unique(ele_ret(lev1["ele"][index], coeff)))
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
            
        else:
            c_w1, c_w2, sp = w1(lev1["ele"][index]), w2(lev1["ele"][index]), pn(lev1["ele"][index])
            i_sc, i_os = in_sc(lev1["ele"][index]), in_os(lev1["ele"][index])
            o_sc, o_os = op_sc(lev1["ele"][index]), op_os(lev1["ele"][index])

            for ix, iv in enumerate(index):
                hl_in = np.concatenate(([1.], (ret_in[iv, 1:] - i_os[ix, :]) * i_sc[ix, :]))
                hl = np.concatenate(([1.], np.tanh(sp[ix] * hl_in.dot(c_w1[ix, :, :]))))
                rpg_dat[prod][ix, :] = np.tanh(sp[ix] * hl.dot(c_w2[ix, :, :])) * o_sc[ix, :] + o_os[ix, :]
                if prod == "water_vapor_vmr":
                    rpg_dat[prod][ix, :] = rpg_dat[prod][ix, :] / 1000.
            rpg_dat[str(prod + "_random_error")] = np.ones(rpg_dat[prod].shape, np.float32) * Fill_Value_Float
            rpg_dat[str(prod + "_systematic_error")] = np.ones(rpg_dat[prod].shape, np.float32) * Fill_Value_Float
         

    elif data_type == "2P02":

        coeff = get_mvr_coeff(site, "tel", lev1["frequency"][:])
        if coeff[0]["rt"] < 2:
            coeff, offset, lin, quad, ran_err, sys_err = get_mvr_coeff(
                site, "tel", lev1["frequency"][:]
            )
        else:
            coeff, in_sc, in_os, op_sc, op_os, w1, w2, pn = get_mvr_coeff(site, "tel", lev1["frequency"][:])
        ret_in = retrieval_input(lev1, coeff)        

        _, freq_ind, _ = np.intersect1d(
            lev1["frequency"][:],
            coeff["freq"],
            assume_unique=False,
            return_indices=True,
        )
        _, freq_bl, _ = np.intersect1d(
            coeff["freq"], coeff["freq_bl"], assume_unique=False, return_indices=True
        )

        ix0 = (
            np.where(
                (lev1["ele"][:] > coeff["ele"][-1] - 0.5)
                & (lev1["ele"][:] < coeff["ele"][-1] + 0.5)
                & (lev1["pointing_flag"][:] == 1)
            )[0]
            + 1
        )
        ibl, tb = (
            np.empty([0, len(coeff["ele"])], np.int32),
            np.ones((len(freq_ind), len(coeff["ele"]), 0), np.float32) * Fill_Value_Float,
        )
        for ix0v in ix0:
            if np.allclose(
                lev1["ele"][ix0v - len(coeff["ele"]) : ix0v], coeff["ele"], atol=0.5
            ):             
                ibl = np.append(ibl, [np.array(range(ix0v - len(coeff["ele"]), ix0v))], axis=0)
                tb = np.concatenate(
                    (
                        tb,
                        np.expand_dims(
                            lev1["tb"][ix0v - len(coeff["ele"]) : ix0v, freq_ind].T, 2
                        ),
                    ),
                    axis=2,
                )               

        if len(ibl) == 0:
            raise RuntimeError(
                ["No suitable data found for processing for data type: " + data_type]
            )
            
        index = ibl[:, -1]
        rpg_dat["altitude"] = coeff["height_grid"][:] + params["station_altitude"]
        rpg_dat["temperature"] = ma.masked_all((len(index), coeff["n_height_grid"]))                  
        
        if coeff["rt"] < 2:
            tb_alg = []
            if len(freq_ind) - len(freq_bl) > 0:
                tb_alg = np.squeeze(tb[0 : len(freq_ind) - len(freq_bl), 0, :])
            for ifq, _ in enumerate(coeff["freq_bl"]):
                tb_alg = np.append(tb_alg, np.squeeze(tb[freq_bl[ifq], :, :]), axis=0)

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
        else:
            for ix in range(ibl.shape[0]):
                ret_array = np.reshape(ret_in[np.ix_(ibl[ix, :], freq_ind+1)], len(coeff["ele"]) * len(freq_ind))
                for i_add in range(ret_in.shape[1] - len(coeff["freq"]) - 1, 0, -1):
                    ret_array = np.concatenate((ret_array, [ma.median(ret_in[ibl[ix, :], -i_add])]))
                hl_in = np.concatenate(([1.], (ret_array - coeff["in_os"][0, :]) * coeff["in_sc"][0, :]))
                hl = np.concatenate(([1.], np.tanh(coeff["np"][0] * hl_in.dot(coeff["w1"][0, :, :]))))
                rpg_dat["temperature"][ix, :] = np.tanh(coeff["np"][0] * hl.dot(coeff["w2"][0, :, :])) * coeff["op_sc"][0, :] + coeff["op_os"][0, :]
            rpg_dat["temperature_random_error"] = np.ones(rpg_dat["temperature"].shape, np.float32) * Fill_Value_Float
            rpg_dat["temperature_systematic_error"] = np.ones(rpg_dat["temperature"].shape, np.float32) * Fill_Value_Float   
            
        # file_list_tpb = get_file_list("/data/obs/site/cgn/foghat/l2/2022/11/23/", "", "", "tpb")
        # rpg_bin = get_rpg_bin(file_list_tpb)            
        # tem = nc.Dataset("/data/obs/site/cgn/foghat/l2/2022/11/23/sups_cgn_mwrBL00_l2_ta_p00_20221123000908.nc")
        # import matplotlib.pyplot as plt
        # plt.plot(rpg_bin.data["time"], rpg_bin.data["T"][:,0],"o")
        # plt.plot(lev1["time"][ibl[:, -1]], rpg_dat["temperature"][:,0],"o")
        # plt.plot(tem["time"][:], tem["ta"][:,0],"o")
        # plt.show()            
        # import pdb
        # pdb.set_trace()
                

    elif data_type in ("2P04", "2P07", "2P08"):

        tem_dat, tem_freq, tem_ang, prod = load_product(
            site, datetime.strptime(lev1.date, "%Y-%m-%d"), "2P02"
        )
        hum_dat, hum_freq, hum_ang, _ = load_product(
            site, datetime.strptime(lev1.date, "%Y-%m-%d"), "2P03"
        )

        coeff, index = {}, []
        coeff["retrieval_frequencies"] = str(
            np.unique(np.sort(np.concatenate([tem_freq, hum_freq])))
        )
        coeff["retrieval_elevation_angles"] = str(
            np.unique(np.sort(np.concatenate([tem_ang, hum_ang])))
        )
        coeff["retrieval_description"] = "derived product from: " + prod + ", 2P03"
        coeff["dependencies"] = prod + ", 2P03"

        hum = interpol_2d(
            hum_dat.variables["time"][:],
            hum_dat.variables["water_vapor_vmr"][:, :],
            tem_dat.variables["time"][:],
        )
        hum_re = interpol_2d(
            hum_dat.variables["time"][:],
            hum_dat.variables["water_vapor_vmr_random_error"][:, :],
            tem_dat.variables["time"][:],
        )
        hum_se = interpol_2d(
            hum_dat.variables["time"][:],
            hum_dat.variables["water_vapor_vmr_systematic_error"][:, :],
            tem_dat.variables["time"][:],
        )

        rpg_dat["altitude"] = tem_dat.variables["altitude"][:]
        pres = np.interp(tem_dat.variables["time"][:], lev1["time"][:], lev1["air_pressure"][:])
        if data_type == "2P04":
            rpg_dat["relative_humidity"] = rel_hum(tem_dat.variables["temperature"][:, :], hum[0])
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
            rpg_dat["potential_temperature"] = pot_tem(
                tem_dat.variables["temperature"][:, :], hum[0], pres, rpg_dat["altitude"]
            )
        if data_type == "2P08":
            rpg_dat["equivalent_potential_temperature"] = eq_pot_tem(
                tem_dat.variables["temperature"][:, :], hum[0], pres, rpg_dat["altitude"]
            )

        _combine_lev1(
            tem_dat.variables,
            rpg_dat,
            np.arange(len(tem_dat.variables["time"][:])),
            data_type,
            params,
        )

    elif data_type == "2S02":

        c_list = get_coeff_list(site, "tbx")
        coeff, _, _, _, _, _ = get_mvr_coeff(site, "tbx", lev1["frequency"][:])
        index = np.where(
            (lev1["ele"][:] > coeff["ele"] - 0.5)
            & (lev1["ele"][:] < coeff["ele"] + 0.5)
            & (lev1["pointing_flag"][:] == 0)
        )[0]
        _, tb_ret = spectral_consistency(lev1, c_list)
        rpg_dat["tb_spectrum"] = tb_ret[index, :]
        rpg_dat["tb"] = lev1["tb"][index, :]
        fields = ["frequency", "receiver_nb", "receiver"]
        for name in fields:
            rpg_dat[name] = lev1[name][:]
        rpg_dat["data_quality"] = np.array(np.zeros(len(index)), dtype=np.int32)

    return rpg_dat, coeff, index


def _combine_lev1(
    lev1: dict, rpg_dat: dict, index: np.ndarray, data_type: str, params: dict
) -> None:
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
            if (ivars == "time_bnds") & (data_type == "2P02"):
                rpg_dat[ivars] = add_time_bounds(lev1["time"][index], params["scan_time"])
            elif (ivars == "time_bnds") & (data_type in ("2P04", "2P07", "2P08")):
                rpg_dat[ivars] = np.ones(lev1[ivars].shape, np.int32) * Fill_Value_Int
            else:
                rpg_dat[ivars] = lev1[ivars][index]


def _add_att(global_attributes: dict, coeff: dict, lev1: dict) -> None:
    "add retrieval and calibration attributes"
    fields = [
        "retrieval_type",
        "retrieval_elevation_angles",
        "retrieval_frequencies",
        "retrieval_auxiliary_input",
        "retrieval_description",
    ]
    for name in fields:
        if name in coeff:
            global_attributes[name] = coeff[name]
        else:
            global_attributes[name] = ""

    fields = [
        "instrument_calibration_status",
        "receiver1_date_of_last_absolute_calibration",
        "receiver1_type_of_last_absolute_calibration",
        "receiver2_date_of_last_absolute_calibration",
        "receiver2_type_of_last_absolute_calibration",
    ]
    for name in fields:
        if hasattr(lev1, name):
            global_attributes[name] = eval("lev1." + name)
        else:
            global_attributes[name] = ""


def load_product(site: str, date_in: str, prod: str):
    "load existing lev2 file for deriving other products"
    file = []
    global_attributes, params = get_site_specs(site, "1C01")
    ID = global_attributes["wigos_station_id"]
    data_out_l2 = params["data_out"] + "level2/" + date_in.strftime("%Y/%m/%d/")
    file_name = data_out_l2 + "MWR_" + prod + "_" + ID + "_" + date_in.strftime("%Y%m%d") + ".nc"

    if os.path.isfile(file_name):
        file = nc.Dataset(file_name)
    # Load single pointing T if no BL scans are performed
    if (not os.path.isfile(file_name)) & (prod == "2P02"):
        file_name = data_out_l2 + "MWR_2P01_" + ID + "_" + date_in.strftime("%Y%m%d") + ".nc"
        if os.path.isfile(file_name):
            file = nc.Dataset(file_name)
            prod = "2P01"
    ret_freq = get_ret_freq(file_name)
    ret_ang = get_ret_ang(file_name)
    return file, ret_freq, ret_ang, prod


def _test_BL_scan(site: str, lev1: dict) -> bool:
    "Check for existing BL scans in lev1 data"
    BL_scan = True
    coeff = get_mvr_coeff(site, "tel", lev1["frequency"][:])
    ix0 = np.where(
        (lev1["ele"][:] > coeff[0]["ele"][0] - 0.5)
        & (lev1["ele"][:] < coeff[0]["ele"][0] + 0.5)
        & (lev1["pointing_flag"][:] == 1)
    )[0]
    if len(ix0) == 0:
        BL_scan = False
    return BL_scan


def ele_ret(ele_arr: np.ndarray, coeff: dict) -> np.ndarray:
    "Extracts elevation angles used in retrieval"
    val = coeff["ele"]
    if val.shape == ():
        val = np.array([val])
    return np.array([val[(np.abs(val - v)).argmin()] for v in ele_arr])


def retrieval_input(lev1: dict, cf: list) -> np.ndarray:
    "Get retrieval input"
    bias = np.ones((len(lev1["time"][:]), 1), np.float32)
    doy = np.ones((len(lev1["time"][:]), 2), np.float32) * Fill_Value_Float
    tz = tzwhere.tzwhere()
    sun = np.ones((len(lev1["time"][:]), 2), np.float32) * Fill_Value_Float
    irt = np.ones((len(lev1["time"][:]), 2), np.float32) * Fill_Value_Float
    for i, time in enumerate(lev1["time"][:].data):
        timezone_str = tz.tzNameAt(lev1["station_latitude"][i], lev1["station_longitude"][i])
        timezone = pytz.timezone(timezone_str)
        dt = datetime.fromtimestamp(time, timezone)
        dy = datetime(dt.year, 12, 31, 0, 0).timetuple().tm_yday
        doy[i, 0] = np.cos(datetime.fromtimestamp(time).timetuple().tm_yday / dy * 2 * np.pi)
        doy[i, 1] = np.sin(datetime.fromtimestamp(time).timetuple().tm_yday / dy * 2 * np.pi)
        sun[i, 0] = np.cos((dt.hour + dt.minute / 60 + dt.second / 3600) / 24 * 2 * np.pi)
        sun[i, 1] = np.sin((dt.hour + dt.minute / 60 + dt.second / 3600) / 24 * 2 * np.pi)

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
                        np.reshape(lev1["air_temperature"][:].data, (len(lev1["time"][:]), 1)),
                    ),
                    axis=1,
                )
            if (field == "HS") & (cf["aux_flg"][i] == 1):
                ret_in = np.concatenate(
                    (
                        ret_in,
                        np.reshape(lev1["relative_humidity"][:].data, (len(lev1["time"][:]), 1)),
                    ),
                    axis=1,
                )
            if (field == "PS") & (cf["aux_flg"][i] == 1):
                ret_in = np.concatenate(
                    (
                        ret_in,
                        np.reshape(lev1["air_pressure"][:].data * 100., (len(lev1["time"][:]), 1)),
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
