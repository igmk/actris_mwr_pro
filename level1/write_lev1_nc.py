from read_specs import get_site_specs
from level1.rpg_bin import get_rpg_bin
from level1.lev1_meta_nc import get_data_attributes
from level1.quality_control import apply_qc
from level1.met_quality_control import apply_met_qc

# from level1.tb_offset import correct_tb_offset
from utils import isbit, epoch2unix, add_time_bounds
import rpg_mwr
import numpy as np
import pandas as pd
import glob
import datetime
from itertools import groupby

Fill_Value_Float = -999.0
Fill_Value_Int = -99


def lev1_to_nc(
    site: str,
    data_type: str,
    path_to_files: str,
    path_to_prev: str,
    path_to_next: str,
    output_file: str,
):
    """This function reads one day of RPG MWR binary files,
    adds attributes and writes it into netCDF file.

    Args:
        site: Name of site.
        data_type: Data type of the netCDF file.
        path_to_files: Folder containing one day of RPG MWR binary files.
        path_to_prev: Folder containing previous day RPG MWR binary files.
        path_to_next: Folder containing next day RPG MWR binary files.
        output_file: Output file name.

    Examples:
        >>> from level1.write_lev1_nc import lev1_to_nc
        >>> lev1_to_nc('site_name', '1B01', '/path/to/files/', '/path/to/previous/day/', '/path/to/next/day/', 'rpg_mwr.nc')
    """

    file_list_hkd = get_file_list(path_to_files, path_to_prev, path_to_next, "hkd")
    global_attributes, params = get_site_specs(site, data_type)
    rpg_bin = prepare_data(path_to_files, path_to_prev, path_to_next, data_type, params)
    if data_type in ("1B01", "1C01"):
        apply_qc(site, rpg_bin.data, params)
        rpg_cal = cal_his(params, global_attributes, rpg_bin.data["time"][0])
        # if rpg_cal != []:
        #     correct_tb_offset(site, rpg_bin.data, params, global_attributes, rpg_cal)
    if data_type in ("1B21", "1C01"):
        apply_met_qc(rpg_bin.data, params)
    hatpro = rpg_mwr.Rpg(rpg_bin.data)
    hatpro.find_valid_times()
    hatpro.data = get_data_attributes(hatpro.data, data_type)
    rpg_mwr.save_rpg(hatpro, output_file, global_attributes, data_type, params)


def get_file_list(path_to_files: str, path_to_prev: str, path_to_next: str, extension: str):

    f_list = sorted(glob.glob(path_to_files + "*." + extension))
    if len(f_list) == 0:
        f_list = sorted(glob.glob(path_to_files + "*." + extension.upper()))
    f_list_p = sorted(glob.glob(path_to_prev + "*." + extension))
    if len(f_list_p) == 0:
        f_list_p = sorted(glob.glob(path_to_prev + "*." + extension.upper()))
    f_list_n = sorted(glob.glob(path_to_next + "*." + extension))
    if len(f_list_n) == 0:
        f_list_n = sorted(glob.glob(path_to_next + "*." + extension.upper()))

    if len(f_list) == 0:
        raise RuntimeError(
            [
                "Error: no binary files with extension "
                + extension
                + " found in directory "
                + path_to_files
            ]
        )
    if (len(f_list_p) > 0) & (len(f_list_n) > 0):
        f_list = [f_list_p[-1], *f_list, f_list_n[0]]
    elif (len(f_list_p) > 0) & (len(f_list_n) == 0):
        f_list = [f_list_p[-1], *f_list]
    elif (len(f_list_p) == 0) & (len(f_list_n) > 0):
        f_list = [*f_list, f_list_n[0]]
    return f_list


def prepare_data(
    path_to_files: str,
    path_to_prev: str,
    path_to_next: str,
    data_type: str,
    params: dict,
) -> dict:
    """Load and prepare data for netCDF writing"""

    if data_type in ("1B01", "1C01"):

        file_list_brt = get_file_list(path_to_files, path_to_prev, path_to_next, "brt")
        rpg_bin = get_rpg_bin(file_list_brt)
        rpg_bin.data["frequency"] = rpg_bin.header["_f"]
        fields = [
            "bandwidth",
            "n_sidebands",
            "sideband_IF_separation",
            "freq_shift",
            "receiver_nb",
            "receiver",
        ]
        for name in fields:
            rpg_bin.data[name] = params[name]
        rpg_bin.data["time_bnds"] = add_time_bounds(rpg_bin.data["time"], params["int_time"])
        rpg_bin.data["pointing_flag"] = np.zeros(len(rpg_bin.data["time"]), np.int32)
        (
            rpg_bin.data["liquid_cloud_flag"],
            rpg_bin.data["liquid_cloud_flag_status"],
        ) = find_lwcl_free(rpg_bin.data, np.arange(len(rpg_bin.data["time"])))

        file_list_hkd = get_file_list(path_to_files, path_to_prev, path_to_next, "hkd")
        rpg_hkd = get_rpg_bin(file_list_hkd)
        # try:
        file_list_blb = get_file_list(path_to_files, path_to_prev, path_to_next, "blb")
        rpg_blb = get_rpg_bin(file_list_blb)
        _add_blb(rpg_bin, rpg_blb, rpg_hkd, params)
        # except:
        #     print(["No binary files with extension blb found in directory " + path_to_files])

        if params["azi_cor"] != Fill_Value_Float:
            _azi_correction(rpg_bin.data, params)

        if data_type == "1C01":
            try:
                file_list_irt = get_file_list(path_to_files, path_to_prev, path_to_next, "irt")
                rpg_irt = get_rpg_bin(file_list_irt)
                rpg_irt.data["irt"][rpg_irt.data["irt"] < 150.0] = Fill_Value_Float
                rpg_bin.data["ir_wavelength"] = rpg_irt.header["_f"]
                rpg_bin.data["ir_bandwidth"] = params["ir_bandwidth"]
                rpg_bin.data["ir_beamwidth"] = params["ir_beamwidth"]
                _add_interpol1(rpg_bin.data, rpg_irt.data["irt"], rpg_irt.data["time"], "irt")
                _add_interpol1(rpg_bin.data, rpg_irt.data["ir_ele"], rpg_irt.data["time"], "ir_ele")
                _add_interpol1(rpg_bin.data, rpg_irt.data["ir_azi"], rpg_irt.data["time"], "ir_azi")
            except:
                print(["No binary files with extension irt found in directory " + path_to_files])

            try:
                file_list_met = get_file_list(path_to_files, path_to_prev, path_to_next, "met")
                rpg_met = get_rpg_bin(file_list_met)
                _add_interpol1(
                    rpg_bin.data,
                    rpg_met.data["air_temperature"],
                    rpg_met.data["time"],
                    "air_temperature",
                )
                _add_interpol1(
                    rpg_bin.data,
                    rpg_met.data["relative_humidity"],
                    rpg_met.data["time"],
                    "relative_humidity",
                )
                _add_interpol1(
                    rpg_bin.data,
                    rpg_met.data["air_pressure"],
                    rpg_met.data["time"],
                    "air_pressure",
                )
                if (int(rpg_met.header["_n_sen"], 2) & 1) != 0:
                    _add_interpol1(
                        rpg_bin.data,
                        rpg_met.data["adds"][:, 0],
                        rpg_met.data["time"],
                        "wind_speed",
                    )
                    rpg_bin.data["wind_speed"] = rpg_bin.data["wind_speed"] / 3.6
                if (int(rpg_met.header["_n_sen"], 2) & 2) != 0:
                    _add_interpol1(
                        rpg_bin.data,
                        rpg_met.data["adds"][:, 1],
                        rpg_met.data["time"],
                        "wind_direction",
                    )
                if (int(rpg_met.header["_n_sen"], 2) & 4) != 0:
                    _add_interpol1(
                        rpg_bin.data,
                        rpg_met.data["adds"][:, 2],
                        rpg_met.data["time"],
                        "rain_rate",
                    )
            except:
                print(["No binary files with extension met found in directory " + path_to_files])

    elif data_type == "1B11":
        file_list_irt = get_file_list(path_to_files, path_to_prev, path_to_next, "irt")
        rpg_bin = get_rpg_bin(file_list_irt)
        rpg_bin.data["ir_wavelength"] = rpg_bin.header["_f"]
        rpg_bin.data["ir_bandwidth"] = params["ir_bandwidth"]
        rpg_bin.data["ir_beamwidth"] = params["ir_beamwidth"]

    elif data_type == "1B21":
        file_list_met = get_file_list(path_to_files, path_to_prev, path_to_next, "met")
        rpg_bin = get_rpg_bin(file_list_met)
        if (int(rpg_bin.header["_n_sen"], 2) & 1) != 0:
            rpg_bin.data["wind_speed"] = rpg_bin.data["adds"][:, 0] / 3.6
        if (int(rpg_bin.header["_n_sen"], 2) & 2) != 0:
            rpg_bin.data["wind_direction"] = rpg_bin.data["adds"][:, 1]
        if (int(rpg_bin.header["_n_sen"], 2) & 4) != 0:
            rpg_bin.data["rain_rate"] = rpg_bin.data["adds"][:, 2]

    else:
        raise RuntimeError(["Data type " + data_type + " not supported for file writing."])

    file_list_hkd = get_file_list(path_to_files, path_to_prev, path_to_next, "hkd")
    _append_hkd(file_list_hkd, rpg_bin, data_type, params)
    rpg_bin.data["station_altitude"] = (
        np.ones(len(rpg_bin.data["time"]), np.float32) * params["station_altitude"]
    )

    return rpg_bin


def _append_hkd(file_list_hkd: list, rpg_bin: dict, data_type: str, params: dict) -> None:
    """Append hkd data on same time grid and perform TB sanity check"""

    hkd = get_rpg_bin(file_list_hkd)

    if all(hkd.data["station_latitude"] == Fill_Value_Float):
        _add_interpol1(
            rpg_bin.data,
            np.ones(len(hkd.data["time"])) * params["station_latitude"],
            hkd.data["time"],
            "station_latitude",
        )
    else:
        idx = np.where(hkd.data["station_latitude"] != Fill_Value_Float)[0]
        _add_interpol1(
            rpg_bin.data,
            hkd.data["station_latitude"][idx],
            hkd.data["time"][idx],
            "station_latitude",
        )
    if all(hkd.data["station_longitude"] == Fill_Value_Float):
        _add_interpol1(
            rpg_bin.data,
            np.ones(len(hkd.data["time"])) * params["station_longitude"],
            hkd.data["time"],
            "station_longitude",
        )
    else:
        idx = np.where(hkd.data["station_longitude"] != Fill_Value_Float)[0]
        _add_interpol1(
            rpg_bin.data,
            hkd.data["station_longitude"][idx],
            hkd.data["time"][idx],
            "station_longitude",
        )

    if data_type in ("1B01", "1C01"):
        _add_interpol1(rpg_bin.data, hkd.data["temp"][:, 0:2], hkd.data["time"], "t_amb")
        _add_interpol1(rpg_bin.data, hkd.data["temp"][:, 2:4], hkd.data["time"], "t_rec")
        _add_interpol1(rpg_bin.data, hkd.data["stab"], hkd.data["time"], "t_sta")

        """check time intervals of +-15 min of .HKD data for sanity checks"""
        rpg_bin.data["status"] = np.zeros(rpg_bin.data["tb"].shape, np.int32)
        for i_time, v_time in enumerate(rpg_bin.data["time"]):
            ind = np.where((hkd.data["time"] >= v_time - 900) & (hkd.data["time"] <= v_time + 900))
            status = hkd.data["status"][ind]
            for bit in range(7):
                # status flags for channel 1 to 7 of the humidity profiler receiver
                if np.any(~isbit(status, bit)):
                    rpg_bin.data["status"][i_time, bit] = 1
                # status flags for channel 1 to 7 of the temperature profiler receiver
                if np.any(~isbit(status, bit + 8)):
                    rpg_bin.data["status"][i_time, bit + 7] = 1
            # receiver 1 (humidity profiler) thermal stability & ambient target stability & noise diode
            if np.any(
                isbit(status, 25)
                | isbit(status, 29)
                | ~isbit(status, 22)
                | (~isbit(status, 24) & ~isbit(status, 25))
            ):
                rpg_bin.data["status"][i_time, 0:6] = 1
            # receiver 2 (temperature profiler) thermal stability & ambient target stability & noise diode
            if np.any(
                isbit(status, 27)
                | isbit(status, 29)
                | ~isbit(status, 23)
                | (~isbit(status, 26) & ~isbit(status, 27))
            ):
                rpg_bin.data["status"][i_time, 7:13] = 1


def _add_interpol1(data0: dict, data1: np.ndarray, time1: np.ndarray, output_name: str) -> None:

    if data1.ndim > 1:
        data0[output_name] = (
            np.ones([len(data0["time"]), data1.shape[1]], np.float32) * Fill_Value_Float
        )
        for ndim in range(data1.shape[1]):
            data0[output_name][:, ndim] = np.interp(data0["time"], time1, data1[:, ndim])
    else:
        data0[output_name] = np.interp(data0["time"], time1, data1)


def find_lwcl_free(lev1: dict, ix: np.ndarray) -> tuple:

    index = np.ones(len(ix)) * np.nan
    status = np.ones(len(ix), dtype=np.int32)
    freq_31 = np.where(np.round(lev1["frequency"][:], 1) == 31.4)[0]
    if len(freq_31) == 1:
        time = lev1["time"][ix]
        tb = np.squeeze(lev1["tb"][ix, freq_31])
        tb[(lev1["pointing_flag"][ix] == 1) | (lev1["ele"][ix] < 89.0)] = np.nan
        tb_df = pd.DataFrame({"Tb": tb}, index=pd.to_datetime(time, unit="s"))
        tb_std = tb_df.rolling("2min", center=True, min_periods=10).std()
        tb_mx = tb_std.rolling("20min", center=True, min_periods=100).max()

        if "irt" in lev1:
            irt = lev1["irt"][ix, 0]
            irt[(lev1["pointing_flag"][ix] == 1) | (lev1["ele"][ix] < 89.0)] = np.nan
            irt_df = pd.DataFrame({"Irt": irt[:]}, index=pd.to_datetime(time, unit="s"))
            irt_mx = irt_df.rolling("20min", center=True, min_periods=100).max()
            index[(irt_mx["Irt"] > 233.15) & (tb_mx["Tb"] > 0.3)] = 1
            status[:] = 0
        else:
            index[(tb_mx["Tb"] > 0.3)] = 1

        df = pd.DataFrame({"index": index}, index=pd.to_datetime(time, unit="s"))
        df = df.fillna(method="bfill", limit=120)
        df = df.fillna(method="ffill", limit=120)
        index = np.array(df["index"])
        index[(tb_mx["Tb"] < 0.3) & (index != 1.0)] = 0.0
        index[(lev1["ele"][ix] < 89.0) & (index != 1.0)] = 2.0

    return np.nan_to_num(index, nan=2).astype(int), status


def _add_blb(brt: dict, blb: dict, hkd: dict, params: dict) -> None:
    """Add boundary-layer scans using a linear time axis"""

    time_add, ele_add, azi_add, rain_add, tb_add = (
        np.empty([0], dtype=np.int32),
        [],
        [],
        [],
        [],
    )
    seqs = [(key, len(list(val))) for key, val in groupby(hkd.data["status"][:] & 2**18 > 0)]
    seqs = np.array(
        [
            (key, sum(s[1] for s in seqs[:i]), len)
            for i, (key, len) in enumerate(seqs)
            if key == True
        ]
    )

    for time_ind, time_blb in enumerate(blb.data["time"]):
        seqi = np.where(np.abs(hkd.data["time"][seqs[:, 1] + seqs[:, 2] - 1] - time_blb) < 2)[0]

        if len(seqi) == 1:
            sq = 0.0  # scan quadrant, 0: 1st, 180: 2nd
            if (not isbit(blb.data["rf_mod"][time_ind], 5)) & (
                isbit(blb.data["rf_mod"][time_ind], 6)
            ):
                sq = 180.0

            time_add = np.concatenate(
                (
                    time_add,
                    np.squeeze(
                        np.linspace(
                            hkd.data["time"][seqs[seqi, 1]]
                            + int(
                                np.floor(
                                    ((time_blb - 2) - hkd.data["time"][seqs[seqi, 1]])
                                    / (blb.header["_n_ang"] - 1)
                                )
                            ),
                            hkd.data["time"][seqs[seqi, 1]]
                            + (blb.header["_n_ang"] - 1)
                            * int(
                                np.floor(
                                    ((time_blb - 2) - hkd.data["time"][seqs[seqi, 1]])
                                    / (blb.header["_n_ang"] - 1)
                                )
                            ),
                            blb.header["_n_ang"] - 1,
                            dtype=np.int32,
                        )
                    ),
                )
            )
            time_add = np.append(time_add, time_blb - 1)
            # time_add = np.append(time_add, hkd.data['time'][seqs[seqi, 1] + seqs[seqi, 2]])

            azi_add = np.concatenate(
                (
                    azi_add,
                    np.ones(blb.header["_n_ang"]) * ((sq + params["const_azi"]) % 360),
                )
            )
            rain_add = np.concatenate(
                (
                    rain_add,
                    np.ones(blb.header["_n_ang"], np.int32)
                    * int(isbit(blb.data["rf_mod"][time_ind], 1)),
                )
            )
            ele_add = np.concatenate((ele_add, blb.header["_ang"]))
            for ang in range(blb.header["_n_ang"]):
                if len(tb_add) == 0:
                    tb_add = blb.data["tb"][time_ind, :, ang]
                else:
                    tb_add = np.vstack((tb_add, blb.data["tb"][time_ind, :, ang]))

    time_bnds_add = add_time_bounds(
        time_add[0:-1], int(np.floor(params["scan_time"] / (blb.header["_n_ang"] - 1)))
    )
    time_bnds_add = np.concatenate(
        (time_bnds_add, add_time_bounds(np.array(time_add[-1], ndmin=1), params["int_time"]))
    )
    pointing_flag_add = np.ones(len(time_add), np.int32)
    liquid_cloud_flag_add = np.ones(len(time_add), np.int32) * 2
    liquid_cloud_flag_status_add = np.ones(len(time_add), np.int32) * Fill_Value_Int
    brt.data["time"] = np.concatenate((brt.data["time"], time_add))
    ind = np.argsort(brt.data["time"])
    brt.data["time"] = brt.data["time"][ind]
    names = [
        "time_bnds",
        "ele",
        "azi",
        "rain",
        "tb",
        "pointing_flag",
        "liquid_cloud_flag",
        "liquid_cloud_flag_status",
    ]
    for var in names:
        brt.data[var] = np.concatenate((brt.data[var], eval(var + "_add")))
        if brt.data[var].ndim > 1:
            brt.data[var] = brt.data[var][ind, :]
        else:
            brt.data[var] = brt.data[var][ind]
    brt.header["n"] = brt.header["n"] + blb.header["n"] * blb.header["_n_ang"]


def _azi_correction(brt: dict, params: dict) -> None:
    """Azimuth correction (transform to "geographical" coordinates)"""

    ind180 = np.where((brt["azi"][:] >= 0) & (brt["azi"][:] <= 180))
    ind360 = np.where((brt["azi"][:] > 180) & (brt["azi"][:] <= 360))
    brt["azi"][ind180] = params["azi_cor"] - brt["azi"][ind180]
    brt["azi"][ind360] = 360.0 + params["azi_cor"] - brt["azi"][ind360]
    brt["azi"][brt["azi"][:] < 0] += 360.0


def cal_his(params: dict, glob_att: dict, time0: int) -> dict:
    """Load and add information from ABSCAL.HIS file"""
    file_list_cal, rpg_cal = [], []
    cal_type = [
        "instrument performs no absolute calibration",
        "liquid nitrogen calibration",
        "sky tipping calibration",
    ]
    try:
        file_list_cal = get_file_list(params["path_to_cal"], " ", " ", "his")
    except:
        print(["No binary files with extension his found in directory " + params["path_to_cal"]])

    if file_list_cal != []:
        rpg_cal = get_rpg_bin(file_list_cal)

        if (
            (rpg_cal.data["cal1_t"] > 0)
            & (rpg_cal.data["cal2_t"] > 0)
            & (
                (time0 - epoch2unix(rpg_cal.data["t1"], rpg_cal.header["_time_ref"])) / 86400.0
                < 183.0
            )
        ):
            glob_att["instrument_calibration_status"] = "calibrated"
        else:
            glob_att["instrument_calibration_status"] = "needs calibration"

        cal_t1 = np.where(epoch2unix(rpg_cal.data["t1"], rpg_cal.header["_time_ref"]) < time0)[0]
        if len(cal_t1) > 0:
            glob_att[
                "receiver1_date_of_last_absolute_calibration"
            ] = datetime.datetime.utcfromtimestamp(
                epoch2unix(rpg_cal.data["t1"][cal_t1[-1]], rpg_cal.header["_time_ref"])
            ).strftime(
                "%Y%m%d"
            )
            glob_att["receiver1_type_of_last_absolute_calibration"] = cal_type[
                int(rpg_cal.data["cal1_t"])
            ]
        cal_t2 = np.where(epoch2unix(rpg_cal.data["t2"], rpg_cal.header["_time_ref"]) < time0)[0]
        if len(cal_t2) > 0:
            glob_att[
                "receiver2_date_of_last_absolute_calibration"
            ] = datetime.datetime.utcfromtimestamp(
                epoch2unix(rpg_cal.data["t2"][cal_t2[-1]], rpg_cal.header["_time_ref"])
            ).strftime(
                "%Y%m%d"
            )
            glob_att["receiver2_type_of_last_absolute_calibration"] = cal_type[
                int(rpg_cal.data["cal2_t"])
            ]

    return rpg_cal
