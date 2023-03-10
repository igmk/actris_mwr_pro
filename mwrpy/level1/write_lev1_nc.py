"""Module for writing Level 1 netCDF files"""
import datetime
from itertools import groupby

import numpy as np
import pandas as pd
from numpy import ma

import rpg_mwr
from level1.lev1_meta_nc import get_data_attributes
from level1.met_quality_control import apply_met_qc
from level1.quality_control import apply_qc
from level1.rpg_bin import get_rpg_bin

from utils import (
    add_interpol1d,
    add_time_bounds,
    epoch2unix,
    get_file_list,
    isbit,
    read_yaml_config,
    update_lev1_attributes,
)

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
        >>> lev1_to_nc('site_name', '1B01', '/path/to/files/', '/path/to/previous/day/',
        '/path/to/next/day/', 'rpg_mwr.nc')
    """

    global_attributes, params = read_yaml_config(site)
    if data_type != "1C01":
        update_lev1_attributes(global_attributes, data_type)
    rpg_bin = prepare_data(path_to_files, path_to_prev, path_to_next, data_type, params, site)
    if data_type in ("1B01", "1C01"):
        apply_qc(site, rpg_bin.data, params)
    if data_type in ("1B21", "1C01"):
        apply_met_qc(rpg_bin.data, params)
    hatpro = rpg_mwr.Rpg(rpg_bin.data)
    hatpro.find_valid_times()
    hatpro.data = get_data_attributes(hatpro.data, data_type)
    rpg_mwr.save_rpg(hatpro, output_file, global_attributes, data_type)


def prepare_data(
    path_to_files: str,
    path_to_prev: str,
    path_to_next: str,
    data_type: str,
    params: dict,
    site: str,
) -> dict:
    """Load and prepare data for netCDF writing"""

    if data_type in ("1B01", "1C01"):

        file_list_brt = get_file_list(path_to_files, path_to_prev, path_to_next, "brt")
        rpg_bin = get_rpg_bin(file_list_brt)
        rpg_bin.data["tb"] = rpg_bin.data["tb"][:, np.argsort(params["bandwidth"])]
        rpg_bin.data["frequency"] = rpg_bin.header["_f"][np.argsort(params["bandwidth"])]
        fields = [
            "bandwidth",
            "n_sidebands",
            "sideband_IF_separation",
            "freq_shift",
            "receiver_nb",
            "receiver",
        ]
        for name in fields:
            rpg_bin.data[name] = np.array(params[name])
        rpg_bin.data["time_bnds"] = add_time_bounds(rpg_bin.data["time"], params["int_time"])
        rpg_bin.data["pointing_flag"] = np.zeros(len(rpg_bin.data["time"]), np.int32)

        if data_type == "1B01":
            (
                rpg_bin.data["liquid_cloud_flag"],
                rpg_bin.data["liquid_cloud_flag_status"],
            ) = find_lwcl_free(rpg_bin.data, np.arange(len(rpg_bin.data["time"])))
        else:
            (
                rpg_bin.data["liquid_cloud_flag"],
                rpg_bin.data["liquid_cloud_flag_status"],
            ) = np.ones(len(rpg_bin.data["time"]), np.int32) * 2, np.ones(
                len(rpg_bin.data["time"]), np.int32
            )

        file_list_hkd = get_file_list(path_to_files, path_to_prev, path_to_next, "hkd")
        rpg_hkd = get_rpg_bin(file_list_hkd)
        rpg_bin.data["status"] = np.zeros(
            (len(rpg_bin.data["time"]), len(params["receiver"])), np.int32
        )
        _, ind_brt, ind_hkd = np.intersect1d(
            rpg_bin.data["time"], rpg_hkd.data["time"], assume_unique=False, return_indices=True
        )
        rpg_bin.data["status"][ind_brt, :] = hkd_sanity_check(
            rpg_hkd.data["status"][ind_hkd], params
        )
        if params["scan_time"] != Fill_Value_Int:
            file_list_bls = []
            try:
                file_list_bls = get_file_list(path_to_files, path_to_prev, path_to_next, "bls")
            except:
                print(["No binary files with extension bls found in directory " + path_to_files])
            if len(file_list_bls) > 0:
                rpg_bls = get_rpg_bin(file_list_bls)
                _add_bls(rpg_bin, rpg_bls, rpg_hkd, params, site)
            else:
                try:
                    file_list_blb = get_file_list(path_to_files, path_to_prev, path_to_next, "blb")
                    rpg_blb = get_rpg_bin(file_list_blb)
                    _add_blb(rpg_bin, rpg_blb, rpg_hkd, params, site)
                except:
                    print(
                        ["No binary files with extension blb found in directory " + path_to_files]
                    )

        if params["azi_cor"] != Fill_Value_Float:
            _azi_correction(rpg_bin.data, params)
        if params["const_azi"] != Fill_Value_Float:
            rpg_bin.data["azi"] = (rpg_bin.data["azi"] + params["const_azi"]) % 360

        if data_type == "1C01":
            if params["ir_flag"]:
                try:
                    file_list_irt = get_file_list(path_to_files, path_to_prev, path_to_next, "irt")
                    rpg_irt = get_rpg_bin(file_list_irt)
                    rpg_irt.data["irt"][rpg_irt.data["irt"] <= 125.5] = Fill_Value_Float
                    rpg_bin.data["ir_wavelength"] = rpg_irt.header["_f"]
                    rpg_bin.data["ir_bandwidth"] = params["ir_bandwidth"]
                    rpg_bin.data["ir_beamwidth"] = params["ir_beamwidth"]
                    add_interpol1d(rpg_bin.data, rpg_irt.data["irt"], rpg_irt.data["time"], "irt")
                    add_interpol1d(
                        rpg_bin.data, rpg_irt.data["ir_ele"], rpg_irt.data["time"], "ir_ele"
                    )
                    add_interpol1d(
                        rpg_bin.data, rpg_irt.data["ir_azi"], rpg_irt.data["time"], "ir_azi"
                    )
                except:
                    print(
                        ["No binary files with extension irt found in directory " + path_to_files]
                    )

            (
                rpg_bin.data["liquid_cloud_flag"],
                rpg_bin.data["liquid_cloud_flag_status"],
            ) = find_lwcl_free(rpg_bin.data, np.arange(len(rpg_bin.data["time"])))

            try:
                file_list_met = get_file_list(path_to_files, path_to_prev, path_to_next, "met")
                rpg_met = get_rpg_bin(file_list_met)
                add_interpol1d(
                    rpg_bin.data,
                    rpg_met.data["air_temperature"],
                    rpg_met.data["time"],
                    "air_temperature",
                )
                add_interpol1d(
                    rpg_bin.data,
                    rpg_met.data["relative_humidity"],
                    rpg_met.data["time"],
                    "relative_humidity",
                )
                add_interpol1d(
                    rpg_bin.data,
                    rpg_met.data["air_pressure"],
                    rpg_met.data["time"],
                    "air_pressure",
                )
                if (int(rpg_met.header["_n_sen"], 2) & 1) != 0:
                    add_interpol1d(
                        rpg_bin.data,
                        rpg_met.data["adds"][:, 0],
                        rpg_met.data["time"],
                        "wind_speed",
                    )
                    rpg_bin.data["wind_speed"] = rpg_bin.data["wind_speed"] / 3.6
                if (int(rpg_met.header["_n_sen"], 2) & 2) != 0:
                    add_interpol1d(
                        rpg_bin.data,
                        rpg_met.data["adds"][:, 1],
                        rpg_met.data["time"],
                        "wind_direction",
                    )
                if (int(rpg_met.header["_n_sen"], 2) & 4) != 0:
                    add_interpol1d(
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
        add_interpol1d(
            rpg_bin.data,
            np.ones(len(hkd.data["time"])) * params["station_latitude"],
            hkd.data["time"],
            "station_latitude",
        )
    else:
        idx = np.where(hkd.data["station_latitude"] != Fill_Value_Float)[0]
        add_interpol1d(
            rpg_bin.data,
            hkd.data["station_latitude"][idx],
            hkd.data["time"][idx],
            "station_latitude",
        )
    if all(hkd.data["station_longitude"] == Fill_Value_Float):
        add_interpol1d(
            rpg_bin.data,
            np.ones(len(hkd.data["time"])) * params["station_longitude"],
            hkd.data["time"],
            "station_longitude",
        )
    else:
        idx = np.where(hkd.data["station_longitude"] != Fill_Value_Float)[0]
        add_interpol1d(
            rpg_bin.data,
            hkd.data["station_longitude"][idx],
            hkd.data["time"][idx],
            "station_longitude",
        )

    if data_type in ("1B01", "1C01"):
        add_interpol1d(rpg_bin.data, hkd.data["temp"][:, 0:2], hkd.data["time"], "t_amb")
        add_interpol1d(rpg_bin.data, hkd.data["temp"][:, 2:4], hkd.data["time"], "t_rec")
        add_interpol1d(rpg_bin.data, hkd.data["stab"], hkd.data["time"], "t_sta")


def hkd_sanity_check(status: np.ndarray, params: dict) -> np.ndarray:
    """Perform sanity checks for .HKD data"""
    status_flag = np.zeros((len(status), len(params["receiver"])), np.int32)
    for irec, nrec in enumerate(np.array(params["receiver"])):
        # status flags for individual channels
        status_flag[~isbit(status, irec + (nrec - 1) * (8 - irec)), irec] = 1
        if nrec == 1:
            # receiver 1 thermal stability & ambient target stability & noise diode
            status_flag[
                isbit(status, 25)
                | isbit(status, 29)
                | ~isbit(status, 22)
                | (~isbit(status, 24) & ~isbit(status, 25)),
                irec,
            ] = 1
        if nrec == 2:
            # receiver 2 thermal stability & ambient target stability & noise diode
            status_flag[
                isbit(status, 27)
                | isbit(status, 29)
                | ~isbit(status, 23)
                | (~isbit(status, 26) & ~isbit(status, 27)),
                irec,
            ] = 1

    return status_flag


def find_lwcl_free(lev1: dict, ix: np.ndarray) -> tuple:
    """Identification of liquid water cloud free periods using TB variability at 31.4 GHz and IRT.
    Uses pre-defined time index and additionally returns status of IRT availability"""

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
            tb_thres = 0.1
            irt = lev1["irt"][ix, :]
            irt[irt == Fill_Value_Float] = np.nan
            irt = np.nanmean(irt, axis=1)
            irt[(lev1["pointing_flag"][ix] == 1) | (lev1["ele"][ix] < 89.0)] = np.nan
            irt_df = pd.DataFrame({"Irt": irt[:]}, index=pd.to_datetime(time, unit="s"))
            irt_mx = irt_df.rolling("20min", center=True, min_periods=100).max()
            index[(irt_mx["Irt"] > 263.15) & (tb_mx["Tb"] > tb_thres)] = 1
            status[:] = 0

        tb_thres = 0.2
        index[(tb_mx["Tb"] > tb_thres)] = 1
        df = pd.DataFrame({"index": index}, index=pd.to_datetime(time, unit="s"))
        df = df.fillna(method="bfill", limit=120)
        df = df.fillna(method="ffill", limit=120)
        index = np.array(df["index"])
        index[(tb_mx["Tb"] < tb_thres) & (index != 1.0)] = 0.0
        index[(lev1["ele"][ix] < 89.0) & (index != 1.0)] = 2.0

    return np.nan_to_num(index, nan=2).astype(int), status


def _add_bls(brt: dict, bls: dict, hkd: dict, params: dict, site: str) -> None:
    """Add BLS boundary-layer scans using a linear time axis"""

    bls.data["time_bnds"] = add_time_bounds(bls.data["time"] + 1, params["int_time"])
    bls.data["status"] = np.zeros((len(bls.data["time"]), len(params["receiver"])), np.int32)

    for time_ind, time_bls in enumerate(bls.data["time"]):
        if np.min(np.abs(hkd.data["time"] - time_bls)) <= params["int_time"] * 2:
            ind_hkd = np.argmin(np.abs(hkd.data["time"] - time_bls))
            bls.data["status"][time_ind, :] = hkd_sanity_check(
                np.array([hkd.data["status"][ind_hkd]], np.int32), params
            )

    bls.data["pointing_flag"] = np.ones(len(bls.data["time"]), np.int32)
    bls.data["liquid_cloud_flag"] = np.ones(len(bls.data["time"]), np.int32) * 2
    bls.data["liquid_cloud_flag_status"] = np.ones(len(bls.data["time"]), np.int32) * Fill_Value_Int
    brt.data["time"] = np.concatenate((brt.data["time"], bls.data["time"]))
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
        "status",
    ]
    for var in names:
        brt.data[var] = np.concatenate((brt.data[var], bls.data[var]))
        if brt.data[var].ndim > 1:
            brt.data[var] = brt.data[var][ind, :]
        else:
            brt.data[var] = brt.data[var][ind]
    brt.header["n"] = len(brt.data["time"])


def _add_blb(brt: dict, blb: dict, hkd: dict, params: dict, site: str) -> None:
    """Add BLB boundary-layer scans using a linear time axis"""

    time_add, time_bnds_add, ele_add, azi_add, rain_add, tb_add, status_add = (
        np.empty([0], dtype=np.int32),
        [],
        [],
        [],
        [],
        [],
        [],
    )
    const_azi = params["const_azi"]
    seqs = [(key, len(list(val))) for key, val in groupby(hkd.data["status"][:] & 2**18 > 0)]
    seqs = np.array(
        [
            (key, sum(s[1] for s in seqs[:i]), len)
            for i, (key, len) in enumerate(seqs)
            if key == True
        ]
    )

    for time_ind, time_blb in enumerate(blb.data["time"]):

        if (
            (site in ["juelich", "cologne"])
            & (
                datetime.datetime.utcfromtimestamp(hkd.data["time"][0])
                >= datetime.datetime(2022, 12, 1)
            )
            & (time_blb + int(params["scan_time"]) < hkd.data["time"][-1])
        ):
            time_blb = time_blb + int(params["scan_time"])
        seqi = np.where(np.abs(hkd.data["time"][seqs[:, 1] + seqs[:, 2] - 1] - time_blb) < 60)[0]
        if len(seqi) != 1:
            continue

        if np.abs(time_blb - hkd.data["time"][seqs[seqi, 1]][0]) >= blb.header["_n_ang"]:
            scan_quadrant = 0.0  # scan quadrant, 0 deg: 1st, 180 deg: 2nd
            if (isbit(blb.data["rf_mod"][time_ind], 1)) & (
                not isbit(blb.data["rf_mod"][time_ind], 2)
            ):
                scan_quadrant = 180.0

            time_add = np.concatenate(
                (
                    time_add,
                    np.squeeze(
                        np.linspace(
                            hkd.data["time"][seqs[seqi, 1]]
                            + int(
                                np.floor(
                                    ((time_blb - 1) - hkd.data["time"][seqs[seqi, 1]])
                                    / (blb.header["_n_ang"])
                                )
                            ),
                            hkd.data["time"][seqs[seqi, 1]]
                            + (blb.header["_n_ang"])
                            * int(
                                np.floor(
                                    ((time_blb - 1) - hkd.data["time"][seqs[seqi, 1]])
                                    / (blb.header["_n_ang"])
                                )
                            ),
                            blb.header["_n_ang"],
                            dtype=np.int32,
                        )
                    ),
                )
            )

            brt_ind = np.where(
                (brt.data["time"] > time_blb - 3600) & (brt.data["time"] < time_blb + 3600)
            )[0]
            brt_azi = ma.median(brt.data["azi"][brt_ind])
            azi_add = np.concatenate(
                (
                    azi_add,
                    np.ones(blb.header["_n_ang"]) * ((scan_quadrant + brt_azi) % 360),
                )
            )

            rain_add = np.concatenate(
                (
                    rain_add,
                    np.ones(blb.header["_n_ang"], np.int32)
                    * int(isbit(blb.data["rf_mod"][time_ind], 0)),
                )
            )
            ele_add = np.concatenate((ele_add, blb.header["_ang"]))

            for ang in range(blb.header["_n_ang"]):
                if len(tb_add) == 0:
                    tb_add = blb.data["tb"][time_ind, :, ang]
                else:
                    tb_add = np.vstack((tb_add, blb.data["tb"][time_ind, :, ang]))

            if len(time_bnds_add) == 0:
                time_bnds_add = add_time_bounds(
                    time_add, int(np.floor(params["scan_time"] / (blb.header["_n_ang"])))
                )
            else:
                time_bnds_add = np.concatenate(
                    (
                        time_bnds_add,
                        add_time_bounds(
                            time_add, int(np.floor(params["scan_time"] / (blb.header["_n_ang"])))
                        ),
                    )
                )

            blb_status = hkd_sanity_check(
                hkd.data["status"][seqs[seqi, 1][0] : seqs[seqi, 1][0] + seqs[seqi, 2][0]], params
            )
            blb_status_add = np.zeros((blb.header["_n_ang"], len(params["receiver"])), np.int32)
            for i_ch, _ in enumerate(params["receiver"]):
                blb_status_add[:, i_ch] = int(np.any(blb_status[:, i_ch]))
            if len(status_add) == 0:
                status_add = blb_status_add
            else:
                status_add = np.concatenate((status_add, blb_status_add))

    if len(time_add) > 0:
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
            "status",
        ]
        for var in names:
            brt.data[var] = np.concatenate((brt.data[var], eval(var + "_add")))
            if brt.data[var].ndim > 1:
                brt.data[var] = brt.data[var][ind, :]
            else:
                brt.data[var] = brt.data[var][ind]
        brt.header["n"] = len(
            brt.data["time"]
        )


def _azi_correction(brt: dict, params: dict) -> None:
    """Azimuth correction (transform to "geographical" coordinates)"""
    ind180 = np.where((brt["azi"][:] >= 0) & (brt["azi"][:] <= 180))
    ind360 = np.where((brt["azi"][:] > 180) & (brt["azi"][:] <= 360))
    brt["azi"][ind180] = params["azi_cor"] - brt["azi"][ind180]
    brt["azi"][ind360] = 360.0 + params["azi_cor"] - brt["azi"][ind360]
    brt["azi"][brt["azi"][:] < 0] += 360.0

