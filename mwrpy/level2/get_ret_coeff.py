"""Module to load in retrieval coefficient files"""
import netCDF4 as nc
import numpy as np
from numpy import ma

from utils import get_coeff_list

Fill_Value_Float = -999.0
Fill_Value_Int = -99


def get_mvr_coeff(site: str, prefix: str, freq: np.ndarray):
    """This function extracts retrieval coefficients for given files.

    Args:
        site: Name of site.
        prefix: Identifier for type of product.
        freq: Frequencies of observations.

    Examples:
        >>> from level2.get_ret_coeff import get_mvr_coeff
        >>> get_mvr_coeff('site_name', 'lwp', [22, 31.4])
    """

    c_list = get_coeff_list(site, prefix)
    coeff = {}
    coeff["ret_type"] = Fill_Value_Int

    if (str(c_list[0][-3:]).lower() == "ret") & (len(c_list) == 1):
        with open(c_list[0], "r", encoding="utf8") as f:
            lines = f.readlines()
            lines = [line.rstrip("\n") for line in lines]
            for line in lines:
                l_split = line.split("=")
                if l_split[0] == "AG":
                    N = len([float(x) for x in l_split[1].split()])

        for i_file, file in enumerate(c_list):
            with open(file, "r", encoding="utf8") as f:
                lines = f.readlines()
                lines = [line.rstrip("\n") for line in lines]

                for ind_l, line in enumerate(lines):
                    l_split = line.split("=")
                    if l_split[0] == "RT":
                        coeff["ret_type"] = int(l_split[1].split("#")[0])
                    if l_split[0] == "AG":
                        coeff["ele"] = np.ones(N) * Fill_Value_Float
                        coeff["ele"][i_file:N] = np.array(
                            [float(x) for x in l_split[1].split()], np.float32
                        )
                    if l_split[0] == "FR":
                        freq_ret = np.array(
                            [float(idx) for idx in l_split[1][:].split()], np.float32
                        )
                        _, freq_ind, freq_coeff = np.intersect1d(
                            freq[:],
                            freq_ret[:],
                            assume_unique=False,
                            return_indices=True,
                        )
                        if len(freq_coeff) < len(freq_ret):
                            raise RuntimeError(
                                ["Instrument and retrieval frequencies do not match."]
                            )
                        coeff["freq"] = np.ones([len(freq), N]) * Fill_Value_Float
                        coeff["freq"][freq_ind, i_file:N] = np.resize(
                            freq_ret[freq_coeff], (len(freq_coeff), N)
                        )
                    if l_split[0] == "AL":
                        coeff["height_grid"] = np.array(
                            [float(x) for x in l_split[1].split()], np.float32
                        )
                        coeff["n_height_grid"] = len(coeff["height_grid"])
                        
                    # Regression
                    if coeff["ret_type"] < 2:
                        coeff["aux"] = ["TS", "HS", "PS", "IS"]
                        if "aux_flg" not in coeff:
                            coeff["aux_flg"] = np.zeros(len(coeff["aux"]), np.int32)
                        for ii, aux_i in enumerate(coeff["aux"]):
                            if l_split[0] == aux_i:
                                coeff["aux_flg"][ii] = int(l_split[1])
                        if l_split[0] == "OS":
                            for jj in range(N):
                                if ":" in lines[il + jj]:
                                    cll = lines[il + jj].split(":")
                                else:
                                    cll = lines[il + jj].split("=")
                                    coeff["offset"][i + jj] = float(cll[1].split("#")[0])
                        if l_split[0] == "TL":
                            for jj in range(N):
                                if ":" in lines[il + jj]:
                                    cll = lines[il + jj].split(":")
                                else:
                                    cll = lines[il + jj].split("=")
                                coeff["coeff_lin"][i + jj, freq_ind] = np.array(
                                [float(idx) for idx in cll[1].split()[0 : len(freq_coeff)]],
                                np.float32,
                                )
                        if (l_split[0] == "TQ") & (coeff["ret_type"] == 1):
                            for jj in range(N):
                                if ":" in lines[il + jj]:
                                    cll = lines[il + jj].split(":")
                                else:
                                    cll = lines[il + jj].split("=")
                                coeff["coeff_quad"][i + jj, freq_ind] = np.array(
                                [float(idx) for idx in cll[1].split()[0 : len(freq_coeff)]],
                                np.float32,
                                )

                    # Neural Network
                    if coeff["ret_type"] == 2:
                        if l_split[0] == "ND":
                            coeff["n_hidden"] = np.array(
                                [int(idx) for idx in l_split[1][:].split()[0:2]],
                                np.int32,
                            )
                        coeff["aux"] = [
                            "TS",
                            "HS",
                            "PS",
                            "ZS",
                            "IR",
                            "I1",
                            "I2",
                            "DY",
                            "SU",
                        ]
                        if "aux_flag" not in coeff:
                            coeff["aux_flag"] = np.zeros(len(coeff["aux"]), np.int32)
                        if l_split[0] in coeff["aux"]:
                            xf = np.where(np.array(coeff["aux"]) == l_split[0])[0]
                            if len(xf) > 0:
                                coeff["aux_flag"][xf] = int(l_split[1])
                        if l_split[0] == "NP":
                            if "factor" not in coeff:
                                coeff["factor"] = (
                                    np.ones(len(coeff["freq"]), np.float32) * Fill_Value_Float
                                )
                                coeff["factor"][0] = float(l_split[1][:].split()[0])
                            else:
                                coeff["factor"] = np.hstack(
                                    (coeff["factor"], float(l_split[1][:].split()[0]))
                                )
                        if l_split[0] == "NS":
                            cll = l_split[1].split()
                            if "output_offset" not in coeff:
                                if prefix == "tpb":
                                    nn = (
                                        len(freq_coeff) * len(coeff["ele"])
                                        + np.sum(coeff["aux_flag"])
                                        + 1
                                    )
                                    coeff["freq_bl"] = coeff["freq"]
                                else:
                                    nn = len(freq_coeff) + np.sum(coeff["aux_flag"]) + 1
                                coeff["input_offset"] = (
                                    np.ones((1, nn), np.float32) * Fill_Value_Float
                                )
                                coeff["input_scale"] = (
                                    np.ones((1, nn), np.float32) * Fill_Value_Float
                                )
                                coeff["output_offset"] = (
                                    np.ones((1, len(coeff["height_grid"])), np.float32)
                                    * Fill_Value_Float
                                )
                                coeff["output_scale"] = (
                                    np.ones((1, len(coeff["height_grid"])), np.float32)
                                    * Fill_Value_Float
                                )
                                coeff["input_offset"][0, :] = np.array(
                                    [float(idx) for idx in cll], np.float32
                                )
                                cll = lines[ind_l + 1].split(":")[1].split()
                                coeff["input_scale"][0, :] = np.array(
                                    [float(idx) for idx in cll], np.float32
                                )
                                cll = lines[ind_l + 2].split(":")[1].split()
                                coeff["output_offset"][0, :] = np.array(
                                    [float(idx) for idx in cll], np.float32
                                )
                                cll = lines[ind_l + 3].split(":")[1].split()
                                coeff["output_scale"][0, :] = np.array(
                                    [float(idx) for idx in cll], np.float32
                                )
                            else:
                                coeff["input_offset"] = np.vstack(
                                    (
                                        coeff["input_offset"],
                                        np.array([float(idx) for idx in cll], np.float32),
                                    )
                                )
                                cll = lines[ind_l + 1].split(":")[1].split()
                                coeff["input_scale"] = np.vstack(
                                    (
                                        coeff["input_scale"],
                                        np.array([float(idx) for idx in cll], np.float32),
                                    )
                                )
                                cll = lines[ind_l + 2].split(":")[1].split()
                                coeff["output_offset"] = np.vstack(
                                    (
                                        coeff["output_offset"],
                                        np.array([float(idx) for idx in cll], np.float32),
                                    )
                                )
                                cll = lines[ind_l + 3].split(":")[1].split()
                                coeff["output_scale"] = np.vstack(
                                    (
                                        coeff["output_scale"],
                                        np.array([float(idx) for idx in cll], np.float32),
                                    )
                                )

                        if l_split[0] == "W1":
                            cll = l_split[1].split()
                            if "weights1" not in coeff:
                                if coeff["aux_flag"][7] == 1:
                                    nn += 1
                                if coeff["aux_flag"][8] == 1:
                                    nn += 1
                                coeff["weights1"] = np.zeros((1, nn, coeff["n_hidden"][0]))
                                coeff["weights1"][0, 0, :] = np.array(
                                    [float(idx) for idx in cll[0 : coeff["n_hidden"][0]]],
                                    np.float32,
                                )
                                for jj in range(nn - 1):
                                    cll = lines[ind_l + jj + 1].split(":")[1].split()
                                    coeff["weights1"][0, jj + 1, :] = np.array(
                                        [float(idx) for idx in cll[0 : coeff["n_hidden"][0]]],
                                        np.float32,
                                    )
                            elif coeff["weights1"].ndim < 3:
                                weights_1 = np.zeros((nn, coeff["n_hidden"][0]))
                                weights_1[0, :] = np.array(
                                    [float(idx) for idx in cll[0 : coeff["n_hidden"][0]]],
                                    np.float32,
                                )
                                for jj in range(nn - 1):
                                    cll = lines[ind_l + jj + 1].split(":")[1].split()
                                    weights_1[jj + 1, :] = np.array(
                                        [float(idx) for idx in cll[0 : coeff["n_hidden"][0]]],
                                        np.float32,
                                    )
                                coeff["weights1"] = np.stack((coeff["weights1"], weights_1))
                            else:
                                weights_1 = np.zeros((nn, coeff["n_hidden"][0]))
                                weights_1[0, :] = np.array(
                                    [float(idx) for idx in cll[0 : coeff["n_hidden"][0]]],
                                    np.float32,
                                )
                                for jj in range(nn - 1):
                                    cll = lines[ind_l + jj + 1].split(":")[1].split()
                                    weights_1[jj + 1, :] = np.array(
                                        [float(idx) for idx in cll[0 : coeff["n_hidden"][0]]],
                                        np.float32,
                                    )
                                coeff["weights1"] = np.vstack(
                                    (coeff["weights1"], weights_1[np.newaxis, ...])
                                )

                        if l_split[0] == "W2":
                            cll = l_split[1].split()
                            if "weights2" not in coeff:
                                coeff["weights2"] = (
                                    np.ones(
                                        (1, coeff["n_hidden"][0] + 1, len(coeff["height_grid"])),
                                        np.float32,
                                    )
                                    * Fill_Value_Float
                                )
                                for ialt, _ in enumerate(coeff["height_grid"]):
                                    cll = lines[ind_l + ialt].split()[1 : coeff["n_hidden"][0] + 2]
                                    coeff["weights2"][0, :, ialt] = np.array(
                                        [float(idx) for idx in cll[0 : coeff["n_hidden"][0] + 1]],
                                        np.float32,
                                    )
                            else:
                                coeff["weights2"] = np.vstack(
                                    (
                                        coeff["weights2"],
                                        np.reshape(
                                            np.array(
                                                [
                                                    float(idx)
                                                    for idx in cll[0 : coeff["n_hidden"][0] + 1]
                                                ],
                                                np.float32,
                                            ),
                                            (
                                                1,
                                                coeff["n_hidden"][0] + 1,
                                                len(coeff["height_grid"]),
                                            ),
                                        ),
                                    )
                                )

                f.close()

        if coeff["ret_type"] < 2:
            
            def f_offset(x):
                return np.array([coeff["offset"][(np.abs(coeff["ele"] - v)).argmin()] for v in x])

            def f_lin(x):
                return np.array(
                    [coeff["coeff_lin"][(np.abs(coeff["ele"] - v)).argmin(), :] for v in x]
                )

            def f_quad(x):
                return np.array(
                    [coeff["coeff_quad"][(np.abs(coeff["ele"] - v)).argmin(), :] for v in x]
                )
            
        else:
        
            def input_scale(x):
                return np.array(
                    [coeff["input_scale"][(np.abs(coeff["ele"] - v)).argmin(), :] for v in x]
                )

            def input_offset(x):
                return np.array(
                    [coeff["input_offset"][(np.abs(coeff["ele"] - v)).argmin(), :] for v in x]
                )

            def output_scale(x):
                return np.array(
                    [coeff["output_scale"][(np.abs(coeff["ele"] - v)).argmin(), :] for v in x]
                )

            def output_offset(x):
                return np.array(
                    [coeff["output_offset"][(np.abs(coeff["ele"] - v)).argmin(), :] for v in x]
                )

            def Weights1(x):
                return np.array(
                    [coeff["weights1"][(np.abs(coeff["ele"] - v)).argmin(), :, :] for v in x]
                )

            def Weights2(x):
                return np.array(
                    [coeff["weights2"][(np.abs(coeff["ele"] - v)).argmin(), :, :] for v in x]
                )

            def factor(x):
                return np.array([coeff["factor"][(np.abs(coeff["ele"] - v)).argmin()] for v in x])


    if str(c_list[0][-3:]).lower() == "ret":
        retrieval_type = ["linear regression", "quadratic regression", "neural network"]
        coeff["retrieval_type"] = retrieval_type[coeff["ret_type"]]
        coeff["retrieval_elevation_angles"] = str(coeff["ele"])
        coeff["retrieval_frequencies"] = str(coeff["freq"][:])
    if np.sum(coeff["aux_flag"]) == 0:
        coeff["retrieval_auxiliary_input"] = "no_surface"
    else:
        coeff["retrieval_auxiliary_input"] = "surface"
    coeff["retrieval_description"] = "Neural network"


    return (
        (coeff, f_offset, f_lin, f_quad)
        if (coeff["ret_type"] < 2)
        else (
            coeff,
            input_scale,
            input_offset,
            output_scale,
            output_offset,
            Weights1,
            Weights2,
            factor,
        )
    )
