import netCDF4 as nc
import numpy as np
from numpy import ma

from utils import get_coeff_list

Fill_Value_Float = -999.0
Fill_Value_Int = -99


def get_mvr_coeff(site: str, prefix: str, freq: np.ndarray):
    """This function extracts retrieval coefficients for given files
    and performs elevation angle interpolation for integrated quantities

    Args:
        site: Name of site.
        prefix: Identifier for type of product.
        freq: Frequencies of observations.

    Examples:
        >>> from level2.get_ret_coeff import get_mvr_coeff
        >>> get_mvr_coeff('site_name', 'lwp', [22, 31.4])
    """

    c_list = get_coeff_list(site, prefix)
    cf = {}
    cf["rt"] = -99

    if (str(c_list[0][-3:]).lower() == "ret") & (len(c_list) == 1):
        with open(c_list[0], "r") as f:
            lines = f.readlines()
            lines = [line.rstrip("\n") for line in lines]
            for line in lines:
                ls = line.split("=")
                if ls[0] == "AG":
                    N = len([float(x) for x in ls[1].split()])
    else:
        N = len(c_list)

    if prefix in ("lwp", "iwv"):

        cf["ele"] = np.ones(N) * Fill_Value_Float
        cf["freq"] = np.ones([len(freq), N]) * Fill_Value_Float
        cf["coeff_lin"] = np.zeros([N, 23])
        cf["coeff_quad"] = np.zeros([N, 23])
        cf["offset"] = np.zeros(N)
        cf["err_ran"] = np.ones(N) * Fill_Value_Float
        cf["err_sys"] = np.ones(N) * Fill_Value_Float

        if str(c_list[0][-3:]).lower() == "ret":
            for i, file in enumerate(c_list):
                with open(file, "r") as f:
                    lines = f.readlines()
                    lines = [line.rstrip("\n") for line in lines]

                    for il, line in enumerate(lines):
                        ls = line.split("=")
                        if ls[0] == "RT":
                            cf["rt"] = int(ls[1].split("#")[0])
                        if ls[0] == "AG":
                            cf["ele"][i:N] = np.array([float(x) for x in ls[1].split()], np.float32)
                        if ls[0] == "FR":
                            freq_ret = np.array(
                                [float(idx) for idx in ls[1][:].split()], np.float32
                            )
                            _, freq_ind, freq_cf = np.intersect1d(
                                freq[:],
                                freq_ret[:],
                                assume_unique=False,
                                return_indices=True,
                            )
                            if len(freq_cf) < len(freq_ret):
                                raise RuntimeError(
                                    ["Instrument and retrieval frequencies do not match."]
                                )
                            cf["freq"][freq_ind, i:N] = np.resize(
                                freq_ret[freq_cf], (len(freq_cf), N)
                            )

                        # Regression
                        if (cf["rt"] != -99) & (cf["rt"] < 2):
                            cf["aux"] = ["TS", "HS", "PS", "IS"]
                            if "aux_flg" not in cf:
                                cf["aux_flg"] = np.zeros(len(cf["aux"]), np.int32)
                            for ii, aux_i in enumerate(cf["aux"]):
                                if ls[0] == aux_i:
                                    cf["aux_flg"][ii] = int(ls[1])
                            if ls[0] == "OS":
                                for jj in range(N):
                                    if ":" in lines[il + jj]:
                                        cll = lines[il + jj].split(":")
                                    else:
                                        cll = lines[il + jj].split("=")
                                    cf["offset"][i + jj] = float(cll[1].split("#")[0])
                            if ls[0] == "TL":
                                for jj in range(N):
                                    if ":" in lines[il + jj]:
                                        cll = lines[il + jj].split(":")
                                    else:
                                        cll = lines[il + jj].split("=")
                                    cf["coeff_lin"][i + jj, freq_ind] = np.array(
                                        [float(idx) for idx in cll[1].split()[0 : len(freq_cf)]],
                                        np.float32,
                                    )
                            if (ls[0] == "TQ") & (cf["rt"] == 1):
                                for jj in range(N):
                                    if ":" in lines[il + jj]:
                                        cll = lines[il + jj].split(":")
                                    else:
                                        cll = lines[il + jj].split("=")
                                    cf["coeff_quad"][i + jj, freq_ind] = np.array(
                                        [float(idx) for idx in cll[1].split()[0 : len(freq_cf)]],
                                        np.float32,
                                    )

                        # Neural Network
                        elif (cf["rt"] != -99) & (cf["rt"] == 2):
                            if ls[0] == "ND":
                                cf["nd"] = np.array(
                                    [int(idx) for idx in ls[1][:].split()[0:2]],
                                    np.int32,
                                )
                            cf["aux"] = [
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
                            if "aux_flg" not in cf:
                                cf["aux_flg"] = np.zeros(len(cf["aux"]), np.int32)
                            if ls[0] in cf["aux"]:
                                xf = np.where(np.array(cf["aux"]) == ls[0])[0]
                                cf["aux_flg"][xf] = int(ls[1])
                            if ls[0] == "NP":
                                if "np" not in cf:
                                    cf["np"] = float(ls[1][:].split()[0])
                                else:
                                    cf["np"] = np.hstack((cf["np"], float(ls[1][:].split()[0])))
                            if ls[0] == "NS":
                                cll = ls[1].split()
                                if "op_os" not in cf:
                                    cf["in_os"] = np.array([float(idx) for idx in cll], np.float32)
                                    cll = lines[il + 1].split(":")[1].split()
                                    cf["in_sc"] = np.array([float(idx) for idx in cll], np.float32)
                                    cf["op_os"] = float(lines[il + 2].split(":")[1].split()[0])
                                    cf["op_sc"] = float(lines[il + 3].split(":")[1].split()[0])
                                else:
                                    cf["in_os"] = np.vstack(
                                        (
                                            cf["in_os"],
                                            np.array([float(idx) for idx in cll], np.float32),
                                        )
                                    )
                                    cll = lines[il + 1].split(":")[1].split()
                                    cf["in_sc"] = np.vstack(
                                        (
                                            cf["in_sc"],
                                            np.array([float(idx) for idx in cll], np.float32),
                                        )
                                    )
                                    cf["op_os"] = np.vstack(
                                        (
                                            cf["op_os"],
                                            float(lines[il + 2].split(":")[1].split()[0]),
                                        )
                                    )
                                    cf["op_sc"] = np.vstack(
                                        (
                                            cf["op_sc"],
                                            float(lines[il + 3].split(":")[1].split()[0]),
                                        )
                                    )
                                    # cf['ns_ta'] = np.vstack((cf['ns_ta'], 
                                    # np.array([float(lines[il+2].split(':')[1].split()[0]), 
                                    # float(lines[il+3].split(':')[1].split()[0])])))
                            if ls[0] == "W1":
                                cll = ls[1].split()
                                if "w1" not in cf:
                                    nn = np.sum(cf["aux_flg"]) + 15
                                    if cf["aux_flg"][7] == 1:
                                        nn += 1
                                    if cf["aux_flg"][8] == 1:
                                        nn += 1
                                    cf["w1"] = np.zeros((nn, cf["nd"][0]))
                                    cf["w1"][0, :] = np.array(
                                        [float(idx) for idx in cll[0 : cf["nd"][0]]],
                                        np.float32,
                                    )
                                    for jj in range(nn - 1):
                                        cll = lines[il + jj + 1].split(":")[1].split()
                                        cf["w1"][jj + 1, :] = np.array(
                                            [float(idx) for idx in cll[0 : cf["nd"][0]]],
                                            np.float32,
                                        )
                                elif cf["w1"].ndim < 3:
                                    w1 = np.zeros((nn, cf["nd"][0]))
                                    w1[0, :] = np.array(
                                        [float(idx) for idx in cll[0 : cf["nd"][0]]],
                                        np.float32,
                                    )
                                    for jj in range(nn - 1):
                                        cll = lines[il + jj + 1].split(":")[1].split()
                                        w1[jj + 1, :] = np.array(
                                            [float(idx) for idx in cll[0 : cf["nd"][0]]],
                                            np.float32,
                                        )
                                    cf["w1"] = np.stack((cf["w1"], w1))
                                else:
                                    w1 = np.zeros((nn, cf["nd"][0]))
                                    w1[0, :] = np.array(
                                        [float(idx) for idx in cll[0 : cf["nd"][0]]],
                                        np.float32,
                                    )
                                    for jj in range(nn - 1):
                                        cll = lines[il + jj + 1].split(":")[1].split()
                                        w1[jj + 1, :] = np.array(
                                            [float(idx) for idx in cll[0 : cf["nd"][0]]],
                                            np.float32,
                                        )
                                    cf["w1"] = np.vstack((cf["w1"], w1[np.newaxis, ...]))
                            if ls[0] == "W2":
                                cll = ls[1].split()
                                if "w2" not in cf:
                                    cf["w2"] = np.array(
                                        [float(idx) for idx in cll[0 : cf["nd"][0]]],
                                        np.float32,
                                    )
                                else:
                                    cf["w2"] = np.vstack(
                                        (
                                            cf["w2"],
                                            np.array(
                                                [float(idx) for idx in cll[0 : cf["nd"][0]]],
                                                np.float32,
                                            ),
                                        )
                                    )

                    f.close()

        else:
            for i, file in enumerate(c_list):
                coeff = nc.Dataset(file)
                cf["ele"][i] = coeff["elevation_predictor"][i]
                _, freq_ind, freq_cf = np.intersect1d(
                    freq[:], coeff["freq"][:], assume_unique=False, return_indices=True
                )
                if len(freq_cf) < len(coeff["freq"][:]):
                    raise RuntimeError(["Instrument and retrieval frequencies do not match."])

                cf["freq"][freq_ind, i] = coeff["freq"][freq_cf]
                cf["coeff_lin"][i, freq_ind] = coeff["coefficient_mvr"][freq_cf]
                if coeff.regression_type == "quadratic":
                    cf["coeff_quad"][i, freq_ind] = coeff["coefficient_mvr"][freq_cf + len(freq_cf)]
                cf["offset"][i] = coeff["offset_mvr"][0]
                cf["err_ran"][i] = coeff["predictand_err"][0]
                cf["err_sys"][i] = coeff["predictand_err_sys"][0]

        if cf["rt"] < 2:

            def f_offset(x):
                return np.array([cf["offset"][(np.abs(cf["ele"] - v)).argmin()] for v in x])

            def f_lin(x):
                return np.array([cf["coeff_lin"][(np.abs(cf["ele"] - v)).argmin(), :] for v in x])

            def f_quad(x):
                return np.array([cf["coeff_quad"][(np.abs(cf["ele"] - v)).argmin(), :] for v in x])

            def e_ran(x):
                return np.array([cf["err_ran"][(np.abs(cf["ele"] - v)).argmin()] for v in x])

            def e_sys(x):
                return np.array([cf["err_sys"][(np.abs(cf["ele"] - v)).argmin()] for v in x])

        else:

            def ns_ta(x):
                return np.array([cf["ns_ta"][(np.abs(cf["ele"] - v)).argmin(), :] for v in x])

            def ns_sc(x):
                return np.array([cf["ns_sc"][(np.abs(cf["ele"] - v)).argmin(), :] for v in x])

            def ns_os(x):
                return np.array([cf["ns_os"][(np.abs(cf["ele"] - v)).argmin(), :] for v in x])

            def W1(x):
                return np.array([cf["w1"][(np.abs(cf["ele"] - v)).argmin(), :, :] for v in x])

            def W2(x):
                return np.array([cf["w2"][(np.abs(cf["ele"] - v)).argmin(), :] for v in x])

            def pn(x):
                return np.array([cf["np"][(np.abs(cf["ele"] - v)).argmin()] for v in x])

    elif prefix in ("tze", "hze"):

        if cf["rt"] < 2:
            coeff = nc.Dataset(c_list[0])
            n_height_grid = coeff.dimensions["n_height_grid"].size

            cf["ele"] = np.ones(N) * Fill_Value_Float
            cf["freq"] = np.ones([len(freq), N]) * Fill_Value_Float
            cf["coeff_lin"] = np.zeros([N, n_height_grid, 23])
            cf["coeff_quad"] = np.zeros([N, n_height_grid, 23])
            cf["offset"] = np.zeros([n_height_grid, N])
            cf["err_ran"] = ma.masked_all((n_height_grid, N))
            cf["err_sys"] = ma.masked_all((n_height_grid, N))
            cf["n_height_grid"] = n_height_grid
            cf["height_grid"] = coeff["height_grid"]

            for i, file in enumerate(c_list):
                coeff = nc.Dataset(file)
                cf["ele"][i] = coeff["elevation_predictor"][i]
                _, freq_ind, freq_cf = np.intersect1d(
                    freq[:], coeff["freq"][:], assume_unique=False, return_indices=True
                )
                if len(freq_cf) < len(coeff["freq"][:]):
                    raise RuntimeError(["Instrument and retrieval frequencies do not match."])

                cf["freq"][freq_ind, i] = coeff["freq"][freq_cf]
                cf["coeff_lin"][i, :, freq_ind] = coeff["coefficient_mvr"][freq_cf, :]
                if coeff.regression_type == "quadratic":
                    cf["coeff_quad"][i, :, freq_ind] = coeff["coefficient_mvr"][
                        freq_cf + len(freq_cf), :
                    ]
                cf["offset"][:, i] = coeff["offset_mvr"][:]
                cf["err_ran"][:, i] = coeff["predictand_err"][:]
                cf["err_sys"][:, i] = coeff["predictand_err_sys"][:]

            def f_offset(x):
                return np.array([cf["offset"][:, (np.abs(cf["ele"] - v)).argmin()] for v in x])

            def f_lin(x):
                return np.array(
                    [cf["coeff_lin"][(np.abs(cf["ele"] - v)).argmin(), :, :] for v in x]
                )

            def f_quad(x):
                return np.array(
                    [cf["coeff_quad"][(np.abs(cf["ele"] - v)).argmin(), :, :] for v in x]
                )

            def e_ran(x):
                return np.array([cf["err_ran"][:, (np.abs(cf["ele"] - v)).argmin()] for v in x])

            def e_sys(x):
                return np.array([cf["err_sys"][:, (np.abs(cf["ele"] - v)).argmin()] for v in x])

    elif prefix == "tel":

        coeff = nc.Dataset(c_list[0])
        _, freq_ind, freq_cf = np.intersect1d(
            freq[:], coeff["freq"][:], assume_unique=False, return_indices=True
        )
        if len(freq_cf) < len(coeff["freq"][:]):
            raise RuntimeError(["Instrument and retrieval frequencies do not match."])

        cf["ele"] = coeff["elevation_predictor"][:]
        cf["height_grid"] = coeff["height_grid"]
        cf["freq"] = coeff["freq"]
        cf["freq_bl"] = coeff["freq_bl"]
        cf["n_height_grid"] = coeff.dimensions["n_height_grid"].size
        f_offset = coeff["offset_mvr"][:]
        f_lin, f_quad = coeff["coefficient_mvr"][:, :], []
        e_ran = coeff["predictand_err"][:]
        e_sys = coeff["predictand_err_sys"][:]

    elif prefix == "tbx":
        cf["ele"] = np.ones((N, 1)) * Fill_Value_Float
        coeff = nc.Dataset(c_list[0])

        cf["ele"] = coeff["elevation_predictor"][:]
        cf["freq"] = freq[:]
        e_ran = coeff["predictand_err"][:]
        e_sys = coeff["predictand_err_sys"][:]
        f_offset, f_lin, f_quad, = (
            [],
            [],
            [],
        )

    else:
        raise RuntimeError(
            ["Prefix " + prefix + " not recognized for retrieval coefficient file(s)."]
        )

    if str(c_list[0][-3:]).lower() == "ret":
        retrieval_type = ["linear regression", "quadratic regression", "neural network"]
        cf["retrieval_type"] = retrieval_type[cf["rt"]]
        cf["retrieval_elevation_angles"] = str(cf["ele"])
        cf["retrieval_frequencies"] = str(cf["freq"][:])
        if np.sum(cf["aux_flg"]) == 0:
            cf["retrieval_auxiliary_input"] = "no_surface"
        else:
            cf["retrieval_auxiliary_input"] = "surface"
        cf["retrieval_description"] = "Neural network"

    else:
        cf["retrieval_type"] = coeff.regression_type
        cf["retrieval_elevation_angles"] = str(cf["ele"])
        cf["retrieval_frequencies"] = str(coeff["freq"][:])
        cf["retrieval_auxiliary_input"] = coeff.surface_mode
        cf["retrieval_description"] = coeff.retrieval_version

    return (
        (cf, f_offset, f_lin, f_quad, e_ran, e_sys)
        if (cf["rt"] < 2)
        else (cf, ns_ta, ns_sc, ns_os, W1, W2, pn)
    )
