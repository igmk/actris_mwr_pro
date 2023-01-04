"""Module for processing."""
import os
import sys
import datetime
from level1.write_lev1_nc import lev1_to_nc
from level2.write_lev2_nc import lev2_to_nc
from plots.generate_plots import generate_figure
from read_specs import get_site_specs

site = sys.argv[1]
dvec = sys.argv[2]
do_plot = int(sys.argv[3])

date = datetime.datetime.strptime(dvec, "%Y%m%d")
prv = date - datetime.timedelta(days=1)
nex = date + datetime.timedelta(days=1)

global_attributes, params = get_site_specs(site, "1C01")
ID = global_attributes["wigos_station_id"]
data_out_l1 = params["data_out"] + "level1/" + date.strftime("%Y/%m/%d/")
data_out_l2 = params["data_out"] + "level2/" + date.strftime("%Y/%m/%d/")
product = [
    ("1C01", ""),
    ("2I01", "lwp"),
    ("2I02", "iwv"),
    ("2P01", "temperature"),
    ("2P02", "temperature"),
    ("2P03", "water_vapor_vmr"),
    ("2P04", "relative_humidity"),
    ("2P07", "potential_temperature"),
    ("2P08", "equivalent_potential_temperature"),
    ("2S02", "tb_spectrum"),
]


def _link_quicklook(link_direc: str, figure_name: str) -> None:
    if not os.path.isdir(link_direc):
        os.makedirs(link_direc)
    link_rr = link_direc + figure_name.rsplit("/")[-1]

    if not os.path.islink(link_rr):
        print('linking: ln -s "%s" "%s"' % (figure_name, link_rr))
        os.symlink(figure_name, link_rr)


for prod, var in product:
    if prod == "1C01":
        link_dir = (
            "/home/hatpro/public_html/quicklooks/"
            + params["data_out"][6:]
            + "level1/"
            + date.strftime("%Y/%m/%d/")
        )
        if not os.path.isdir(data_out_l1):
            os.makedirs(data_out_l1)
        lev1_data = data_out_l1 + "MWR_" + prod + "_" + ID + "_" + date.strftime("%Y%m%d") + ".nc"
        lev1_to_nc(
            site,
            prod,
            params["data_in"] + date.strftime("%Y/%m/%d/"),
            params["data_in"] + prv.strftime("%Y/%m/%d/"),
            params["data_in"] + nex.strftime("%Y/%m/%d/"),
            lev1_data,
        )
        if do_plot:
            fig_name = generate_figure(
                lev1_data,
                ["tb"],
                ele_range=[89.0, 91.0],
                save_path=data_out_l1,
                image_name="tb",
            )
            _link_quicklook(link_dir, fig_name)
            fig_name = generate_figure(
                lev1_data, ["ele", "azi"], save_path=data_out_l1, image_name="sen"
            )
            _link_quicklook(link_dir, fig_name)
            fig_name = generate_figure(
                lev1_data,
                ["quality_flag"],
                save_path=data_out_l1,
                image_name="quality_flag",
            )
            _link_quicklook(link_dir, fig_name)
            fig_name = generate_figure(
                lev1_data,
                ["met_quality_flag"],
                save_path=data_out_l1,
                image_name="met_quality_flag",
            )
            _link_quicklook(link_dir, fig_name)
            fig_name = generate_figure(
                lev1_data,
                ["t_amb", "t_rec", "t_sta"],
                save_path=data_out_l1,
                image_name="hkd",
            )
            _link_quicklook(link_dir, fig_name)
            fig_name = generate_figure(
                lev1_data,
                [
                    "air_temperature",
                    "relative_humidity",
                    "rain_rate",
                ],
                save_path=data_out_l1,
                image_name="met",
            )
            _link_quicklook(link_dir, fig_name)
            fig_name = generate_figure(
                lev1_data,
                [
                    "air_pressure",
                    "wind_direction",
                    "wind_speed",
                ],
                save_path=data_out_l1,
                image_name="met2",
            )
            _link_quicklook(link_dir, fig_name)
            if params["ir_beamwidth"] != -999.0:
                fig_name = generate_figure(
                    lev1_data,
                    ["irt"],
                    ele_range=[89.0, 91.0],
                    save_path=data_out_l1,
                    image_name="irt",
                )
                _link_quicklook(link_dir, fig_name)

    elif os.path.isdir("site_config/" + site + "/coefficients/"):
        link_dir = (
            "/home/hatpro/public_html/quicklooks/"
            + params["data_out"][6:]
            + "level2/"
            + date.strftime("%Y/%m/%d/")
        )
        if not os.path.isdir(data_out_l2):
            os.makedirs(data_out_l2)
        lev1_data = data_out_l1 + "MWR_1C01_" + ID + "_" + date.strftime("%Y%m%d") + ".nc"
        lev2_to_nc(date.strftime("%Y%m%d"), site, prod, lev1_data, data_out_l2)
        if prod in ("2I01", "2I02"):
            elevation = [89.0, 91.0]
        else:
            elevation = [0.0, 91.0]
        if (do_plot) & (
            os.path.isfile(
                data_out_l2 + "MWR_" + prod + "_" + ID + "_" + date.strftime("%Y%m%d") + ".nc"
            )
        ):
            fig_name = generate_figure(
                data_out_l2 + "MWR_" + prod + "_" + ID + "_" + date.strftime("%Y%m%d") + ".nc",
                [var],
                ele_range=elevation,
                save_path=data_out_l2,
                image_name=var,
            )
            _link_quicklook(link_dir, fig_name)
