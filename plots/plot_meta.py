"""Metadata for plotting module."""
from typing import NamedTuple, Optional, Sequence, Tuple, Union


class PlotMeta(NamedTuple):
    name: Optional[str] = None
    cbar: Optional[Union[str, Sequence[str]]] = None
    clabel: Optional[Union[str, Sequence[Tuple[str, str]]]] = None
    ylabel: Optional[str] = None
    plot_range: Optional[Tuple[float, float]] = None
    plot_type: Optional[str] = None
    source: Optional[str] = None
    cbar_ext: Optional[str] = None
    nlev: Optional[int] = None
    ele: Optional[Tuple[float, float]] = None


_K = "K"
_KGM2 = "kg m$^{-2}$"
_GM2 = "g m$^{-2}$"
_KGM3 = "kg m$^{-3}$"
_GM3 = "g m$^{-3}$"
_PP = "%"
_DEG = "DEG"
_HPA = "hPa"
_MMH = "mm h$^{-1}$"
_MS = "m s$^{-1}$"

_COLORS = {
    "green": "#3cb371",
    "darkgreen": "#253A24",
    "lightgreen": "#70EB5D",
    "yellowgreen": "#C7FA3A",
    "yellow": "#FFE744",
    "orange": "#ffa500",
    "pink": "#B43757",
    "red": "#F57150",
    "shockred": "#E64A23",
    "seaweed": "#646F5E",
    "seaweed_roll": "#748269",
    "white": "#ffffff",
    "lightblue": "#6CFFEC",
    "blue": "#209FF3",
    "skyblue": "#CDF5F6",
    "darksky": "#76A9AB",
    "darkpurple": "#464AB9",
    "lightpurple": "#6A5ACD",
    "purple": "#BF9AFF",
    "darkgray": "#2f4f4f",
    "lightgray": "#ECECEC",
    "gray": "#d3d3d3",
    "lightbrown": "#CEBC89",
    "lightsteel": "#a0b0bb",
    "steelblue": "#4682b4",
    "mask": "#C8C8C8",
}

# Labels (and corresponding data) starting with an underscore are NOT shown:

_CLABEL = {
    "quality_flag_0": (
        ("receiver1_sanity_failed", _COLORS["seaweed_roll"]),
        ("receiver2_sanity_failed", _COLORS["lightbrown"]),
        ("rain_detected", _COLORS["blue"]),
        ("sun_in_beam", _COLORS["yellow"]),
    ),
    "quality_flag_1": (
        ("missing_tb", _COLORS["darkgray"]),
        ("spectral_consistency_failed", _COLORS["shockred"]),
    ),
    "quality_flag_2": (
        ("missing_tb", _COLORS["darkgray"]),
        ("tb_below_threshold", _COLORS["darkpurple"]),
        ("tb_above_threshold", _COLORS["orange"]),
    ),
    "quality_flag_3": (
        ("missing_tb", _COLORS["mask"]),
        ("tb_offset_above_threshold", _COLORS["lightbrown"]),
    ),
    "met_quality_flag": (
        ("low_quality_air_temperature", _COLORS["shockred"]),
        ("low_quality_relative_humidity", _COLORS["darkpurple"]),
        ("low_quality_air_pressure", _COLORS["lightbrown"]),
        ("low_quality_rain_rate", _COLORS["blue"]),
        ("low_quality_wind_direction", _COLORS["orange"]),
        ("low_quality_wind_speed", _COLORS["seaweed_roll"]),
    ),
}

_CBAR = {"bit": (_COLORS["white"], _COLORS["steelblue"])}

ATTRIBUTES = {
    "lwp": PlotMeta(
        name="Retrieved column-integrated liquid water path",
        cbar="Blues",
        ylabel=_KGM2,
        plot_range=(-0.05, 1.0),
        plot_type="bar",
        source="int",
        ele=(89.0, 91.0),
    ),
    "iwv": PlotMeta(
        name="Retrieved column-integrated water vapour",
        cbar="Blues",
        ylabel=_KGM2,
        plot_range=(0, 50),
        plot_type="bar",
        source="int",
    ),
    "water_vapor_vmr": PlotMeta(
        name="Retrieved water vapour profile",
        cbar="Spectral_r",
        clabel=_KGM3,
        plot_range=(-0.00000000000001, 0.02),
        plot_type="mesh",
        cbar_ext="max",
        nlev=11,
    ),
    "temperature": PlotMeta(
        # name='Retrieved temperature profile',
        cbar="RdBu_r",
        clabel=_K,
        plot_range=(245.0, 305.0),
        plot_type="mesh",
        cbar_ext="both",
        nlev=21,
    ),
    "potential_temperature": PlotMeta(
        cbar="inferno",
        clabel=_K,
        plot_range=(260.0, 320.0),
        plot_type="mesh",
        cbar_ext="both",
        nlev=31,
    ),
    "equivalent_potential_temperature": PlotMeta(
        cbar="turbo",
        clabel=_K,
        plot_range=(270.0, 350.0),
        plot_type="mesh",
        cbar_ext="both",
        nlev=41,
    ),
    "relative_humidity": PlotMeta(
        name="Relative humidity",
        cbar="viridis",
        clabel=_PP,
        plot_range=(0.0, 100.1),
        plot_type="mesh",
        cbar_ext="neither",
        source="met",
        nlev=11,
    ),
    "air_temperature": PlotMeta(
        name="Air and dewpoint temperature",
        ylabel=_K,
        plot_range=(240, 310),
        plot_type="bar",
        source="met",
    ),
    "air_pressure": PlotMeta(
        name="Air pressure",
        ylabel=_HPA,
        plot_range=(960, 1040),
        plot_type="bar",
        source="met",
    ),
    "rain_rate": PlotMeta(
        name="Precipitation amount",
        ylabel=_MMH,
        plot_range=(0, 50),
        plot_type="bar",
        source="met",
    ),
    "wind_direction": PlotMeta(
        name="Wind direction",
        ylabel=_DEG,
        plot_range=(0, 360),
        plot_type="bar",
        source="met",
    ),
    "wind_speed": PlotMeta(
        name="Wind speed",
        ylabel=_MS,
        plot_range=(0, 30),
        plot_type="bar",
        source="met",
    ),
    "met_quality_flag": PlotMeta(
        name="MET quality flag", clabel=_CLABEL["met_quality_flag"], source="mqf"
    ),
    "tb": PlotMeta(
        name="Microwave brightness temperatures",
        ylabel=_K,
        plot_type="bar",
        source="tb",
    ),
    "tb_spectrum": PlotMeta(
        name="Microwave brightness temperature spectrum",
        ylabel=_K,
        plot_type="bar",
        source="tb",
    ),
    "ele": PlotMeta(
        name="Sensor elevation angle",
        ylabel=_DEG,
        plot_range=(-1, 93),
        plot_type="bar",
        source="sen",
    ),
    "azi": PlotMeta(
        name="Sensor azimuth angle",
        ylabel=_DEG,
        plot_range=(-5, 363),
        plot_type="bar",
        source="sen",
    ),
    "quality_flag": PlotMeta(name="Quality flag", source="qf"),
    "quality_flag_0": PlotMeta(
        name="Sanity quality flag", clabel=_CLABEL["quality_flag_0"], source="qf"
    ),
    "quality_flag_1": PlotMeta(
        name="Spectral consistency quality flag",
        clabel=_CLABEL["quality_flag_1"],
        source="qf",
    ),
    "quality_flag_2": PlotMeta(
        name="TB quality flag", clabel=_CLABEL["quality_flag_2"], source="qf"
    ),
    "quality_flag_3": PlotMeta(
        name="TB offset quality flag", clabel=_CLABEL["quality_flag_3"], source="qf"
    ),
    "irt": PlotMeta(
        name="Infrared brightness temperatures",
        ylabel=_K,
        plot_range=(170, 310),
        plot_type="bar",
        source="irt",
    ),
    "t_amb": PlotMeta(
        name="Ambient target temperature",
        ylabel=_K,
        plot_type="bar",
        source="hkd",
    ),
    "t_rec": PlotMeta(
        name="Receiver physical temperature",
        ylabel=_K,
        plot_type="bar",
        source="hkd",
    ),
    "t_sta": PlotMeta(
        name="Receiver temperature stability",
        ylabel=_K,
        plot_type="bar",
        source="hkd",
    ),
}
