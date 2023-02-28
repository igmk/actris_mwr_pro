#!/bin/bash
#!/usr/bin/env python3
"""A wrapper script for calling data processing functions."""

import argparse
import sys
import os
import utils
import process_mwrpy

"""All modules MUST have an add_arguments function which adds the subcommand to the subparser."""
modules = {
    "process": process_mwrpy,
}


def main(args):
    os.chdir("/home/hatpro/mwr_pro_actris/")
    args = _parse_args(args)
    cmd = args.cmd
    modules[cmd].main(args)


def _parse_args(args):
    parser = argparse.ArgumentParser(
        description="MWRpy processing main wrapper.", epilog="Have fun!"
    )
    subparsers = parser.add_subparsers(
        title="Command", help="Command to execute.", required=True, dest="cmd"
    )
    for module in modules.values():
        subparsers = module.add_arguments(subparsers)
    group = parser.add_argument_group(title="General options")
    group.add_argument(
        "-s", "--site", required=True, help="Site to process data from, e.g. juelich", type=str
    )
    group.add_argument(
        "-p",
        "--products",
        help="Products to be processed, e.g., 1C01, 2I01, 2P02.\
                        Default is all regular products.",
        type=lambda s: s.split(","),
        default=[
            "1C01",
            "2I01",
            "2I02",
            "2P01",
            "2P02",
            "2P03",
            "2P04",
            "2P07",
            "2P08",
            "2S02",
        ],
    )
    group.add_argument(
        "--start",
        type=str,
        metavar="YYYY-MM-DD",
        help="Starting date. Default is current day - 1 (included).",
        default=utils.get_date_from_past(1),
    )
    group.add_argument(
        "--stop",
        type=str,
        metavar="YYYY-MM-DD",
        help="Stopping date. Default is current day + 1 (excluded).",
        default=utils.get_date_from_past(-1),
    )
    group.add_argument(
        "-d", "--date", type=str, metavar="YYYY-MM-DD", help="Single date to be processed."
    )
    return parser.parse_args(args)


if __name__ == "__main__":
    main(sys.argv[1:])
