"""Command-line interface for the viewer."""

import argparse


def build_parser():
    p = argparse.ArgumentParser(
        prog="orca_vib_viewer.py",
        description="View ORCA vibrational modes and frontier-orbital character.",
    )
    p.add_argument("freq_file", nargs="?", metavar="FREQ.out",
                   help="ORCA output with a VIBRATIONAL FREQUENCIES block; "
                        "opens in the Vibrational Modes tab")
    p.add_argument("--pop", metavar="POP.log",
                   help="ORCA output with Loewdin populations (or the overlap and "
                        "MO matrices); loads into the Orbital Analysis tab at startup")
    p.add_argument("--groups", metavar="GROUPS.json",
                   help="group definitions saved from the Orbital Analysis tab")
    p.add_argument("--tab", choices=["modes", "orbital"],
                   help="tab to show first (default: orbital if --pop or --groups "
                        "is given without FREQ.out, else modes)")
    return p


def parse_args(argv=None):
    args = build_parser().parse_args(argv)
    if args.tab is None:
        args.tab = "orbital" if (args.pop or args.groups) and not args.freq_file else "modes"
    return args
