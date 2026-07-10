"""Command-line argument parsing."""
import argparse


def build_parser():
    parser = argparse.ArgumentParser(description="Pressure-binned statistics in file.")
    parser.add_argument("filename", help="Path to the log file")

    parser.add_argument("--scale", type=str, default="null",
                         help="Scale by obs/hofx0/etc (default: False)")
    parser.add_argument("--obtype", type=str, default="mls55_aura",
                         help="Observation type (default: mls55_aura)")
    parser.add_argument("--var", type=str, default="auto",
                         help="Variable name (default: auto)")
    parser.add_argument("--tarname", type=str, default="none",
                         help="Tar file name (default: none)")
    parser.add_argument("--satid", type=int, default=-999,
                         help="Observation type (default: -999 (all))")
    parser.add_argument("--nbins", type=int, default=40,
                         help="Number of bins (default: 40 (all))")
    parser.add_argument("--qc", type=int, default=0,
                         help="Quality mark (default: 0 (used obs))")
    parser.add_argument("--xGSI", action="store_true",
                         help="Comp with GSI(default: True)")
    parser.add_argument("--common", action="store_true",
                         help="When xGSI, use only common used in comp: True)")
    parser.add_argument("--bias", action="store_true",
                         help="When xGSI, comp bias as opposed to sigo: True)")
    parser.add_argument("--fig", type=str, default="none",
                         help="give a filename to save plot (default: none)")

    return parser


def parse_args(argv=None):
    return build_parser().parse_args(argv)
