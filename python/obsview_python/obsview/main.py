#!/usr/bin/env python3
"""Entry point: parse args, read the file, compute+plot stats, save/show."""
import matplotlib.pyplot as plt

from . import cli
from . import config
from .io.readers import ioda_from_tarball, is_ods
from .stats.ods_stats import ods_pressure_binned, ods_channel
from .stats.ioda_stats import ioda_pressure_binned, ioda_channel
from .stats.compare_stats import jediXgsi_channel, jediXgsi_pressure_binned


def run(args):
    cfg = config.resolve(args.obtype, args.var, args.scale)

    print(f"obtype: {args.obtype}")
    print(f"variable: {cfg.varname}")
    print(f"kt: {cfg.kt}")
    print(f"bins: {args.nbins}")
    print(f"scaled by: {cfg.scaleby}")

    if is_ods(args.filename):
        if args.obtype == "radiance":
            ods_channel(args.filename, cfg.varname, cfg.levlim, args.nbins, args.satid, cfg.kt, args.qc)
        else:
            ods_pressure_binned(
                args.filename, cfg.varname, cfg.levlim, args.nbins, cfg.scaleby,
                args.satid, cfg.kt, args.qc,
            )
    else:
        nc = ioda_from_tarball(args.tarname, args.filename)
        if args.xGSI:
            if args.obtype == "radiance":
                jediXgsi_channel(
                    nc, cfg.varname, cfg.levlim, args.nbins, cfg.scaleby,
                    args.satid, args.qc, args.bias, args.common,
                )
            else:
                jediXgsi_pressure_binned(
                    nc, cfg.varname, cfg.levlim, args.nbins, cfg.scaleby, args.satid, args.qc
                )
        else:
            if args.obtype == "radiance" or args.obtype == "aero":
                ioda_channel(nc, cfg.varname, cfg.levlim, args.nbins, cfg.scaleby, args.satid, args.qc)
            else:
                ioda_pressure_binned(nc, cfg.varname, cfg.levlim, args.nbins, cfg.scaleby, args.satid, args.qc)

    if args.fig == "none":
        plt.show()
    else:
        plt.savefig(args.fig, dpi=300, orientation="landscape", format="png")


def main():
    args = cli.parse_args()
    run(args)


if __name__ == "__main__":
    main()
