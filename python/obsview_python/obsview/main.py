#!/usr/bin/env python3
#Main script
#Used to initiate command line interface, read ODS or IODA files, compute stats, make plots

"""Entry point: parse args, read the file, compute+plot stats, save/show."""
import matplotlib.pyplot as plt

from . import cli
from . import config
from .io.readers import ioda_from_tarball, is_ods, is_ioda
from .stats.ods_stats import ods_pressure_binned, ods_channel
from .stats.ioda_stats import ioda_pressure_binned, ioda_channel
from .stats.compare_stats import jediXgsi_channel, jediXgsi_pressure_binned

#Print meta data of the file: observation type, variable name, data type (kt), number of bins, if there is scaling factor 
#Check if file is an ODS or IODA file
#Choose whether to show plot or save plot to a file

#def print_meta():

def run(args)-> None:
    cfg = config.resolve(args.obtype, args.var, args.scale)

    print(f"obtype: {args.obtype}")
    print(f"variable: {cfg.varname}")
    print(f"kt: {cfg.kt}")
    print(f"bins: {args.nbins}")
    print(f"scaled by: {cfg.scaleby}")

#Check if file is an ODS or IODA file
    if is_ods(args.filename):
        if args.obtype == "radiance": #Only radiance observation type uses channel
            ods_channel(args.filename, cfg.varname, cfg.levlim, args.nbins, args.satid, cfg.kt, args.qc)
        else: #Use pressure levels (these will be binned)
            ods_pressure_binned(
                args.filename, cfg.varname, cfg.levlim, args.nbins, cfg.scaleby,
                args.satid, cfg.kt, args.qc,
            )
#If not ODS, it must be IODA file            
    elif is_ioda(args.filename):
        #Stores file contents from either tarball or IODA file into variable
        nc = ioda_from_tarball(args.tarname, args.filename)
        #Check if user wants to make comparison plots
        if args.xGSI:
            #Use channel function if observation type is radiance
            if args.obtype == "radiance":
                jediXgsi_channel(
                    nc, cfg.varname, cfg.levlim, args.nbins, cfg.scaleby,
                    args.satid, args.qc, args.bias, args.common,
                )
            #Everything else uses binned pressure levels
            else:
                jediXgsi_pressure_binned(
                    nc, cfg.varname, cfg.levlim, args.nbins, cfg.scaleby, args.satid, args.qc
                )
        #If user does not want to make comparison plots, then make regular obsview plots
        else:
            if args.obtype == "radiance" or args.obtype == "aero": #Radiance and aerosol observation types both use channel levels
                ioda_channel(nc, cfg.varname, cfg.levlim, args.nbins, cfg.scaleby, args.satid, args.qc)
            else: #Everything else uses binned pressure levels
                ioda_pressure_binned(nc, cfg.varname, cfg.levlim, args.nbins, cfg.scaleby, args.satid, args.qc)
    else:
        print(f"The provided file is neither an ODS or IODA file, please try again...\n")

    #Check if user wants to save the figure to a file
    if args.fig == "none":
        plt.show()
    else:
        plt.savefig(args.fig, dpi=300, orientation="landscape", format="png")

#Passes command line arguments into args object
#Executes run() using said arguments
def main()-> None:
    args = cli.parse_args()
    run(args)


if __name__ == "__main__":
    main()
