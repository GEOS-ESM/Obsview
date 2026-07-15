#Functions to parse(break down) command line arguments

"""Command-line argument parsing."""
import argparse

#Create the parser and set up arguments
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Pressure-binned statistics in file.") #Create argument parser object
    parser.add_argument("filename", help="Path to the log file")
    #Scales stats for levels with low magnitude values
    parser.add_argument("--scale", type=str, default="null",                            #Only applies to radio occultation(RO) data from GPS satellites 
                         help="Scale by obs/hofx0/etc (default: False)")  
    #Observation type: User sets what type of observation should get plotted              
    parser.add_argument("--obtype", type=str, default="mls55_aura",
                         help="Observation type (default: mls55_aura)")
    #Similar to observation type but user inputs variable name(e.g. ozoneProfile instead of msl55_aura)
    parser.add_argument("--var", type=str, default="auto",
                         help="Variable name (default: auto)")
    #Used when user wants to read IODA files from tarball
    parser.add_argument("--tarname", type=str, default="none",
                         help="Tar file name (default: none)")
    #SatID == Observation ID for IODA files, table found in ostats.rc file
    parser.add_argument("--satid", type=int, default=-999,
                         help="Observation type (default: -999 (all))")
    #Set number of bins for pressure or channel levels
    parser.add_argument("--nbins", type=int, default=40,
                         help="Number of bins (default: 40 (all))")
    #See observations used, passive, unused
    parser.add_argument("--qc", type=int, default=0,
                         help="Quality mark (default: 0 (used obs))")
    #Make plots comparing JEDI(newer) stats to GSI(older)
    parser.add_argument("--xGSI", action="store_true",
                         help="Compare with GSI(default: True)")
    #Compare only the levels that contain both JEDI and GSI observation
    parser.add_argument("--common", action="store_true",                                #Only used when --xGSI is specified
                         help="When xGSI, use only common used in comp: True)")
    #Compare bias instead of sigo
    parser.add_argument("--bias", action="store_true",                                  #Only used when --xGSI is specified
                         help="When xGSI, compare bias as opposed to sigo: True)")
    #Specify whether the program displays plots or saves them to a file
    parser.add_argument("--fig", type=str, default="none",
                         help="give a filename to save plot (default: none)")

    return parser

#Return Namespace object containing parsed data
def parse_args(argv=None) -> argparse.Namespace:
    return build_parser().parse_args(argv)
