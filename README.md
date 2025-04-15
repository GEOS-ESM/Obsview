# OBS2HTML
## Creating HTML with Observing System Summary Plots from netcdf5 ods files

##  SYNOPSIS
     obs2html -expid expid [...] odsfile
     obs2html -expid expid -begdate yyyymmddhh -enddate yyyymmddhh [...]

##  DESCRIPTION

     Creates an HTML document with observing system summary.
     Must be run in a writeable directory.
     Requires MATLAB.

     Use the first form to process one or more ods files, create
     data coverage plots, and compute/store various statistics.

     Use the second form to compute/plot time-averaged statistics
     over the period specified, to create time series plots, and
     to create html.

##  ARGUMENTS

      odsfile             name of ods file containing observations
     -expid   expid       identification string
     -begdate yyyymmddhh  start date for time series and statistics
     -enddate yyyymmddhh  end date for time series and statistics
    [-nocoverage]         do not plot data coverage
    [-noseries]           do not plot time series
    [-nostats]            do not plot statistics
    [-dirname dirname]    directory for storing results
                         (Default: obsmon)
    [-grpfile grpfile]    .m file containing data group definitions
                         (Default: obsgroups)
    [-grpnum  grpnum ]    Process only group number grpnum
                         (Default: process all groups)
    [-plotfmt plotfmt]    graphics format plots
                         (Default: png)
    [-plotdpi plotdpi]    graphics resolution, in dots per inch
                         (Default: 120)


## EXAMPLE
     obs2html -expid 'myexp' "*.ods"     (note the double quotes)
     obs2html -expid 'myexp' -begdate 2004080100 -enddate 2004080718
     netscape obsmon/index.html

## NOTE
     If the input filename contains wild cards (*)
     then it must be surrounded with double quotes.

