#Module for creating the 4 panel plot (called the 'Statistics plot' in Obsview)
#This plot will be a monthly averaged plot in the future when this program can read multiple files

"""Individual subplot panels shared by every stats-plotting entry point.

Each function draws into the *current* matplotlib axes (via plt.subplot,
called by the caller beforehand) - they don't create figures themselves.
"""
import numpy as np
import matplotlib.pyplot as plt

#Create an x-axis label for specific variable names
def _xlabel_for(varname):
    if varname == "ozoneProfile":
        return "Ozone Value (mol mol$^{-1}$)"
    if varname == "virtualTemperature":
        return "Virtual Temperature (K)"
    if varname == "bendingAngle":
        return "Bending Angle"
    return None

#Function to create the residual stats panel
def show_resstats_panel(
    varname, mean_ombg, mean_oman, rms_ombg, rms_oman, bin_centers, bin_heights, labels, radiance
):
    bar_width = bin_heights * 0.4
    offsets = [-1.5, -0.5, 0.5, 1.5]

    labeled = {"mean_ombg": False, "mean_oman": False, "rms_ombg": False, "rms_oman": False}
    for i in range(len(bin_centers)):
        y = bin_centers[i]

        plt.barh(y + offsets[0] * bar_width[i], mean_ombg[i], height=bar_width[i],
                 color="cyan", label=labels[0] if not labeled["mean_ombg"] else "")
        labeled["mean_ombg"] = True

        plt.barh(y + offsets[2] * bar_width[i], mean_oman[i], height=bar_width[i],
                 color="orange", label=labels[1] if not labeled["mean_oman"] else "")
        labeled["mean_oman"] = True

        plt.barh(y + offsets[1] * bar_width[i], rms_ombg[i], height=bar_width[i],
                 color="blue", label=labels[2] if not labeled["rms_ombg"] else "")
        labeled["rms_ombg"] = True

        plt.barh(y + offsets[3] * bar_width[i], rms_oman[i], height=bar_width[i],
                 color="red", label=labels[3] if not labeled["rms_oman"] else "")
        labeled["rms_oman"] = True

    if not radiance:
        plt.yscale("log")
        plt.gca().invert_yaxis()
    xlabel = _xlabel_for(varname)
    if xlabel:
        plt.xlabel(xlabel)
    if radiance:
        plt.ylabel("Channel")
        plt.title("Mean & RMS of Obs Residuals vs Channel")
    else:
        plt.ylabel("Pressure (hPa)")
        plt.title("Mean & RMS of Obs Residuals vs Pressure")
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.legend()

#Function to create the jo panel
def show_jo_panel(varname, mean_job, mean_joa, bin_centers, bin_heights, labels, radiance):
    bar_width = bin_heights * 0.4
    offsets = [-1.5, -0.5, 0.5, 1.5]

    labeled = {"mean_job": False, "mean_joa": False}
    for i in range(len(bin_centers)):
        y = bin_centers[i]

        plt.barh(y + offsets[0] * bar_width[i], mean_job[i], height=bar_width[i],
                 color="blue", label=labels[0] if not labeled["mean_job"] else "")
        labeled["mean_job"] = True

        plt.barh(y + offsets[1] * bar_width[i], mean_joa[i], height=bar_width[i],
                 color="red", label=labels[1] if not labeled["mean_joa"] else "")
        labeled["mean_joa"] = True

    if not radiance:
        plt.yscale("log")
        plt.gca().invert_yaxis()
    xlabel = _xlabel_for(varname)
    if xlabel:
        plt.xlabel(xlabel)
    if radiance:
        plt.ylabel("Channel")
        plt.title("Jo/p vs Channel")
    else:
        plt.ylabel("Pressure (hPa)")
        plt.title("Jo/p vs Pressure")
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.legend()

#Function to make the sigo panel
def show_sigo_panel(varname, mean_sigo, mean_esigo, mean_esigb, bin_centers, bin_heights, labels, radiance):
    bar_width = bin_heights * 0.4
    offsets = [-1.5, -0.5, 0.5, 1.5]

    labeled = {"mean_sigo": False, "mean_esigo": False, "mean_esigb": False}
    for i in range(len(bin_centers)):
        y = bin_centers[i]

        plt.barh(y + offsets[0] * bar_width[i], mean_sigo[i], height=bar_width[i],
                 color="cyan", label=labels[0] if not labeled["mean_sigo"] else "")
        labeled["mean_sigo"] = True

        plt.barh(y + offsets[1] * bar_width[i], mean_esigo[i], height=bar_width[i],
                 color="orange", label=labels[1] if not labeled["mean_esigo"] else "")
        labeled["mean_esigo"] = True

        if len(labels) > 2:
            plt.barh(y + offsets[2] * bar_width[i], mean_esigb[i], height=bar_width[i],
                     color="black", label=labels[2] if not labeled["mean_esigb"] else "")
            labeled["mean_esigb"] = True

    if not radiance:
        plt.yscale("log")
        plt.gca().invert_yaxis()
    xlabel = _xlabel_for(varname)
    if xlabel:
        plt.xlabel(xlabel)
    if radiance:
        plt.ylabel("Channel")
        plt.title("Prescribed & Estimated Errors vs Channel")
    else:
        plt.ylabel("Pressure (hPa)")
        plt.title("Prescribed & Estimated Errors vs Pressure")
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.ticklabel_format(axis="x", style="sci", scilimits=(0, 10))
    plt.legend()

#Function to create the number of observations panel
def show_nobs_panel(varname, sum_nobs, sum_nonobs, bin_centers, bin_heights, usrqc, show2, radiance):
    bar_width = bin_heights * 0.4
    offsets = [-0.5, -0.25, 0.5, 1.5]

    labeled = {"sum_nobs": False, "sum_nonobs": False}
    for i in range(len(bin_centers)):
        y = bin_centers[i]

        plt.barh(y + offsets[0] * bar_width[i], sum_nobs[i], height=bar_width[i],
                 color="green", label="used" if not labeled["sum_nobs"] else "")
        labeled["sum_nobs"] = True

        if show2:
            plt.barh(y + offsets[2] * bar_width[i], sum_nonobs[i], height=bar_width[i],
                     color="red", label="not used" if not labeled["sum_nonobs"] else "")
            labeled["sum_nonobs"] = True

    if radiance:
        plt.ylabel("Channel")
        plt.title("Observation Count vs Channel")
    else:
        plt.yscale("log")
        plt.gca().invert_yaxis()
        plt.ylabel("Pressure (hPa)")
        plt.title("Observation Count vs Pressure")

    xlabel = _xlabel_for(varname)
    if xlabel:
        plt.xlabel(xlabel)

    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.legend()

#Function to create the number of observations panel(only when comparing JEDI and GSI)
def comp_nobs_panel(varname, jsum_nobs, gsum_nobs, bin_centers, bin_heights, usrqc, radiance):
    bar_width = bin_heights * 0.4
    offsets = [-1.5, -0.5, 0.5, 1.5]

    labeled = {"JEDI nobs": False, "GSI nobs": False}
    for i in range(len(bin_centers)):
        y = bin_centers[i]

        plt.barh(y + offsets[0] * bar_width[i], jsum_nobs[i], height=bar_width[i],
                 color="green", label="JEDI nobs" if not labeled["JEDI nobs"] else "")
        labeled["JEDI nobs"] = True

        plt.barh(y + offsets[1] * bar_width[i], gsum_nobs[i], height=bar_width[i],
                 color="red", label="GSI nobs" if not labeled["GSI nobs"] else "")
        labeled["GSI nobs"] = True

    if not radiance:
        plt.yscale("log")
        plt.gca().invert_yaxis()
    xlabel = _xlabel_for(varname)
    if xlabel:
        plt.xlabel(xlabel)
    if radiance:
        plt.ylabel("Channel")
        plt.title("Observation Count vs Channel")
    else:
        plt.ylabel("Pressure (hPa)")
        plt.title("Observation Count vs Pressure")
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.legend()
