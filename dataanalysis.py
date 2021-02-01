#dependencies
import os
import csv
import math
import statistics
import numpy as np
import matplotlib.pyplot as plt
from csv import writer

#data locations; replace with coordinates of LC DAT directory and metadata txt file, respectively
statsfile = r"C:\Users\micha\Documents\DataAnalysis\stats.csv"
metadata = r"C:\Users\micha\Documents\DataAnalysis\BLG_metadata.txt"
LCdir = r"C:\Users\micha\Documents\DataAnalysis\I"
outputdir = r"C:\Users\micha\Documents\DataAnalysis\Graphs"


#functions
def deviationvscolor(statsfile):
    # loadtxt read into numpy ndarrays
    deviations = np.loadtxt(statsfile, skiprows=0, delimiter=',', usecols=0)
    colors = np.loadtxt(statsfile, dtype="str", skiprows=0, delimiter=',', usecols=1)

    #colors[colors == "No color photometric data"] = 0
    for count, color in enumerate(colors):

        if color == "No color photometric data":
            #sentinel values for mask
            colors[count] = "99999"
            deviations[count] = 99999
        #print(count)
        #print(colors[count])
    colors = np.asfarray(colors ,float)
    colors = np.ma.masked_equal(colors, 99999)
    deviations = np.ma.masked_equal(deviations, 99999)
    scatterplot(deviations, colors, "RMS Deviation", "Bp_Rp Color", "RMS Deviation vs Color", 0.2)

def deviationvsbrightness(statsfile, metadata):
    # loadtxt read into numpy ndarrays
    deviations = np.loadtxt(statsfile, skiprows=0, delimiter=',', usecols=0)
    avgbrightnesses = np.loadtxt(metadata, skiprows=7, usecols=6)
    fallbackbrightnesses = np.loadtxt(metadata, skiprows=7, usecols=6)

    for count, brightness in enumerate(avgbrightnesses):
        if (brightness == "RRc"):
            avgbrightnesses[count] = fallbackbrightnesses[count]
    avgbrightnesses = np.ma.masked_equal(avgbrightnesses, -99.99)
    #plot
    scatterplot(deviations, avgbrightnesses, "RMS Deviation", "Avg Magnitude", "RMS Deviation vs Avg Magnitude", 0.1)
    #scatterplot(deviations, colors, "RMS Deviation", "Bp_Rp Color", "RMS Deviation vs Color", 0.3)


def periodvsmagnitude(metadata):
    #loadtxt read into numpy ndarrays
    periods = np.loadtxt(metadata, skiprows=7, usecols=8)
    fallbackperiods = np.loadtxt(metadata, skiprows=7, usecols=7)
    ibands = np.loadtxt(metadata, skiprows=7, usecols=6)
    for count, period in enumerate(periods):
        if period < 0.05:
            periods[count] = fallbackperiods[count]
#        if ibands[count] < 0:
#            print("less than 0")
#            print(ibands[count])
#            np.delete[ibands, count]
#            np.delete[periods, count]
    ibands = np.ma.masked_equal(ibands, -99.99)
    scatterplot(periods, ibands, "Period Length", "Mean I-Band Magnitude", "Period Lengths vs Mean Magnitude", 0.25)

def amountofdeviation():
    # loadtxt read into numpy ndarrays
    # deviation first column in CSV
    deviations = np.loadtxt(statsfile, skiprows=0, delimiter=',', usecols=0)
    #values
    #number = len(deviations)
    histplot(deviations, 20, "RMS Deviation", "Frequency", "RMS Deviation Frequency Histogram")

#plotters
def scatterplot(x, y, xtitle, ytitle, plotname, size):
    fig = plt.figure(figsize=(10, 5))
    fig.suptitle(str(plotname), fontsize=20, fontweight='bold', y = 0.96)
    ax = fig.add_subplot(111)
    fig.subplots_adjust(left=0.1, right=0.98, top=0.87, bottom=0.1)
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)
    #plot
    ax.scatter(x, y, s = size)
    # show plot if testing in IDE
    #plt.show()
    # save plot
    plt.savefig((outputdir + "\\" + (plotname) + ".png"), format="png")
    plt.close("all")

def histplot(x, bins, xtitle, ytitle, plotname):
    fig = plt.figure(figsize=(10, 5))
    fig.suptitle(str(
        plotname), fontsize=20, fontweight='bold', y=0.96)
    ax = fig.add_subplot(111)
    fig.subplots_adjust(left=0.1, right=0.98, top=0.87, bottom=0.1)
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)
    #plt.xlim(0.005, 0.15)
    #plt.xticks([0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21])
    plt.xticks(np.arange(0, 0.2, 0.02))
    # plot
    ax.hist(x, bins="auto")
    # show plot if testing in IDE
    plt.show()
    # save plot
    #plt.savefig((outputdir + "\\" + (plotname) + ".png"), format="png")
    plt.close("all")


#execution
periodvsmagnitude(metadata)
#amountofdeviation()
deviationvscolor(statsfile)
deviationvsbrightness(statsfile, metadata)
