#dependencies
import os
from csv import writer
import math
import statistics
import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt
import astropy.units as u
from astropy import coordinates as coords
from astroquery.gaia import Gaia

#data locations; replace with coordinates of LC DAT directory and metadata txt file, respectively
LCdir = r"C:\Users\micha\Documents\BLG_data_processing\OGLE-ATLAS-RR-c\I"
outputdir = r"C:\Users\micha\Documents\BLG_data_processing\PROCESSED"
stardata = r"C:\Users\micha\Documents\BLG_data_processing\OGLE-ATLAS-RR-c\BLG_metadata.txt"
lastsortedloc = r"C:\Users\micha\Documents\BLG_data_processing\OGLE-ATLAS-RR-c\processed.txt"
colordeviationstats = r"C:\Users\micha\Documents\BLG_data_processing\OGLE-ATLAS-RR-c\stats.csv"
debugloc = r"C:\Users\micha\Documents\BLG_data_processing\OGLE-ATLAS-RR-c\debug.txt"

#get file directory
LCfilelist = os.listdir(LCdir)
LCfilelength = len(LCfilelist)


#reading all star data from DAT file (names and periods for wave fitting, RA and dec for Gaia query)
names = np.loadtxt(stardata, dtype=str, skiprows=7, usecols=0)
periods = np.loadtxt(stardata, skiprows=7, usecols=8)
allRA = np.loadtxt(stardata, dtype = "str", skiprows = 7, usecols = 3)
alldec = np.loadtxt(stardata, dtype = "str", skiprows = 7, usecols = 4)

#lastsaved
lastsorted = open(lastsortedloc, "r")
laststar = int(lastsorted.read())
lastsorted.close()





#implementing Gauss-Newton Algorithm for curve fitting non linear regression; in this case, of the form Asin(Bx+C)+D
def gaussnewton(phases, brightnesses, average, amplitude, mid_phase):
    #in case of Gauss-Newton fitting error, we may fallback onto initial guess values for fit
    fallback_sin = (amplitude * np.sin(((phases - mid_phase) * 2 * math.pi))) + average
    fallback_rms = math.sqrt((1 / (fallback_sin.size)) * np.sum(abs(np.subtract(fallback_sin, brightnesses)) ** 2))
    #in case of phase flipped midpoint mis-fit
    phaseflip_fallback_sin = (-1 * fallback_sin) + (2 * average)
    phaseflip_rms = math.sqrt((1 / (phaseflip_fallback_sin.size)) * np.sum(abs(np.subtract(phaseflip_fallback_sin, brightnesses)) ** 2))

    if fallback_rms > phaseflip_rms:
        fallback_sin = phaseflip_fallback_sin
        fallback_rms = phaseflip_rms


    #Gauss-Newton iterations and damping parameters
    iter = 400
    damping = 0.02

    #PDEs in Gauss-Newton
    def pdA(x, b, c):
        return np.sin(b * x + c)
    def pdB(x, a, b, c):
        return a * x * np.cos(b * x + c)
    def pdC(x, a, b, c):
        return a * np.cos(b * x + c)
    def pdD(x):
        return 1

    #least squares
    def leastSquares(x, y, a, b, c, d):
        return y - (a * np.sin(b * x + c) + d)

    #standard method io
    x = phases
    y = brightnesses

    #initial guesses for A, B, C, D in sine
    B = np.matrix([[amplitude], [2 * math.pi], [(mid_phase) * 2 * math.pi], [average]])
    #jacobian matrix for diff eq
    J = np.zeros((x.size, 4))
    #least square vector
    r = np.zeros((x.size, 1))

    for _ in range(0, iter):
        for i in range(0, x.size):
            #compute each value of r in this iteration
            r[i, 0] = leastSquares(x[i], y[i], B[0], B[1], B[2], B[3])

            #calculate the Values of Jacobian matrix on this iteration
            J[i, 0] = pdA(x[i], B[1], B[2])
            J[i, 1] = pdB(x[i], B[0], B[1], B[2])
            J[i, 2] = pdC(x[i], B[0], B[1], B[2])
            J[i, 3] = pdD(x[i])

        Jt = J.T
        B += damping * np.dot(np.dot(inv(np.dot(Jt, J)), Jt), r)

    #test code to compare fitting performance between forced 2pi (full phase plot) vs B[1] (gauss newton determined % of plot)
    #2pisin = np.array([B[0] * np.sin(2 * math.pi * x + B[2]) + B[3]])
    #2pirms = math.sqrt((1 / (ysin.size)) * np.sum(abs(np.subtract(ysin, y)) ** 2))

    #generate Gauss-Newton fitted sine curve, and calculate RMS
    gaussnewton_sin = np.array([B[0] * np.sin((B[1] * x) + B[2]) + B[3]])
    rms = math.sqrt((1 / (gaussnewton_sin.size)) * np.sum(abs(np.subtract(gaussnewton_sin, brightnesses)) ** 2))

    #fallback in case of fitting failure, use initial values
    if rms > fallback_rms:
        print("gaussnewton failed, fallback to standard")
        gaussnewton_sin = fallback_sin
        rms = fallback_rms
        #read Gauss-Newton failrate
        debugfiler = open(debugloc, "r")
        GNfail = int(lastsorted.read())
        debugfiler.close()
        #write Gauss-Newton fail to debug
        debugfilew = open(debugloc, "w")
        debugfilew.write(str(GNfail + 1))
        debugfilew.close()

    #shorten RMS for better visibility
    roundedrms = float(rms)
    roundedrms = round(roundedrms, 4 - int(math.floor(math.log10(abs(roundedrms)))) - 1)

    #return fitted curve and RMS
    return (gaussnewton_sin, str(rms), str(roundedrms))





#implementing method of querying Gaia data for the bp-rp color of the star
def gaiaquery(starnumber, allRA, alldec):
    #reading coords of the star
    RA = allRA[starnumber]
    dec = alldec[starnumber]

    #setting up coords query, height and width search precision
    coord = coords.SkyCoord(ra = RA, dec = dec, unit = (u.hourangle, u.deg), frame = "icrs")
    height = u.Quantity(1, u.arcsec)
    width = u.Quantity(1, u.arcsec)

    #query
    star = Gaia.query_object(coordinate=coord, width=width, height=height, columns=["source_id, ra, dec, bp_rp"])
    #star is a table
    #print("star:")
    #print(star)
    if (star["bp_rp"].size == 0):
        color = "No color photometric data"
        roundedcolor = color
    else:
        color = str(star["bp_rp"][0])
        if color == "--":
            color = "No color photometric data"
            roundedcolor = color
        else:
            #shorten color for better visibility
            roundedcolor = float(color)
            roundedcolor = round(roundedcolor, 5 - int(math.floor(math.log10(abs(roundedcolor))) - 1))


    #return coordinates, color
    return (RA, dec, color, roundedcolor)





def plot(phases, brightnesses, gaussnewton_sin, RA, dec, rms, color, name, outputdir):
    fig = plt.figure(figsize = (10,5))
    fig.suptitle(str(name), fontsize = 22, fontweight = 'bold')
    ax = fig.add_subplot(111)
    fig.subplots_adjust(left = 0.1, right = 0.98, top = 0.87, bottom = 0.1)
    ax.set_xlabel("Phase")
    ax.set_ylabel("Magnitude")
    #plotting original LC, and fitted sine curve
    ax.scatter(phases, brightnesses)
    ax.scatter(phases, gaussnewton_sin)
    #adding legend to graph for readability
    ax.legend(["Original LC", "Fitted Sine Curve"])
    #plotting the RMS deviation onto the graph
    ax.text(0.94, 1.03, ("RMS = " + str(rms)), verticalalignment = 'bottom', horizontalalignment = 'center', transform = ax.transAxes, color = 'purple', fontsize = 10)
    #plotting coords, color onto the graph
    ax.text(0.01, 1.09, ("RA = " + str(RA)), verticalalignment = 'bottom', horizontalalignment = 'center', transform = ax.transAxes, color = 'green', fontsize = 8)
    ax.text(0.01, 1.06, ("dec = " + str(dec)), verticalalignment = 'bottom', horizontalalignment = 'center', transform = ax.transAxes, color = 'green', fontsize = 8)
    if color == "No color photometric data":
        ax.text(0.01, 1.03, (str(color)), verticalalignment = 'bottom', horizontalalignment = 'center', transform = ax.transAxes, color = 'orange', fontsize = 7)
    else:
        ax.text(0.01, 1.03, ("Color = " + str(color)), verticalalignment = 'bottom', horizontalalignment = 'center', transform = ax.transAxes, color = 'orange', fontsize = 8)
    #show plot if testing in IDE
    #plt.show()
    #save plot
    plt.savefig((outputdir + "\\" + (name) + ".png"), format = "png")
    plt.close("all")





#driver to iterate through all the stars, starting from where we last left off
for countLC in range(laststar, LCfilelength, 1):
    #specifying the LC data file, without iterating through listdir in outer for loop
    file = LCfilelist[countLC]
    #reading LC data from LC files (dates and brightness)
    dates = np.loadtxt(LCdir + "\\" + file, delimiter=" ", usecols=0)
    brightnesses = np.loadtxt(LCdir + "\\" + file, delimiter=" ", usecols=1)
    #grabbing relevant star data for current star (name, starting time, period) from DAT file
    name = names[countLC]
    starting_date = dates[0]
    period = periods[countLC]

    #simple progress indicator, tells us which star program is up to
    print(name)

    #simple phasing calculation to convert dates to 0 to 1 of a complete phase
    phases = ((dates - starting_date) / period) % 1

    #determining approximate values to fit the sine curve (Amplitude, Average Height Shift, Phase Shift), for initial guess values of Gauss-Newton algorithm.
    #mean value taken as the arithmetic mean of all y-values
    average = (max(brightnesses) + min(brightnesses))/2
    #better amplitude is the mean of the differences between average and max/min
    amplitude = statistics.mean([(max(brightnesses) - average), average - min(brightnesses)])
    #check for the closest value to the mean value of brightness (mid_brightness), find its corresponding x value (mid_index, mid_phase)
    brightness_diff = lambda checkbrightness: abs(checkbrightness - average)
    mid_brightness = min(brightnesses, key=brightness_diff)
    mid_index = np.where(brightnesses == mid_brightness)[0][0]
    mid_phase = phases[mid_index]

    #Gauss-Newton fit, return the y values of the fitted gauss-newton curve
    gaussnewton_sin, rms, roundedrms = gaussnewton(phases, brightnesses, average, amplitude, mid_phase)

    #query Gaia for color
    RA, dec, color, roundedcolor = gaiaquery(countLC, allRA, alldec)

    #temp star stats data to write to CSV later
    tempstatsarray = [rms, color]

    #plotting
    #basic setup
    plot(phases, brightnesses, gaussnewton_sin, RA, dec, roundedrms, roundedcolor, name, outputdir)

    #autosaver and resume
    lastsorted = open(lastsortedloc, "w")
    lastsorted.write(str(countLC + 1))
    lastsorted.close()

    #star stats for final processing
    with open(colordeviationstats, "a+", newline = "") as statsfile:
        csv_writer = writer(statsfile)
        csv_writer.writerow(tempstatsarray)
