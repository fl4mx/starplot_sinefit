#dependencies
import os
import math
import statistics
import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt

#data locations; replace with coordinates of LC DAT directory and metadata txt file, respectively
LCdir = r"C:\Users\micha\Documents\SciExt_RR-c\OGLE-ATLAS-RR-c\sample_DATs"
outputdir = r"C:\Users\micha\Documents\SciExt_RR-c\PROCESSED"
stardata = r"C:\Users\micha\Documents\SciExt_rr-c\OGLE-ATLAS-RR-c\star_data.txt"

#get file directory
LCfilelist = os.listdir(LCdir)

#reading all star data from DAT file (names and periods)
names = np.loadtxt(stardata, dtype=str, skiprows=7, usecols=0)
periods = np.loadtxt(stardata, skiprows=7, usecols=8)




#implementing Gauss-Newton Algorithm for curve fitting non linear regression; in this case, of the form Asin(Bx+C)+D
def gaussnewton(name, phases, brightnesses, average, amplitude, mid_phase):
    #Gauss-Newton iterations and damping parameters
    iter = 400
    damping = 0.1

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
    # print(x)
    # print(y)

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
    gaussnewton_sin = np.array([B[0] * np.sin(B[1] * x + B[2]) + B[3]])
    rms = math.sqrt((1 / (gaussnewton_sin.size)) * np.sum(abs(np.subtract(gaussnewton_sin, brightnesses)) ** 2))

    #print the RMS
    print("RMS deviation = " + str(rms))

    #return fitted curve and RMS
    return (gaussnewton_sin, str(rms))




#driver to iterate through all the stars and plot
for countLC, file in enumerate(LCfilelist):
    #reading LC data from LC files (dates and brightness)
    dates = np.loadtxt(LCdir + "\\" + file, delimiter=" ", usecols=0)
    brightnesses = np.loadtxt(LCdir + "\\" + file, delimiter=" ", usecols=1)
    #grabbing relevant star data for current star (name, starting time, period) from DAT file
    name = names[countLC]
    starting_date = dates[0]
    period = periods[countLC]

    #simple phasing calculation to convert dates to 0 to 1 of a complete phase
    phases = ((dates - starting_date) / period) % 1

    #determining approximate values to fit the sine curve (Amplitude, Average Height Shift, Phase Shift), for initial guess values of Gauss-Newton algorithm.
    #mean value taken as the arithmetic mean of all y-values
    average = statistics.mean(brightnesses)
    #better amplitude is the mean of the differences between average and max/min
    amplitude = statistics.mean([(max(brightnesses) - average), average - min(brightnesses)])
    #check for the closest value to the mean value of brightness (mid_brightness), find its corresponding x value (mid_index, mid_phase)
    brightness_diff = lambda checkbrightness: abs(checkbrightness - average)
    mid_brightness = min(brightnesses, key=brightness_diff)
    mid_index = np.where(brightnesses == mid_brightness)[0][0]
    mid_phase = phases[mid_index]

    #print star name
    print("name = " + str(name))

    #Gauss-Newton fit, return the y values of the fitted gauss-newton curve
    gaussnewton_sin, rms = gaussnewton(name, phases, brightnesses, average, amplitude, mid_phase)

    #plotting
    #basic setup
    fig = plt.figure(figsize = (10,5))
    fig.suptitle(str(name), fontsize = 14, fontweight = 'bold')
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top = 0.93)
    ax.set_xlabel("Phase")
    ax.set_ylabel("Magnitude")
    #plotting original LC, and fitted sine curve
    ax.scatter(phases, brightnesses)
    ax.scatter(phases, gaussnewton_sin)
    #plotting the RMS deviation onto the graph too
    ax.text(0.90, 0.01, ("RMS = " + str(rms)), verticalalignment = 'bottom', horizontalalignment = 'right', transform = ax.transAxes, color = 'purple', fontsize = 10)
    #show plot if testing in IDE
    #plt.show()
    #save plot
    plt.savefig((outputdir + "\\" + (name) + ".png"), format = "png")
