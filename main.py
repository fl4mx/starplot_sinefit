import os
import math
import statistics
import numpy as np
import matplotlib.pyplot as plt

#data locations
LCdir = r"C:\Users\micha\Documents\OGLE-ATLAS-RR-c\sample_DATs"
stardata = r"C:\Users\micha\Documents\OGLE-ATLAS-RR-c\star_data.txt"

#reading all star data from DAT file (names and periods)
names = np.loadtxt(stardata, dtype=str, skiprows=7, usecols=0)
periods = np.loadtxt(stardata, skiprows=7, usecols=8)

LCfilelist = os.listdir(LCdir)
#iterate through star LC
for countLC, file in enumerate(LCfilelist):
    #reading LC data from LC files (dates and brightness)
    dates = np.loadtxt(LCdir + "\\" + file, delimiter=" ", usecols=0)
    brightnesses = np.loadtxt(LCdir + "\\" + file, delimiter=" ", usecols=1)
    #grabbing relevant star data for current star (name, starting time, period) from DAT
    name = names[countLC]
    starting_date = dates[0]
    period = periods[countLC]

    #simple phasing calculation
    phases = ((dates - starting_date) / period) % 1
    #mean value taken as the arithmetic mean of all y-values
    average = statistics.mean(brightnesses)
    #better amplitude is the mean of the differences between average and max/min
    amplitude = statistics.mean([(max(brightnesses) - average), average - min(brightnesses)])
    #check for the closest value to the mean value of brightness (mid_brightness), find its corresponding x value (mid_index, mid_phase)
    brightness_diff = lambda checkbrightness: abs(checkbrightness - average)
    mid_brightness = min(brightnesses, key=brightness_diff)
    mid_index = np.where(brightnesses == mid_brightness)[0][0]
    mid_phase = phases[mid_index]

    #trig transform onto the plot
    sine_brightnesses = (amplitude * np.sin(((phases - mid_phase)* 2 * math.pi))) + average

    #calculate root mean square deviation from sine curve of LC
    rms_diff = math.sqrt((1 / len(sine_brightnesses)) * sum(abs(sine_brightnesses - brightnesses) ** 2))

    #jerry rigged solution to prevent erroneous fitting of the wrong midpoint (half-phase off)
    if rms_diff > 0.10:
        sine_brightnesses = (-1 * sine_brightnesses) + (2 * average)

    #calculate root mean square deviation from sine curve of LC (check, after jerry rigging)
    rms_diff = math.sqrt((1 / len(sine_brightnesses)) * sum(abs(sine_brightnesses - brightnesses) ** 2))
    #print RMS diff
    print(str(rms_diff))

    #plot the LC, sine curve
    plt.scatter(phases, brightnesses)
    plt.xlabel("Phase")
    plt.ylabel("Magnitude")
    plt.title(name)
    plt.scatter(phases, sine_brightnesses)
    plt.show()
