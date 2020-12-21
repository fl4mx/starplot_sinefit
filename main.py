import os
import math
import statistics
import numpy as np
import matplotlib.pyplot as plt

#data locations
LCdir = r"C:\Users\micha\Documents\OGLE-ATLAS-RR-c\sample_DATs"
stardata = r"C:\Users\micha\Documents\OGLE-ATLAS-RR-c\star_data.txt"

#get file directory
LCfilelist = os.listdir(LCdir)

#reading all star data from DAT file (names and periods)
names = np.loadtxt(stardata, dtype=str, skiprows=7, usecols=0)
periods = np.loadtxt(stardata, skiprows=7, usecols=8)



# generate sine function
def sine(average, amplitude, mid_phase, phases):
    sine_brightnesses = (amplitude * np.sin(((phases - mid_phase) * 2 * math.pi))) + average

    #calculate root mean square deviation from sine curve of LC
    rms_diff = math.sqrt((1 / len(sine_brightnesses)) * sum(abs(sine_brightnesses - brightnesses) ** 2))

    #jerry rigged solution to prevent erroneous fitting of the wrong midpoint (half-phase off)
    if rms_diff > 0.10:
        sine_brightnesses = (-1 * sine_brightnesses) + (2 * average)
    #calculate root mean square deviation from sine curve of LC (check, after jerry rigging)
    rms_diff = math.sqrt((1 / len(sine_brightnesses)) * sum(abs(sine_brightnesses - brightnesses) ** 2))

    return (sine_brightnesses, rms_diff)
    return rms_diff

#mid_phase +- value tester to determine optimization direction
def mid_phase_shift_test(mid_phase, shift):
    plusdelta = mid_phase + shift
    minusdelta = mid_phase - shift
    rms_diff_diff = sine(average, amplitude, plusdelta, phases)[1] - sine(average, amplitude, minusdelta, phases)[1]
    if rms_diff_diff < 0:
        direction = "right is better"
    else:
        direction = "left is better"
    return (direction, rms_diff_diff)



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
    sine_brightnesses, rms_diff = sine(average, amplitude, mid_phase, phases)

    #note how when mid_phase is increased,the curve is further shifted to the right. By testing two values, +- 0.1, we may find which direction of shifting the sine curve is better optimized, around the minmax point.
    direction, rms_diff_diff = mid_phase_shift_test(mid_phase, 0.1)


    #print star name
    print("name = " + str(name))
    #print RMS diff, smaller = better initial fit
    print("rms = " + str(rms_diff))
    #print optimization details (how much +- 0.1 was off of each other), smaller = better initial fit
    print("plusminusdelta diff = " + str(np.abs(rms_diff_diff)))
    #determine the optimization direction to shift the sine curve
    print("opti direction = " + str(direction))

    #mid_phase shift parameters go here:
    """
    ~good code~
    """

    #plot the LC, sine curve
    plt.scatter(phases, brightnesses)
    plt.xlabel("Phase")
    plt.ylabel("Magnitude")
    plt.title(name)
    plt.scatter(phases, sine_brightnesses)
    plt.show()
