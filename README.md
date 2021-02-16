# starplot_sinefit
<h2>Michael Chen 2020 - 2021
<h3>An Exploration of the Effects of the Deviation of RR-C Type Variables’ Light Curves from Sinusoidal Waveforms.</h3>
<p>RR-c Lightcurves (/folder/OGLE-ATLAS-RR-c/LCdir) and metadata (/folder/OGLE-ATLAS-RR-c/star_data) are all sourced from the OGLE-III Catalog of Variable Stars, available at http://ogledb.astrouw.edu.pl/~ogle/CVS/ <em>(Soszyński et al., 2009a, Soszyński et al., 2010c, Soszyński et al., 2011a)</em>.<br><br>A Sine Wave is fitted using the Gauss-Newton non-linear regression algorithm onto the Lightcurve, and calculates the RMS deviation of the curves from there.<br><br>The code also queries for the bp_rp color of the star using Astroquery, in the Gaia space telescope data, available at https://gea.esac.esa.int/archive/ <em>(Gaia Collaboration et al. 2016b, Gaia Collaboration et al., 2020b)</em>.<br><br>The fitted Sine Wave, original Lightcurve, RMS, RA/dec coordinates, and bp_rp color, from the main driver script (main.py) are all represented on a plot which is exported as a PNG to a sub-directory in the main folder, (/folder/PROCESSED/)<br><br>The plots from the data analysis and plotter script (dataanalysis.py) are represented below as well, and are exported to a sub-directory in its own main folder, (/folder/Graphs)<br><br></p>
<h3 align="center">Example Plot</h3>
<p align="center"><img src="https://i.imgur.com/UosouhG.png" width="350" title="OGLE-BLG-RRLYR-00238"></p>
<p align="center">Example Plot of RR-C star, OGLE-BLG-RRLYR-01088</p>
<br>
<h3 align="center">Data Analyses</h3>
<p align="center"><img src="https://i.imgur.com/7SV1o0y.png" width="350" title="Period Lengths vs Mean Magnitude"></p>
<p align="center">Generated plot of Period Lengths vs Mean Magnitude</p>

<p align="center"><img src="https://i.imgur.com/Zry0rEj.png" width="350" title="RMS Deviation Frequency Histogram"></p>
<p align="center">Generated plot of RMS Deviation Frequency Histogram</p>

<p align="center"><img src="https://i.imgur.com/X6458aU.png" width="350" title="RMS Deviation vs Color"></p>
<p align="center">Generated plot of RMS Deviation vs Color</p>

<p align="center"><img src="https://i.imgur.com/XrUagjy.png" width="350" title="RMS Deviation vs Mean Magnitude"></p>
<p align="center">Generated plot of RMS Deviation vs Mean Magnitude</p>
