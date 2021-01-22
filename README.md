# starplot_sinefit
<h2>Code for Michael Chen's 2020 - 2021 Science Extension Project
<h3>An Exploration of the Effects of the Deviation of RR-C Type Variables’ Light Curves from Sinusoidal Waveforms.</h3>
<p>RR-c Lightcurves (/mainfolder/OGLE-ATLAS-RR-c/LCdir) and metadata (/mainfolder/OGLE-ATLAS-RR-c/star_data) are all sourced from the OGLE-III Catalog of Variable Stars, available at http://ogledb.astrouw.edu.pl/~ogle/CVS/ <em>(Soszyński et al., 2009a, Soszyński et al., 2010c, Soszyński et al., 2011a)</em>.<br><br>A Sine Wave is fitted using the Gauss-Newton non-linear regression algorithm onto the Lightcurve, and calculates the RMS deviation of the curves from there.<br><br>The code also queries for the bp_rp color of the star using Astroquery, in the Gaia space telescope data, available at https://gea.esac.esa.int/archive/ <em>(Gaia Collaboration et al. 2016b, Gaia Collaboration et al., 2020b)</em>.<br><br>The fitted Sine Wave, original Lightcurve, RMS, RA/dec coordinates, and bp_rp color are all represented on a plot which is exported as a PNG to a sub-directory in the main folder, /mainfolder/PROCESSED/<br><br></p>
<h3 align="center">Example Plot</h3>
<p align="center"><img src="https://i.imgur.com/UosouhG.png" width="350" title="OGLE-BLG-RRLYR-00238"></p>
<p align="center">Example Plot of RR-C star, OGLE-BLG-RRLYR-01088</p>
