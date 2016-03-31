# Orbital-element-plotter

Codes for plotting the orbital elements of imaged companions, as functions of their unknown 
line of sight coordinates.

This repository contains three similar Python codes, which each perform the same function but
allow the user to input data in different formats. Given a companion's sky plane position at
two epochs, these programs produce contour plots of the companion's possible orbital elements 
as functions of its unknown line of sight (z, vz) coordinates at the first observation epoch.
For further details see Pearce, Wyatt & Kennedy 2015. The method assumes that the companion's
motion is close to linear between the two epochs, so it is best suited to analysis of long
-period companions. The three versions of the code allow the user to input the sky plane 
coordinates of the companion relative to the primary in one of three different formats: in 
terms of position angle and separation (orb_elmnt_contours_sep_PA.py), in terms of the 
companion's North and East offset from the primary (orb_elmnt_contours_N_E.py), and finally
in terms of the dimensionless B and phi parameters described in the above paper
(orb_elmnt_contours_B_phi.py). Instructions on how to use these codes are given in the 
comments at the top of each Python file. The codes produce Python plots of the contours, 
and also output all of the necessary data in .txt format for the user's plotting program
of choice.
