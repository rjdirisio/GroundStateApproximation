This code takes a DMC wave function and calculates the vibrational frequencies and intensities using the ground
state approximation outlined by McCoy et. al here (https://doi.org/10.1021/jp811352c).

Some of this is inhereted code from a previous member of our group, most of the infrastructure of this code was not written by me (DMCClusters.py, CalculateSpectrum.py, 
molecularInfo.py) however I have made serious modifications to basically all of the code, including vectorization and extensions like calculating a reduced-dimensional Hamiltonian and overlap matrix.  In the future, I may restructure this code to become more generalizable.

In addition, all the other .py files are my own work within the context of this original code.
