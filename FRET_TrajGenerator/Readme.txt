-----------------------------------------------------------------------------------------
Coded with MATLAB 2014b

2016

-----------------------------------------------------------------------------------------
Instructions for installation:

If you are using the MATLAB source, extract to any directory you like. 

Cite:
A Two-Step Method for smFRET Data Analysis
Jixin Chen, Joseph R. Pyle, Kurt Waldo Sy Piecco, Anatoly B. Kolomeisky, and Christy F. Landes
The Journal of Physical Chemistry B 2016 120 (29), 7128-7132
DOI: 10.1021/acs.jpcb.6b0569

----------------------------------------------------------------------------------------
Instructions for use:
1. FRETTrajGenerator.m
 This code generate smFRET trajectory and analyze the resulting rate constant based on thresholding method. The resulting rate constants is in structure "ResultingRate" and saved to "Results_xxx.mat".
 The bin-free and noise-free trajectories have been saved in separated files if the config.CheckFigure is set to 'true'. The smFRET trajectories are saved in for "DataXXX_XXX.mat" files.
 The correction method one-point-data reassignment is by default set to be 'false'.


------------------------------------------------------------------------------------------
Functions:
2. addAPDNoise.m
   Add noise to a FRET trajectory.

3. JCawgn.m
   Generate a Gaussian noise.

4. binning.m
   Bin the FRET trajetory.

5. CountTransitions_difflength.m
   Rate constant analysis to generate a matrix of rate constants from the FRET trajectory.

6. exprnd.m
   Generate exponentially distributed random number.

7. MCsimulation.m
   Run the Monte Carlo Simulation and give the resulting rate constants. Time trajectories are saved. 

7. poissrnd.m
   Generate a Poisson noise.

-------------------------------------------------------------------------------------------
The code is in 3-state mode and one can change it to a system with any number of states by changing the parameters accordingly.

----------------------
2016 by Dr. Jixin Chen @ Ohio University, Chemistry&Biochemistry.  chenj@ohio.edu
