-----------------------------------------------------------------------------------------
Coded with MATLAB 2014b

Sep. 2015

-----------------------------------------------------------------------------------------
Instructions for installation:

If you are using the MATLAB source, extract to any directory you like. 

----------------------------------------------------------------------------------------
Instructions for use:
1. postFRET.m
The code is set to guess 5-state rate constants from a set of input rate constant matrix. The parameters are to represent the simulation in the supporting information of the paper.

The code initialize parameters and then uses Monte Carlo simulation to guess the best fit for each rate constant. It loads the below functions

The results are shown in value ResultRate and ResultRateStd, showing the mean value of the best rate constants and the standard deviation of multiple simulations.

------------------------------------------------------------------------------------------
Functions:
2. findRateError.m
   Run a Monte Carlo simulation and find the average absolute error of the results from the target.

3. awgn.m
   Generate a Gaussian noise.

4. CountTransitions_difflength.m
   Rate constant analysis to generate a matrix of rate constants from the FRET trajectory. One-data-point reassignment is an option but set to 'false' by default.

5. exprnd.m
   Generate exponentially distributed random number.

6. findRateMin.m
   Compare the screened rate and find the one to give a minimum error score.

7. findRateMin2.m and findRateMin3.m
   Fit the error score vs rate to determine the rate that gives minimum error score.

8. poissrnd.m
   Generate a Poisson noise.

-------------------------------------------------------------------------------------------
The code is in five-state mode and one can change it to a system with any number of states by changing the parameters accordingly.


--------------------------
2015 by Dr. Jixin Chen @ Ohio University, Chemistry&Biochemistry.  chenj@ohio.edu
