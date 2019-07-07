-----------------------------------------------------------------------------------------
Version 2.0 
Coded with MATLAB 2014b
July, 2019
by
Jixin Chen @ Ohio University, Department of Chemistry and Biochemistry, Athens, Ohio 45701
chenj@ohio.edu

-----------------------------------------------------------------------------------------
Instructions for installation:

Extract to any directory you like and change the pathes in S1_loadData.m if needed. 


CITATION:

Jixin Chen,* Joseph R Pyle, Kurt Waldo Sy Piecco, Anatoly B Kolomeisky, Christy F Landes, A Two-Step Method for smFRET Data Analysis. J. Phys. Chem. B 2016, 120(29), 7128 - 7132 


METHOD:

The idea is to fit the experimental single-molecule FRET (smFRET, Wikipedia link) data once using any given method, used in the codes Ver 1.0 and Ver 2.0 the simplest thresholding method, i.e. set a threshold (e.g. the FRET value in the middel of two states) to distinguish two states. The analyzed results contains two major errors: (1) state miss-assignement due to the noise, (2) state miss-assignment due to camera blurring. Then simulate >hundreds of virtal data hope that one can find one of more trajectories that look just like the experimental data using the same analysing method, e.g. the thresholding method. Because we know the ground truth of the simulated data, we assume that the hidden truth of the real experimental data is the same as the simulated data that look the same (minimizing L1-norm, percentage error, as the judging standard in Ver 1.0 and 2.0).

The guesssing algorithm of the simulated data used in the codes is semi-exhausive searching algorithm called JCFit (GitHub link), a fitting algorithm that searches a parameter in an equation (model) within a defined boundary. The searching spacing is exponentially distributed away from the inital guess to the boundary. E.g. -10 to 10 are the boundaries and 1.0 is the initial guess, and 0.1 is the searching acuracy and ln(2) is the exponential factor, then the searching space is [1.0, 1.1, 1.3, 1.7, 2.5, 4.1, 7.3, 10] going up, and [1.0, 0.9, 0.7, 0.3, -0.5, -2.1, -5.3, -10] going down. The boundaries in version 2.0 is set mobile among searching iterations.

 

MATLAB CODES (Ver 2.0):

Step 1: S1_loadData.m

           Load raw data. Data is adapted from kinSoftChallenge 2019 with the following format: (1) each trace has been corrected from photoblinking, photobleaching, and backgrounds; (2) each trace has five columns, 1 time, 2 donor photocounts, 3 acceptor photocounts,  4 sum of donor and acceptor photocounts, 5 FRET values.

Step 2: S2_loadKey.m

           If the answer of the system is known, type in the answer, that is the transition rates between states. If not, give it a random guess based on the number of states one would like to fit (manually determined in postFRET).

Step 3: S3_findNoise.m

            Version 1.0 use a noise model, signal-dependent Gaussian noise pluss signal-independent Poisson noise. Version 2.0 use a signal-independent Gaussian noise model used in the kinSoftChallenge 2019. One has to custumize the noise model in findRateSimuNoise.m.

Step 4: S4_ConfigAndFRETAnalysis.m

            Set the fitting parameters and analyze the experimental data.

Step 5: S5_postFRET.m

            Simulate trajectories and find the best match.


--------------------------
2019 by Dr. Jixin Chen @ Ohio University, Chemistry&Biochemistry.  chenj@ohio.edu
