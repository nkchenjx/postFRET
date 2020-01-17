% The photocount-dependent Gaussian model
% by Jixin Chen@ Ohio University 1/16/2020
% Assume the sigma of the Gaussian is linearly dependent on the average
% photon counts of a state.

% trajIn is a noise free trajector
% noiseMdl is the equation to calculate the noise
% noiseP is the parameter used in the noiseMdl
% datatimestep is the time step resolution used to get the noise model
% simutimestep is the time step resolution of trajIn
% trajOut returns the trajectory with noise

function trajOut = noiseAdding(trajIn, noiseMdl, noiseP, datatimestep, simutimestep)


noiseSigma = noiseMdl(noiseP, trajIn);
noiseSigma = noiseSigma/sqrt(datatimestep/simutimestep); 
% for a faster time resolution, the signal drops 1/x but the noise/signal ratio goes up sqrt(x). 
% thus, the amplitude of noise drop sqrt(x).


trajOut = JCawgn(trajIn, noiseSigma);

end