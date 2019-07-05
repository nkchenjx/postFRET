% Writen by Jixin Chen @ Ohio University 2016

clear all;
close all;

ResultingRate = struct('rate', [], 'rateconstantsfit1', [],  'rateconstantsfit2',[], 'rateconstantsfit3', [], 'rateconstantsfit4', []);

config.RemoveOnePointData = 'false'; % if state A-B-C and B only last one data point, B will be merged to A if 'true'

% config.FRETgap = 0.4; %gap between two adjecent states
config.numState = 3;  % total number of states
config.FRETStates = [0.9, 0.5, 0.1];  % the state FRET values
% rate = zeros(config.numState); % initial state coloum
% %   0   2to1 3to1 ...
% %  1to2   0  3to2 ...
% %  1to3 2to3  0   ...          % (finalstate, initial state)
% rate = [ 0   k21;
%          k12   0];
%{ 
config.rate = [ 0   10   1  10  10   0;
                1   0   10  10  10  10;
                10  10   0  0.5  3   10;
                5   2   10   0  10   5;
                10  10  10  10   0  10;
                0   10  10  10  10   0;];
%}
config.rate = [ 0  15  25;
                5  0   30;
                10 20  0];
    
% rateuplimit = 1/(config.bintime); %set the upper limit
% ratelowlimit = 1/(config.BleachLifetime); % set the lowwer limit

config.delta_t = 0.001; % time step (s). Simulate data with delta_t at 10x finer than bin time to represent the real.
config.bintime = 0.01; % (s)
config.binsize = round(config.bintime/config.delta_t); % data points per bin
config.numMol = 50; % default number of total molecules set to 100
config.BleachLifetime = 5; % (s), the everage life time of the dye before bleaching. The lifetime is set to be exponential decay
display(['total simulation time ~',num2str(config.numMol*config.BleachLifetime), ' s...']);
config.BleachtimeThres = 2; % (s), if the lifetime is less than this value, remove it.
% snr = 2; % signal to noise ratio before bin. give std 1/sqrt(snr).
config.totalsignal = 20*config.delta_t/0.001; % Total photon per measuring time, typical values we use is 20 for 1 ms data
config.backgroundnoise = 0.7/sqrt(config.delta_t/0.001); % Poisson distribution, roughly average 0.7-0.8 counts typically for 1 ms data
config.amplifiernoise = 0.37/sqrt(config.delta_t/0.001); % std(signal-mean)/mean, Gaussian, typical value 0.37+-0.03 for 1 ms data.
%%%%%%check this value again%%%%%%
% ampsnr = 1/amplifiernoise^2; % signal to noise ratio before bin. give std 1/sqrt(snr).
config.checkAddnoise = 'true';  % 'true' add noise with snr. 'false' no noise 
config.CheckFigure = 'true'; % 'true' testing figures are drawn. 'false' no figures are drawn.

save 'config';

ResultingRate = MCsimulation(config, config.rate); 

FileName=sprintf('Results_%s.mat',datestr(now,'yyyymmdd_HHMMSS'));
save(FileName);
%% END (by Jixin.Chen)
