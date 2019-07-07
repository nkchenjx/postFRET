% 2016 by 
% Dr. Jixin Chen 
% Department of Chemistry and Biochemistry
% Ohio University,
% Athens, Ohio 45701 
% chenj@ohio.edu
% 
%    clear all;
%    close all;

figure; hist(data(:,5),100);
prompt = {'masually tell the number of state:'};
dlgtilte = 'numstate';
dims = [1 35];
definput = {'2'};
answer = inputdlg(prompt,dlgtilte,dims,definput);
numstate = str2double(answer{1});

prompt = {'state 1', 'state 2', 'state 3', 'state 4'};
dlgtilte = 'state value, 0 means no state';
dims = [1 35];
definput = {'0', '0', '0', '0'};
answer = inputdlg(prompt,dlgtilte,dims,definput);
state1 = str2double(answer{1});
state2 = str2double(answer{2});
state3 = str2double(answer{3});
state4 = str2double(answer{4});


%% initialize config
%config.rateTrue = config.rateTrue;
config.RemoveOnePointData = 'true'; % if state A-B-C and B only last one data point, B will be merged to A if 'true'
config.WAMSearching = 'false'; % sequential searching the rate constant or search the state that gives the maximum wL1_SE score (WAM).
config.WL1thresh = 0.01; % the convergence wL1 score (WL1thresh*100)% for a "whack-a-mole" search to stop.  
                         % This value adjusts final accuracy. 
                         % This value has different limits for different conditions. 
                         % Ctrl+C to break during the simulation if it takes too long.                      
  display(['Set the wL1_SE target = ', num2str(config.WL1thresh)]);
  
  config.BleachLifetime = mean(config.bleachtime);

  
  config.numState = numstate;  % total number of states
%   0   2to1 3to1 ...
%  1to2   0  3to2 ...
%  1to3 2to3  0   ...          % (finalstate, initial state)

% config.rateTarget = [ 0     1.8;  % get later from the initial fitting
%                       4.6   0  ];

 switch numstate
        case 2
           config.FRETStates = [state1, state2];
        case 3
           config.FRETStates = [state1, state2, state3];
        case 4
           config.FRETStates = [state1, state2, state3, state4];
end


config.delta_t = config.datatimestep/10; % time step (s). Simulate data with delta_t at 10x finer than bin time to represent the real.
config.bintime = config.datatimestep; % (s)
config.binsize = round(config.bintime/config.delta_t); % data points per bin
config.numMol = length(datalength); % default number of total molecules set to 100
config.totalsignal = config.sumIaa*config.delta_t/config.bintime; % Total photon per simulated time step
   display(['simulation time of each run ~',num2str(config.numMol*config.BleachLifetime), ' s data points...']);
config.BleachtimeThres = min(config.bleachtime); % (s), if the lifetime is less than this value, remove it.


% snr = 2; % signal to noise ratio before bin. give std 1/sqrt(snr).
% config.backgroundnoise = 0.7/sqrt(config.delta_t/0.001); % Poisson distribution, roughly average 0.7-0.8 counts typically for 1 ms data
% config.amplifiernoise = 0.37/sqrt(config.delta_t/0.001); % std(signal-mean)/mean, Gaussian, typical value 0.37+-0.03 for 1 ms data.

% config.backgroundnoise = 0; % Poisson distribution, roughly average 0.7-0.8 counts typically for 1 ms data
% config.noise = config.rawnoise*sqrt(config.binsize);
% config.amplifiernoise = config.sumnoise/config.sumIaa*sqrt(config.binsize); % std(signal-mean)/mean, Gaussian, typical value 0.37+-0.03 for 1 ms data.

  x = [-0.5, -0.2, -0.1, -0.045, -0.02, -0.005, 0, 0.006, 0.03, 0.05 0.11, 0.21, 0.51]; % for multiple factor for each rate constant  
  y = exp(x);
config.rateMulti = y ; % times the input value as the guess. 0.1 to 10
config.numoftry = length(config.rateMulti); % the length of the above vector.
numState = config.numState;
config.normFreq = ones(numState,numState);
  for i = 1:numState
      for j = 1:numState
          config.normFreq(i,j) = ceil(abs(i-j));
          if j==i;
              config.normFreq(i,j) =1;
          end
      end
  end
  
    %set energy barier for WAM state searching. The probability of
%searching a non-maximum state is   exp(-EnergyBarrier/n), where n is the number of times before this search for the same state. 
config.energyBarrier = 4; 

config.CheckFigure = 'false'; % 'true' testing figures are drawn. 'false' no figures are drawn.
config.CheckSaveTraj = 'false'; % if 'true', the experimental trajectory will be saved into a matlab file.

% end of initializing parameters


% data2 = data(:,[1,2,3,5]);

% analyze raw data for the target rate
[rateconstantsfit, numtransitionfit, ratesum] =  CountTransitions_difflength(data(:,[1,2,3,5]),datalength, config.FRETStates,config.RemoveOnePointData);
config.rateTarget = rateconstantsfit;
%--------------------over Writing the single simulation results to average
%results of five simulations:
% config.rateTarget = [0,11.0490484041813,6.99900675780571,4.62117144460293,3.12378208684813,1.77988941058650;18.5350496509884,0,14.4624683203641,9.51391094906703,5.67664047804632,4.14623983875533;14.1908731764572,19.2691468377841,0,17.5685738498993,11.7096255248044,7.99941682980906;9.16397389645158,12.7451014313252,18.4262656751761,0,19.8864730421171,15.3623301340084;4.91208151300297,6.05741911733557,9.27719742196829,14.9463379325173,0,19.4879063883545;1.76809060166991,2.64680105448343,4.00628653127624,6.89471292033787,11.0523423737904,0];
%--------------------

fprintf('\n\n--------------------------------------------------\n');
fprintf('True rates =\n');
disp(num2str(config.rateTrue));

fprintf('\nExperimental rates =\n');
disp(num2str(config.rateTarget));

fprintf('\nwL1_AT score for raw data = ');
wL1 = findWL1(config.rateTarget, config.rateTrue);
fprintf([num2str(wL1),'\n']);
fprintf('\n--------------------------------------------------\n\n');
