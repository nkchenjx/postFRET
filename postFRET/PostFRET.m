% 2016 by 
% Dr. Jixin Chen 
% Department of Chemistry and Biochemistry
% Ohio University,
% Athens, Ohio 45701 
% chenj@ohio.edu

   clear all;
   close all;

fprintf('\n------------------------------------------------\n');

config.RemoveOnePointData = 'true'; % if state A-B-C and B only last one data point, B will be merged to A if 'true'
config.WAMSearching = 'false'; % sequential searching the rate constant or search the state that gives the maximum wL1_SE score (WAM).
config.WL1thresh = 0.04; % the convergence wL1 score (WL1thresh*100)% for a "whack-a-mole" search to stop.  
                         % This value adjusts final accuracy. 
                         % This value has different limits for different conditions. 
                         % Ctrl+C to break during the simulation if it takes too long.                      
  display(['Set the wL1_SE target = ', num2str(config.WL1thresh)]);
config.maxIteration = 500; %set the maximum iterations in a searching. (disabled now)


config.numState = 3;  % total number of states
config.FRETStates = [0.90, 0.5, 0.1];  % the state FRET values
config.BleachLifetime = 500; % (s), the everage life time of the dye before bleaching. The lifetime is set to be exponential decay
config.repeatTime = 3; % for each run repeat times to get the average

%   0   2to1 3to1 ...
%  1to2   0  3to2 ...
%  1to3 2to3  0   ...          % (finalstate, initial state)
config.ratePreset = [0 15 25;
                     5  0 30;
                     10 20 0];

config.delta_t = 0.001; % time step (s). Simulate data with delta_t at 10x finer than bin time to represent the real.
config.bintime = 0.01; % (s)
config.binsize = round(config.bintime/config.delta_t); % data points per bin
config.numMol = 1; % default number of total molecules set to 100
config.totalsignal = 20*config.delta_t/0.001; % Total photon per measuring time, typical values we use is 20 for 1 ms data
   display(['simulation time of each run ~',num2str(config.numMol*config.BleachLifetime), ' s data points...']);
config.BleachtimeThres = 2; % (s), if the lifetime is less than this value, remove it.

% set the scale: rate = rateTarget*rateMulti
% x = [-1.5, -1, -0.7, -0.5:0.1:-0.3, -0.2:0.05:0.2, 0.3:0.1:0.5, 0.7, 1, 1.5]; % for multiple factor for each rate constant  
  x = [-0.5, -0.2, -0.1, -0.045, -0.02, 0.03, 0.05 0.11, 0.21, 0.51]; % for multiple factor for each rate constant  
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

%% Generate experimental trajectory
fprintf('Generate the experimental trajectory.\n   Preset rates =\n');
disp(num2str(config.ratePreset));

rateTrue = zeros(numState, numState, config.repeatTime);
rateTarget = zeros(numState, numState, config.repeatTime);
for k = 1:config.repeatTime
    fprintf('.');
    [rateTrue(:,:,k), rateTarget(:,:,k)] = MCsimulation(config, config.ratePreset);
end
config.rateTrue = mean(rateTrue,3);
config.rateTarget = mean(rateTarget,3);
fprintf('Generate the experimental trajectory.\n   True rates =\n');
disp(num2str(config.rateTrue));
fprintf('\nAnalyze the rate constants of the binned trajectory.\n   Experimental rates =\n');
disp(num2str(config.rateTarget));

%--------------------over Writing the single simulation results to average
%results of five simulations:
% config.rateTarget = [0,11.0490484041813,6.99900675780571,4.62117144460293,3.12378208684813,1.77988941058650;18.5350496509884,0,14.4624683203641,9.51391094906703,5.67664047804632,4.14623983875533;14.1908731764572,19.2691468377841,0,17.5685738498993,11.7096255248044,7.99941682980906;9.16397389645158,12.7451014313252,18.4262656751761,0,19.8864730421171,15.3623301340084;4.91208151300297,6.05741911733557,9.27719742196829,14.9463379325173,0,19.4879063883545;1.76809060166991,2.64680105448343,4.00628653127624,6.89471292033787,11.0523423737904,0];
%--------------------

%% postFRET Monte Carlo simulating to guess a trajectory that matches the experimental trajectory the best
   fprintf('\n------------------------------------------------\n');
fprintf('\npostFRET Monte Carlo simulating to guess a trajectory that matches the experimental trajectory the best.\n');
tic

numState = config.numState;
% Iteratively screen each rate constant to find the best fit
if strcmp(config.WAMSearching, 'true')
%----------Whack-A-Mole searching------------
    initialGuess = config.rateTarget; 
    config.rate(:,:,1) = initialGuess;
    % minScore = findRateError(config, initialGuess);
    bestGuess = initialGuess; 
    rateScore =[];  
    rateSimu = zeros(numState, numState, config.repeatTime);
    for k = 1:config.repeatTime
          rateSimu(:,:,k) = findRateSimu(config, bestGuess);
          wl1Score(k) = findWL1(rateSimu(:,:,k),config.rateTarget);
    end
    wl1Scoremean = mean(wl1Score);
    minScore = wl1Scoremean;
    rateSimumean = mean(rateSimu,3);
    [mf,mi] = findIndMaxDiff(rateSimumean, config.rateTarget);
    rateHistory = rateSimumean;
    wl1ScoreHistory = minScore;
    scanHistory = [mi,mf];
    scanFrequency = ones(numState, numState);
       for i = 1:numState
           scanFrequency(i,i) = 0;
       end
    scanFrequency(mf,mi) = scanFrequency(mf,mi)+1;
    

    fprintf(['Initial wL1_SE score = ', num2str(minScore) ,' \n']);
    while(minScore>config.WL1thresh)
        %search all transitions invove these two states for a better guess
%        if wl1Scoremean <= 0.05
%            config.BleachLifetime = 5000;
%            config.repeatTime = 1;
%        elseif wl1Scoremean <= 0.8
%            config.BleachLifetime = 3000;
%            config.repeatTime = 1;
%        elseif wl1Scoremean <= 0.1
%            config.BleachLifetime = 1500;
%            config.repeatTime = 2;
%        elseif wl1Scoremean <= 0.2
%            config.BleachLifetime = 1000;
%            config.repeatTime = 1;
%        elseif wl1Scoremean <= 0.25
%            config.BleachLifetime = 200;
%            config.repeatTime = 1;
%        else 
%            config.BleachLifetime = 100; 
%            config.repeatTime = 1;
%        end
        rate = bestGuess;
        rateScore =[];
        rateScan = [];
        rateScanSimu = [];

       fprintf(['Scanning transtion ', num2str(mi),'  to ', num2str(mf), ' (column to row) \n']);
       for j = 1:config.numoftry
            rate(mf,mi) = config.rateMulti(j).*bestGuess(mf,mi);
            for k = 1:config.repeatTime
              fprintf('.');
              rateSimu(:,:,k) = findRateSimu(config, rate);
              wl1Score(k) = findWL1(rateSimu(:,:,k), config.rateTarget);
            end
            wl1Scoremean = mean(wl1Score);
            rateSimumean = mean(rateSimu,3);
            rateScan{j} = rate;
            rateScanSimu{j} = rateSimumean;
            rateScore(j) = wl1Scoremean;
        end
        [minScore,ind] = min(rateScore);
          fprintf(['\n score = ', num2str(minScore), '\n']);
        bestGuess = rateScan{ind};
        disp(num2str(bestGuess));
        bestGuessSimu = rateScanSimu{ind};
        rateHistory = cat(3, rateHistory, bestGuess);
        wl1ScoreHistory = [wl1ScoreHistory, minScore];

      %find a state to check and decide if it has been checked too many times.
        [mf,mi] = findIndMaxDiff(bestGuessSimu, config.rateTarget); %find which rate constant is mostly off.
%         l = size(scanHistory, 1); %iterations
%         count = 1;
%         for i = l:1
%             if sum(scanHistory(i,:) == [mf,mi]) == 2
%                 count = count+1;
%             else
%                 break
%             end
%         end
%         if rand<=exp(-config.energyBarrier/count)
%             mi = ceil(rand*numState);
%             mi(mi==0) = 1;
%             mf = ceil(rand*numState);
%             while (mf == mi)
%                 mf = ceil(rand*numState);
%             end
%             mf(mf==0) = 1;
%               fprintf(['\n !!!!!A random jump to transtion ', num2str(mi),' to ', num2str(mf), '\n']);
%         end
       % if a state has not been scanned for long time, it has the chance to be scanned specially. 
        tempsf = scanFrequency./config.normFreq;
        for i = 1:numState
            tempsf(i,i) = inf;
        end
        f = min(tempsf(:));
        for i = 1:numState
            tempsf(i,i) = 0;
        end
        if f~=tempsf(mf,mi) && rand<=exp(-config.energyBarrier/(max(tempsf(:)-f)))
              [a, b] = find(tempsf == f);
              i = ceil(rand*length(a));
              mf = a(i);
              mi = b(i);
        end
        scanHistory = [scanHistory; mi, mf];
        scanFrequency(mf,mi) = scanFrequency(mf,mi)+1;
        
%         if size(scanHistory,1)>= config.maxIteration 
%             fprintf('\n!!!Maximum iteration number reached!');
%             break;
%         end
    end
    
%--------------------Sequential search----------------------
else
    
     % screen order among the rate constants
    config.sequence2 = [ 1  2;     
                         2  1];   
    config.sequence3 = [1 3 2 2 3 1;  % generate number state = 3 screen sequence for rate constants order
                        2 2 1 3 1 3];
    config.sequence4 = [3 2 3 1 4 2 2 3 1 4 1 4;                
                        2 3 1 3 2 4 1 4 2 3 4 1];
    config.sequence5 = [5 1 4 2 2 4 5 1 3 3 2 5 1 4 1 5 3 3 2 4 1 5 4 2;  
                        1 5 2 4 1 5 3 3 4 2 3 1 5 3 2 4 5 1 5 1 4 2 1 4]; % one can generage screen seqence for 6-state or above
    config.sequence6 = [1 6 2 5 2 4 1 6 3 5 3 4 1 6 3 4 2 5 4 3 6 1 5 2 6 1 4 3 2 5 5 2 6 1 6 1;
                        6 1 6 1 3 3 2 5 4 4 2 5 5 2 6 1 5 2 6 1 4 3 3 4 3 4 2 5 1 6 1 6 1 6 2 5];
    config.sequence7plus = [];
                    for i = 1 : numState
                        for j = 1 : numState
                            if i~=j;
                                a = [i;j];
                                config.sequence7plus = [config.sequence7plus,a];
                            end
                        end
                    end % if number of state >= 6, a default searching order is provided. One can customize this order. 

    switch numState
        case 2
           config.sequence = config.sequence2;
        case 3
           config.sequence = config.sequence3;
        case 4
           config.sequence = config.sequence4;
        case 5
           config.sequence = config.sequence5;
        case 6
           config.sequence = config.sequence6;
        otherwise
           config.sequence = config.sequence7plus;
    end

    
    bestGuess = config.rateTarget;
    save 'config'
 
    % minScore = findRateError(config, initialGuess);
    bestGuess = config.rateTarget;  
      

    rateHistory = [];
    wl1ScoreHistory = [];
    scanHistory = [];
    
    rateSimu = findRateSimu(config, bestGuess);
    minScore = findWL1(rateSimu, config.rateTarget);
    
    while minScore > config.WL1thresh
        % initialize the scales based on previous iterations
        
       for seq = 1:size(config.sequence,2)
            mi = config.sequence(1,seq);
            mf = config.sequence(2,seq);
              fprintf(['Transition from state ' num2str(mi),' to ', num2str(mf), ' is being scanned (col to row).\n']);
            rate = bestGuess;
            rateScore =[];
            rateScan = [];
            rateScanSimu = [];
            for j = 1:config.numoftry
                rate(mf,mi) = config.rateMulti(j)*bestGuess(mf,mi);
                for k = 1:config.repeatTime
                  fprintf('.');
                  rateSimu(:,:,k) = findRateSimu(config, rate);
                  wl1Score(k) = findWL1(rateSimu(:,:,k), config.rateTarget);
                end
                wl1Scoremean = mean(wl1Score);
                rateSimumean = mean(rateSimu,3);
                rateScan{j} = rate;
                rateScanSimu{j} = rateSimumean;
                rateScore(j) = wl1Scoremean;
            end
            [minScore,ind] = min(rateScore);
                fprintf(['\n score = ', num2str(minScore), '\n']);
            bestGuess = rateScan{ind};
                disp(num2str(bestGuess));
            bestGuessSimu = rateScanSimu{ind};
            rateHistory = cat(3, rateHistory, bestGuess);
            wl1ScoreHistory = [wl1ScoreHistory, minScore];
       end
       
       
       % show the result after each cycle
        fprintf('\nThe best guess is: \n');
        disp(bestGuess);
    end

end

%% Summary
fprintf('\n--------------------------------------------------\n');
fprintf('True rates =\n');
disp(num2str(config.rateTrue));

fprintf('\nExperimental rates =\n');
disp(num2str(config.rateTarget));

fprintf('\nGuess rates =\n');
disp(num2str(bestGuess));

fprintf('\nwL1_AT score = ');
wL1 = findWL1(bestGuess, config.rateTrue);
fprintf([num2str(wL1),'\n']);

display(['Note: wL1_SE target was set to ', num2str(config.WL1thresh)]);
if strcmp(config.RemoveOnePointData, 'true')
    fprintf('Note: one-point-data reassigned.\n');
else
    fprintf('Note: one-point-data not corrected.\n');
end

toc
%% save results

FileName = sprintf('Results_%s.mat',datestr(now,'yyyymmdd_HHMMSS'));
save(FileName); 
%}