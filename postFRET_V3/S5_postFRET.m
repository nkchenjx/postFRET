
% coded by Jixin Chen @Ohio University 2016
% modified 2019 for KinSoftChallenge

% postFRET Monte Carlo simulating to guess a trajectory that matches the experimental trajectory the best

%% initialize config
%config.rateTrue = config.rateTrue;
config.WAMSearching = 'false'; % sequential searching the rate constant or search the state that gives the maximum wL1_SE score (WAM).
config.BleachLifetime = mean(config.bleachtime);
config.BleachtimeThres = min(config.bleachtime); % (s), if the lifetime is less than this value, remove it.

config.delta_t = config.datatimestep; % time step (s). Simulate data with delta_t at 10x finer than bin time to represent the real.
config.bintime = config.datatimestep; % (s)
config.binsize = round(config.bintime/config.delta_t); % data points per bin
config.numMol = length(datalength); % default number of total molecules set to 100
   display(['simulation length of each run ~',num2str(config.numMol*config.BleachLifetime), ' s long data points...']);

   
% snr = 2; % signal to noise ratio before bin. give std 1/sqrt(snr).
% config.backgroundnoise = 0.7/sqrt(config.delta_t/0.001); % Poisson distribution, roughly average 0.7-0.8 counts typically for 1 ms data
% config.amplifiernoise = 0.37/sqrt(config.delta_t/0.001); % std(signal-mean)/mean, Gaussian, typical value 0.37+-0.03 for 1 ms data.

% config.backgroundnoise = 0; % Poisson distribution, roughly average 0.7-0.8 counts typically for 1 ms data
% config.noise = config.rawnoise*sqrt(config.binsize);
% config.amplifiernoise = config.sumnoise/config.sumIaa*sqrt(config.binsize); % std(signal-mean)/mean, Gaussian, typical value 0.37+-0.03 for 1 ms data.

  x = [-4, -2, -0.5, -0.2, -0.1, -0.045, -0.02, -0.005, 0, 0.006, 0.03, 0.05 0.11, 0.21, 0.51, 2.1, 4.1]; % for multiple factor for each rate constant  
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

%------------------------ maximum iteration number:
config.maxIteration = 60; %set the maximum iterations in a searching.
config.repeatTime = 3; % for each run repeat times to get the average
config.WL1thresh = 0.001; % the convergence wL1 score (WL1thresh*100)% for a "whack-a-mole" search to stop.  
                         % This value adjusts final accuracy. 
                         % This value has different limits for different conditions. 
                         % Ctrl+C to break during the simulation if it takes too long.                      
  display(['Set the wL1_SE target = ', num2str(config.WL1thresh)]);

% end of initializing parameters

%% Start simulation
fprintf('\n------------------------------------------------\n');
fprintf('\npostFRET Monte Carlo simulating to guess a trajectory that matches the experimental trajectory the best.\n');
%tic


numState = config.numState;
config.computTime = [];
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
          rateSimu(:,:,k) = MCsimulation(config, bestGuess);
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
        tic
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
            parfor k = 1:config.repeatTime
              fprintf('.');
              rateSimu(:,:,k) = MCsimulation(config, rate);
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
        
        if length(wl1ScoreHistory)>= config.maxIteration 
            fprintf('\n!!!Maximum iteration number reached!');
            break;
        end
        t = toc;
        config.computTime = [config.computTime; t];
        fprintf(['\n----------------- End of ', num2str(length(wl1ScoreHistory)),'/', num2str(config.maxIteration), ' iteration -----------------this iteration costs time = ', num2str(t), ' s\n']);
        
    end
    
%% --------------------Sequential search----------------------
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

    
  % save 'config'
 
    % minScore = findRateError(config, initialGuess);
    bestGuess = config.rateTarget;  
      

    rateHistory = [];
    wl1ScoreHistory = [];
    scanHistory = [];
    
    rateSimu = MCsimulation(config, bestGuess);
    minScore = findWL1(rateSimu, config.rateTarget);
    
    %-----------------------------------------------------------
    % searching algorithm: JCFit (https://github.com/nkchenjx/JCFit) :
    while minScore > config.WL1thresh 
        % initialize the scales based on previous iterations
       tic
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
                parfor k = 1:config.repeatTime
                  rateSimu(:,:,k) = MCsimulation(config, rate);
                  wl1Score(k) = findWL1(rateSimu(:,:,k), config.rateTarget);
                end
                wl1Scoremean = mean(wl1Score);
                rateSimumean = mean(rateSimu,3);
                rateScan{j} = rate;
                rateScanSimu{j} = rateSimumean;
                rateScore(j) = wl1Scoremean;
                fprintf('.');
            end
            [minScore,ind] = min(rateScore);
                fprintf(['\n score = ', num2str(minScore), '\n']);
            bestGuess = rateScan{ind};
                disp(num2str(bestGuess));
            bestGuessSimu = rateScanSimu{ind};

       end
       
       rateHistory = cat(3, rateHistory, bestGuess);
       wl1ScoreHistory = [wl1ScoreHistory, minScore];       
       % show the result after each cycle
        fprintf('\nThe best guess is: \n');
        disp(bestGuess);
        if length(wl1ScoreHistory)>= config.maxIteration 
            fprintf('\n!!!Maximum iteration number reached!');
            break;
        end
        t = toc
        config.computTime = [config.computTime; t];
        fprintf(['\n----------------- End of ', num2str(length(wl1ScoreHistory)),'/', num2str(config.maxIteration), ' iteration ------------ this iteration costs time = ', num2str(t), ' s\n']);
    end
end

if minScore <= config.WL1thresh 
    fprintf('\n!!!WL1 threshold reached (converged)! \n');
end

config.bestGuess = bestGuess;

%% Summary
fprintf('\n\n--------------------------------------------------\n');
fprintf('True rates = (cols: start states; rows: end states)\n');
disp(num2str(config.rateTrue));

fprintf('\nExperimental rates =(cols: start states; rows: end states)\n');
disp(num2str(config.rateTarget));

fprintf('\nGuess rates = (cols: start states; rows: end states)\n');
disp(num2str(bestGuess));

fprintf('\nwL1_AT score = ');
wL1 = findWL1(bestGuess, config.rateTrue);
fprintf([num2str(wL1),'\n']);

fprintf(['Note: wL1_SE target was set to ', num2str(config.WL1thresh)]);
if(wL1 < config.WL1thresh)
    fprintf(': achieved. \n');
else
    fprintf(': has not achieved. \n');
end

if strcmp(config.RemoveOnePointData, 'true')
    fprintf('Note: one-point-data reassigned.\n');
else
    fprintf('Note: one-point-data not corrected.\n');
end

figure; plot(wl1ScoreHistory); title('score history. See rate history in rateHistory');

fprintf('\n\n--------------------------------------------------\n');


fprintf(['\n total computational time = ', num2str(sum(config.computTime)), ' s\n']);
fprintf(['\n each iteration computational time = ', num2str(mean(config.computTime)), ' s\n']);

[rateSimu, DataBin, DataBinLength] = MCsimulation(config, bestGuess); % simulate and save one trajectory using the last bestGuess


%% display final results
figure; plot(data(:,1), data(:,5)); hold on; plot(DataBin(:,1), DataBin(:,4)+ 1); title('lower: norm. raw data; upper: simulated data');
[Nr,Xr] = hist(data(:,5),1000);
[Ns,Xs] = hist(DataBin(:,4),1000);
figure;  plot(Xr, Nr/max(Nr)); hold on; 
plot(Xs, ones(length(Xs),1)); plot(Xs, Ns/max(Ns)+1); title('lower raw data, upper simulated data FRET histogram');


fprintf('\n average rates = ');
mr = mean(rateHistory,3)

fprintf('\n standard deviation of the rates = ');
rs = std(rateHistory, 0, 3)*sqrt(config.repeatTime)


%% save results
cd(config.resultpath);
FileName = sprintf('Results_2tates_%s.mat',datestr(now,'yyyymmdd_HHMMSS'));
save(FileName); 

%}