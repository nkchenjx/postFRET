
% coded by Jixin Chen @Ohio University 2016
% modified 2019 for KinSoftChallenge

%% postFRET Monte Carlo simulating to guess a trajectory that matches the experimental trajectory the best
fprintf('\n------------------------------------------------\n');
fprintf('\npostFRET Monte Carlo simulating to guess a trajectory that matches the experimental trajectory the best.\n');
tic


config.maxIteration = 5; %set the maximum iterations in a searching. Typical value 5-50. the computational time is propotional to this value.
config.repeatTime = 3; % for each run repeat times to get the average. Typical value 3-5.


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
          rateSimu(:,:,k) = findRateSimuNoise(config, bestGuess);
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
            parfor k = 1:config.repeatTime
              fprintf('.');
              rateSimu(:,:,k) = findRateSimuNoise(config, rate);
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
        fprintf(['\n----------------- End of ', num2str(length(wl1ScoreHistory)), ' iteration.-----------------\n']);
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

    
    bestGuess = config.rateTarget;
  % save 'config'
 
    % minScore = findRateError(config, initialGuess);
    bestGuess = config.rateTarget;  
      

    rateHistory = [];
    wl1ScoreHistory = [];
    scanHistory = [];
    
    rateSimu = findRateSimuNoise(config, bestGuess);
    minScore = findWL1(rateSimu, config.rateTarget);
    
    %-----------------------------------------------------------
    % searching algorithm: JCFit (https://github.com/nkchenjx/JCFit) :
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
                parfor k = 1:config.repeatTime
                  fprintf('.');
                  rateSimu(:,:,k) = findRateSimuNoise(config, rate);
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
        fprintf(['\n----------------- End of ', num2str(length(wl1ScoreHistory)), ' iteration.-----------------\n']);
    end
end


config.bestGuess = bestGuess;

%% Summary
fprintf('\n\n--------------------------------------------------\n');
fprintf('True rates =\n');
disp(num2str(config.rateTrue));

fprintf('\nExperimental rates =\n');
disp(num2str(config.rateTarget));

fprintf('\nGuess rates of the last iteration =\n');
bestGuess = rateHistory(:,:, end);
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

toc

[rateSimu, DataBin, DataBinLength] = MCsimulation(config, bestGuess);

%% save results
%cd(config.resultpath);
cd('.\exampleResults\');
FileName = sprintf('Results_%s.mat',datestr(now,'yyyymmdd_HHMMSS'));
save(FileName); 


figure; plot(data(:,1), data(:,5)); hold on; plot(DataBin(:,1), DataBin(:,4)+ 1); title('lower: norm. raw data; upper: one of the simulated data');
[Nr,Xr] = hist(data(:,5),100);
[Ns,Xs] = hist(DataBin(:,4),100);
figure;  plot(Xr, Nr/max(Nr)); hold on; 
plot(Xs, ones(length(Xs),1)); plot(Xs, Ns/max(Ns)+1); title('lower raw data, upper one of the simulated data FRET histogram');


%}