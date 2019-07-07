function [rateError] = findRateError(config, rate)

% initialize the molecules
 numState = config.numState;

 % rate = config.rate;
 % rate = rate.*Ratemulti(ratemulti);
 FRETStates = config.FRETStates;
 % transProb = config.transProb;
 
 transProb = zeros(numState);  % if a state is going to trnasfer to another, the probility of transfering to different states is lined up [k1/(sum kn); (k1+k2)/(sum kn),...]
            % (y,x): x initial state, y final state
for x = 1:numState % initialize transProb
    transProb(1,x) = rate(1,x)/sum(rate(1:numState,x));
    for y = 2:numState
        transProb(y, x) = transProb(y-1,x) + rate(y,x)/sum(rate(1:numState,x));
    end
end
for x = 1:numState
    transProb(x,x) = 0;
end   

stateLifetime = 1./sum(rate,1);
stateConcerntation = stateLifetime./sum(stateLifetime);


delta_t = config.delta_t; % time step (s). Simulate data with delta_t at 10x finer than bin time to represent the real.
bintime = config.bintime; % (s)
binsize = round(bintime/delta_t); % data points per bin
numMol = config.numMol; 
% numMol = numMols(numMolsnum); % changing data length 
BleachLifetime = config.BleachLifetime; % (s), the everage life time of the dye before bleaching. The lifetime is set to be exponential decay
BleachtimeThres = config.BleachtimeThres; % (s), if the lifetime is less than this value, remove it.
%totaltime = 1;  % seconds (s)

totaltime = zeros(1,numMol); % total simulation time for each molecule.
i = 1;
while i<=numMol
    % temptime = exprnd(BleachLifetime);
    temptime = BleachLifetime;
    if temptime >= BleachtimeThres
        totaltime(i) = temptime;
        i = i+1;
    end
end
cycle = ceil(totaltime./delta_t); % Mount Carlo cycle number; Step length of each molecule


%% main Mount Carlo function

molecules = struct('traject', []);
for m = 1: numMol % initialize each trajectory
%------- initial FRET state randomized with ratio based on the equilium concerntation of the states 
%    molecules(m).traject = ones(1,cycle(m)); %initialize all molecules in the whole time in state 1
    P = rand(1);
    for f = numState:-1:1
       if P <= sum(stateConcerntation(1:f)) 
          molecules(m).traject = ones(1,cycle(m))*f;
       end
    end
%------- end initial FRET states
    for i = 2:cycle(m)  % loop for simulation cycles/steps
        newmolecules = molecules(m).traject(i-1); % parent states of all molecules.
        istate = newmolecules;
        p1 = rand(1); % dwell time of a state judge
        if p1 <= 1-exp(-sum(rate(1:numState,istate))*delta_t)     % need transfer to another state
                % A = A0*(1-exp(-sum(kn-f)*t)
            p2 = rand(1);   % generate another random number to check with final state to transfer
                            % generate probability matrix.
            % consider pregenerate for every states.
            for f = 1:numState  % final state
                if f ~= istate && p2 <= transProb(f, istate)
                    newmolecules = f; % convert to final state
                    break
                end
            end
        end
        molecules(m).traject(i) = newmolecules; % record the states for the molecule.
    end

end

%plot the summary of all the trajectories
totalTrajectory = [];
for m = 1:numMol
    totalTrajectory = [totalTrajectory, molecules(m).traject]; %%%% total trajectar
end

eachTrajLength = cycle; %%%%%% traject length of each molecule

% concentration_t = 0:delta_t:totaltime-delta_t; 
% % plot concentration change over time (s)
% if CheckFigure
%     figure;  
%     for state = 1:numState
%         hold on;
%         subplot(1,numState,state); plot(0:delta_t:totaltime-delta_t,concentration(:,state));
%     end
% end

% find the real number of transitions
% split the vector into channeachTrajLength with each channel contain one state
% h=[];
% hmmden = totalTrajectory;
% for i = 1:numStates
%     h(i, :) = (hmmden == FRETstates(i));  %split the vector into numStates channel,with each channel contain one state
% end

%% initialize raw data matrix
totalTrajectory = [(1:size(totalTrajectory,2)).*delta_t; totalTrajectory];
%eachTrajLength = ones(1,numMol).*(totaltime/delta_t); 
eachTrajLength = cycle;
for i = 1:size(totalTrajectory,2)
        f = FRETStates(totalTrajectory(2, i)); %FRET value
        t = config.totalsignal;  % donor + acceptor total counts
        totalTrajectory(3,i) = t - f*t ;  % donor counts
        totalTrajectory(4,i) = f*t ; % acceptor counts
end
totalTrajectory(2,:) = totalTrajectory(4,:) ./ (totalTrajectory(4,:) + totalTrajectory(3,:));

Data = zeros(4,size(totalTrajectory,2));
Data(1,:) = totalTrajectory(1,:); %time
Data(2,:) = totalTrajectory(3,:); %donor
Data(3,:) = totalTrajectory(4,:); %acceptor
Data(4,:) = totalTrajectory(2,:); %FRET
%clm 1 time (s); clm 2 donor; clm 3 acceptor; clm 4 FRET;
Datalength = eachTrajLength;

if strcmp(config.CheckFigure, 'true')
    figure; plot(Data(1,:),Data(4,:));
end


%% 1. analyse the data: noise no;  bin no.
% [rateconstantsfit1, numtransitionfit1, ratesum1] =CountTransitions_difflength(Data, Datalength, FRETStates, config.RemoveOnePointData);

%% 2. analyse the data: noise no;  bin yes.

%% 4. analyse the data: noise yes; bin yes

%% add noise
noise = zeros(1,size(Data,2));
% noise = awgn(noise,ampsnr,'measured', 'linear')-noise; %create white noise
noise = awgn(noise,config.amplifiernoise);
Data(2,:) = Data(2,:)+ noise.*Data(2,:); %donor
Data(2,:) = round(abs(Data(2,:)));
noise = zeros(1,size(Data,2));
% noise = awgn(noise,ampsnr,'measured', 'linear')-noise; %create white noise
noise = awgn(noise,config.amplifiernoise);
Data(3,:) = Data(3,:) + noise.*Data(3,:); %acceptor
Data(3,:) = round(abs(Data(3,:)));
%disp('noise added to both doner channel and acceptor channel');

% add background noise
bknoise = ones(1,size(Data,2)).*config.backgroundnoise;
bknoise = poissrnd(bknoise)-config.backgroundnoise;
Data(2,:) = Data(2,:) + bknoise; % donor
bknoise = ones(1,size(Data,2)).*config.backgroundnoise;
bknoise = poissrnd(bknoise)-config.backgroundnoise;
Data(3,:) = Data(3,:) + bknoise; %acceptor
Data(4,:) = Data(3,:) ./ (Data(2,:) + Data(3,:));

%% bin data
% function [DataBin,DataBinlength]=binData(Data,Datalength)
DataBin = [];
DataRange = Datalength;
DataBinLength = zeros(1, numMol);
for i = 2:numMol
DataRange(i) = sum(Datalength(1:i));
end
for i = 1:numMol
    if i == 1
        DataMol = Data(:, 1:Datalength(i));
    else
        DataMol = Data(:, DataRange(i-1):DataRange(i));
    end
    DataBinLength(i) = floor(size(DataMol,2)/binsize); % new length of the data
    DataMolChop = DataMol(:, 1:DataBinLength(i)*binsize); % choped into multi of binsize
    
    DataBinMol = (1:DataBinLength(i)).*bintime; % first line being time
    
    bintemp = DataMolChop(2,:); %donor
    bintemp = reshape(bintemp, binsize, DataBinLength(i));
    DataBinMol(2,:) = sum(bintemp,1);  %donor
    bintemp = DataMolChop(3,:); %acceptor
    bintemp = reshape(bintemp, binsize, DataBinLength(i));
    DataBinMol(3,:) = sum(bintemp,1);  %acceptor
    DataBinMol(4,:) = DataBinMol(3,:)./(DataBinMol(2,:)+DataBinMol(3,:));
    %%% row 1 time, row 2 donor, row 3 acceptor, row 4 FRET
    DataBin = [DataBin,DataBinMol];
end
DataBin(1,:) = (1:size(DataBin,2)).*bintime; % connect the time for each molecule

if strcmp(config.CheckFigure, 'true')
    figure;
    plot(DataBin(1,:), DataBin(2,:)); title('after binning, donor blue, acceptor red');
    hold on
    plot(DataBin(1,:), DataBin(3,:), 'red');
    
    figure;
    s = 1; 
    for i = 1:numMol
        e = sum(DataBinLength(1:i));
        if rem(i,2)
            c = 'red';
        else
            c = 'blue';
        end
        plot(DataBin(1,s:e), DataBin(4,s:e), c);
        hold on
        s = e;
    end
    title('FRET after binning');
end
% data analysis
% vbFRET
[rateconstantsfit4, numtransitionfit4, ratesum4] = CountTransitions_difflength(DataBin, DataBinLength, FRETStates, config.RemoveOnePointData);
% figure; imagesc(numtransitionfit4); % map the number of transition(probability) from a state to another


%% find the rateError as the sum of errors. Percentage if smaller than 1, flip it.
%{
rateError = rateconstantsfit2(ratex,ratey)/config.rateTarget(ratex,ratey);
if rateError <1
    if rateError ~=0
    rateError = 1/rateError;
    end
end
rateError = rateError - 1;
%}
%
rateErrorRaw = rateconstantsfit4./config.rateTarget;
for i = 1:numState
    for j = 1:numState
        if i ~= j
            if rateErrorRaw(i,j) <1
                rateErrorRaw(i,j) = 1 - rateErrorRaw(i,j);
            else
                rateErrorRaw(i,j) = rateErrorRaw(i,j) - 1;
            end
        else 
            rateErrorRaw(i,j) = 0;
        end
    end
end

rateError = sum(sum(rateErrorRaw))/(numState*numState-numState);
%}