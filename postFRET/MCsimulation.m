%by Jixin Chen @Ohio University

function [ratefit, rateBin] = MCsimulation(config, rate)

numState = config.numState;
FRETStates = config.FRETStates;
delta_t = config.delta_t; % time step (s). Simulate data with delta_t at 10x finer than bin time to represent the real.
bintime = config.bintime; % (s)
binsize = round(bintime/delta_t); % data points per bin
numMol = config.numMol; 
BleachLifetime = config.BleachLifetime; % (s), the everage life time of the dye before bleaching. The lifetime is set to be exponential decay
BleachtimeThres = config.BleachtimeThres; % (s), if the lifetime is less than this value, remove it.


%% Initialize simulation conditoins
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


totalTrajectory = [];
for m = 1:numMol
    totalTrajectory = [totalTrajectory, molecules(m).traject]; %%%% total trajectar
end

eachTrajLength = cycle; %%%%%% traject length of each molecule




%% initialize raw data matrix
totalTrajectory = [(1:size(totalTrajectory,2)).*delta_t; totalTrajectory];
eachTrajLength = cycle;
for i = 1:size(totalTrajectory,2)
        f = FRETStates(totalTrajectory(2, i)); %FRET value
        t = config.totalsignal;  % donor + acceptor total counts
        totalTrajectory(3,i) = t - f*t ;  % donor counts
        totalTrajectory(4,i) = f*t ; % acceptor counts
end
totalTrajectory(2,:) = totalTrajectory(4,:) ./ (totalTrajectory(4,:) + totalTrajectory(3,:));

Data = zeros(4,size(totalTrajectory,2));
Data(1,:) = totalTrajectory(1,:); %time, unit s
Data(2,:) = totalTrajectory(3,:); %donor,  unit count
Data(3,:) = totalTrajectory(4,:); %acceptor,  unit count
Data(4,:) = totalTrajectory(2,:); %FRET, A.U.
%row 1 time (s); row 2 donor; row 3 acceptor; row 4 FRET;
Datalength = eachTrajLength;

if strcmp(config.CheckFigure, 'true')
    figure; plot(Data(1,:),Data(4,:));
end

%% Raw data backup
trjData = Data;


%% 1. analyse the non-binned data
% skip vbFRET because there is no noise and bin
[ratefit, numtransitionfit1, ratesum1] =CountTransitions_difflength(trjData, Datalength, FRETStates, config.RemoveOnePointData);


%% 2. analyse the binned data
%% bin data

DataBin = [];
DataRange = Datalength;
DataBinLength = zeros(1, numMol);
for i = 2:numMol
DataRange(i) = sum(Datalength(1:i));
end
for i = 1:numMol % bin each molecule
    if i == 1
        DataMol = Data(:, 1:Datalength(i));
    else
        DataMol = Data(:, DataRange(i-1)+1:DataRange(i));
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
    figure; plot(DataBin(1,:), DataBin(4,:), 'red'); title('FRET after binning');
end

trjDataBin = DataBin;

% analysis
[rateBin, numtransitionfit2, ratesum2] =CountTransitions_difflength(trjDataBin, DataBinLength, FRETStates, config.RemoveOnePointData);

%% 3. Save the data
if strcmp(config.CheckSaveTraj, 'true')
    FileName=sprintf('exptraj_%s.mat',datestr(now,'yyyymmdd_HHMMSS'));
    save(FileName, 'trjDataBin', 'DataBinLength', 'rateBin', 'trjData', 'Datalength', 'ratefit', 'config');
end
