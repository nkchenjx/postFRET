%% Data Analysis to get the transition rate constants
% Writen by Jixin Chen @ Rice University  jixin.chen@rice.edu  10/2013
%%%% rate constant k_n-f = N_n-f / t_n
% load 'eachTrajLength.mat' % the length of each molecule
% load '6Mureapath.dat' % the fitted data
% %load fitted states value
% hmmden = X6Mureapath(:,5);
% states = sort(states, 'descend');
function [rateconstantsfit, numtransitionfit, ratesum] =  CountTransitions_difflength(Data,Datalength,FRETStates,RemoveOnePointData)
CheckFigure = 0;

FRETTraj = Data(4,:);
bin_t = Data(1,:);
states = FRETStates;
numStates =  size(FRETStates,2);
tempstate = zeros(numStates,size(FRETTraj,2));
FRETTrajArray = tempstate;
bintime = bin_t(2) - bin_t(1);
for i = 1:size(FRETTraj,2)
    tempstate(:,i) = states';  
end
for i = 1:numStates
    FRETTrajArray(i,:) = FRETTraj;
end
FRETTrajArray = abs(FRETTrajArray - tempstate);
hmmden = zeros(1,size(FRETTraj,2)); 
for i = 1:size(FRETTraj,2)
    [x,ind] = min(FRETTrajArray(:,i));
    hmmden(i) = states(ind); 
end

clear tempstate FRETTrajArray;
% hmmden = round(FRETTraj);  %find the states by just round it to the
% nearest state.

if CheckFigure
    figure; plot(bin_t, FRETTraj); title('state identified');
    hold on; plot(bin_t, hmmden, 'red');
end
%hmmden = [totalTrajectary, totalTrajectary2];
% bintime = 0.01;

 ratesum = zeros(2, numStates);
 
  
%% Merge any state that last just one data point into its previous state

if strcmp(RemoveOnePointData, 'true')
    for j = 2:length(hmmden)-1
       if hmmden(j) ~= hmmden(j-1) & hmmden(j) ~= hmmden(j+1) & hmmden(j-1)~= hmmden(j+1)
           hmmden(j) = hmmden(j-1);
       end
    end

% merge any state that last two data points into its previous state
%{
for j = 2:length(hmmden)-2
   if hmmden(j) ~= hmmden(j-1) & hmmden(j) == hmmden(j+1) & hmmden(j) ~= hmmden(j+2) & hmmden(j-1)~hmmden(j+2)
       hmmden(j) = hmmden(j-1); hmmden(j+1) = hmmden(j-1);
   end
end
%}
end

%% split the vector into channeachTrajLength with each channel contain one state

h=[];
for i = 1:numStates
    h(i, :) = (hmmden == states(i));  %split the vector into numStates channel,with each channel contain one state
end

%% set the end of each molecule to 0 so the transition will not be counted
eachTrajLength = Datalength;
numMol = size(eachTrajLength,2);
eachTrajLength2 = eachTrajLength; % 
for i = 2 : numMol
    eachTrajLength2(i) = eachTrajLength2(i-1) + eachTrajLength(i); % eachTrajLength2 is accumulation of eachTrajLength that representing the end of each molecule
end
h(:,eachTrajLength2) = 0;

%% count the transition
dwellHist = struct('dh', [ ],'dhcumu',[], 'x', [], 'y', [], 'xy', []); %cumuHist
dwellhist2 = struct('dh', [ ], 'x', [], 'y', [], 'xy', []); % hist
rateconstants = [];
numdwell = [];
numtransition = [];
for n = 1:numStates
    for f = 1:numStates
        dwellhist = [];
        N=0;
        tn= sum(h(n,:));
        dwell = 0;
        for i = 1:size(h,2)-1
            if h(n,i) == 1
                dwell= dwell+1;
                if h(n,i+1) == 0
                     if h(f,i+1) == 1
                         if ~ismember(i-dwell, eachTrajLength2) %remove the starting first state of a molecule
                             N =N+1;
                             dwellhist = [dwellhist;dwell];
                         end
                    end
                    dwell = 0;
                end
            end
        end
        % N
        % t1
        % k = N/tn /bintime
      % disp(['n =', num2str(n), '  f =', num2str(f),' time=', num2str(sum(dwellhist))]);
        numtransition(f,n) = N;
        numdwell(f,n) = sum(sum(dwellhist));    %%%%% the last state is excluded from the analysis
        dwellHist(f,n).dh = dwellhist*bintime;
        dwellhist2(f,n).dh = dwellhist*bintime;
    end
end

for i = 1: size(dwellHist,1) %cumulative histgram of each transition
   for j = 1: size(dwellHist,2)
    if ~isempty(dwellHist(i,j).dh)
        dwellHist(i,j).dh = sort(dwellHist(i,j).dh, 'descend');
        dwellHist(i,j).dhcumu = (1:size(dwellHist(i,j).dh))';
        for k = 1:size(dwellHist(i,j).dh)-1
            if dwellHist(i,j).dh(k+1) ~= dwellHist(i,j).dh(k)
                dwellHist(i,j).x = [dwellHist(i,j).x; dwellHist(i,j).dh(k)];
                dwellHist(i,j).y = [dwellHist(i,j).y; dwellHist(i,j).dhcumu(k)];
            end
        end
        dwellHist(i,j).x = [dwellHist(i,j).x; dwellHist(i,j).dh(k+1)];
        dwellHist(i,j).y = [dwellHist(i,j).y; dwellHist(i,j).dhcumu(k+1)];
        dwellHist(i,j).xy = [dwellHist(i,j).x,dwellHist(i,j).y];
    end
    end
end

for i = 1: size(dwellhist2,1)  %histgram of each transition
    for j = 1: size(dwellhist2,2)
    if ~isempty(dwellhist2(i,j).dh)
         [dwellhist2(i,j).y,dwellhist2(i,j).x] = hist(dwellhist2(i,j).dh,10);
         dwellhist2(i,j).xy = [dwellhist2(i,j).x',dwellhist2(i,j).y'];       
    end
    end
end
 
 mol = [0,eachTrajLength2];
 molstates = zeros(numMol,numStates);
 for i = 1:numMol
     for j = 1:numStates
         temp = hmmden(mol(i)+1:mol(i+1));
         if size(temp(temp == states(j)),1)
            molstates(i,j) = 1;
         end
     end
 end
molstates(:,numStates+1)= sum(molstates(:,1:numStates),2); %count number of states covered by each molecule

%%%% fit the histogram to get the number of transitions.
lifetime = zeros(numStates);
numtransitionfit = zeros(numStates);
rateconstantsfit = zeros(numStates); % removed 2% tail data before fitting
% figure;
for i = 1: numStates %cumulative histgram of each transition
   for j = 1: numStates
    if ~isempty(dwellHist(i,j).xy)
        xy = dwellHist(i,j).xy;
        %%%remove the last 2% data before fitting
     %   xy = dwellHist(i,j).xy(dwellHist(i,j).xy(:,2)./max(dwellHist(i,j).xy(:,2))>0.02,:);
    
       %{
        if size(xy,1)>5
            p = polyfit(xy(:,1),log10(xy(:,2)),1); % log10 scale
            numtransitionfit(i,j) = 10^(p(2));
       %     numtransitionfit(i,j) = 10^(polyval(p,bintime));
            lifetime(i,j) = 1/(-p(1).*log(10));
    %         f = polyval(p,dwellHist(i,j).xy(:,1));
    %         hold on;  subplot(numStates,numStates, (i-1)*numStates+j);
    %         plot(dwellHist(i,j).xy(:,1),log10(dwellHist(i,j).xy(:,2))); % Log10 scale
    %         %axis([0 0.4 0 4]);
    %         hold on; plot(dwellHist(i,j).xy(:,1),f, 'red');
        else % if the cummulative distribution has too few elements to be fitted, take the average lifetime.
      %}
            numtransitionfit(i,j) = max(xy(:,2)); 
            lifetime(i,j) = sum(xy(:,1).*xy(:,2))/sum(xy(:,2));
  %          disp('cumulative distribution data points not enough for a fitting. Averaging instead!');
      %  end
    elseif i~=j
         numtransitionfit(i,j) = 0; 
         lifetime(i,j) = Inf;
    end
    
   end
end
     
%%%% calculate the rate constants.
numtransition(numStates+1, :) = sum(numdwell,1).*bintime;  % number of transitions between each state, last line is the total time of each state exclude` the final dewell of each molecule
    %%%Based on the number of transitions counted
rateconstants = zeros(numStates,numStates);
for i = 1:numStates
rateconstants(i,:) = numtransition(i,:)./numtransition(numStates+1,:);  %rate constants of each transition
end
stateSumm = [states;numtransition(numStates+1,:)]; stateSumm = stateSumm';  % probability of each state
    %%%Based on the number of transitions fitted from the dwell hist

for i = 1:numStates
%rateconstantsfit(i,:) = numtransitionfit(i,:)./numtransition(numStates+1,:);  %rate constants of each transition
  rateconstantsfit(i,:) = numtransitionfit(i,:)./sum(numtransitionfit.*lifetime,1); % rate constants of each transition
end

%%% calculate the sum dwell time of each state to all others
for i = 1:numStates
    dwellHistState(i) =struct('dh', [ ], 'dhcumu',[], 'x', [], 'y', [], 'xy', [], 'fitpoly',[],'lifetime', []);
    for f = 1:numStates
        dwellHistState(i).dh = [dwellHistState(i).dh; dwellHist(f,i).dh];
    end
end
for i = 1: numStates %cumulative histgram of each transition
    if ~isempty(dwellHistState(i).dh)
        dwellHistState(i).dh = sort(dwellHistState(i).dh, 'descend');
        dwellHistState(i).dhcumu = (1:size(dwellHistState(i).dh))';
        for k = 1:size(dwellHistState(i).dh)-1
            if dwellHistState(i).dh(k+1) ~= dwellHistState(i).dh(k)
                dwellHistState(i).x = [dwellHistState(i).x; dwellHistState(i).dh(k)];
                dwellHistState(i).y = [dwellHistState(i).y; dwellHistState(i).dhcumu(k)];
            end
        end
        dwellHistState(i).x = [dwellHistState(i).x; dwellHistState(i).dh(k+1)];
        dwellHistState(i).y = [dwellHistState(i).y; dwellHistState(i).dhcumu(k+1)];
        dwellHistState(i).xy = [dwellHistState(i).x,dwellHistState(i).y];
    end
end
%figure;
%{
for i = 1: numStates %cumulative histgram of each transition
   if ~isempty(dwellHistState(i).xy)
        p = polyfit(dwellHistState(i).xy(:,1),log10(dwellHistState(i).xy(:,2)),1);
        ratesum(:,i)= [-p(1); 10^(p(2))];
%         f = polyval(p,dwellHistState(i).xy(:,1));
%         hold on;  subplot(1,numStates, i);
%         plot(dwellHistState(i).xy(:,1),log10(dwellHistState(i).xy(:,2))); %  Log10 scale
%         axis([0 0.4 0 4]);
%         hold on; plot(dwellHistState(i).xy(:,1),f, 'red');
    end
end
%}

%% data to save

 %% make figures for histogram of each transition
%{
% if CheckFigure
     figure; 
     for i = 1: numStates %cumulative histgram of each transition
       for j = 1: numStates
        if ~isempty(dwellHist(i,j).xy)
            %%%remove the last 2% data before fitting
            xy = dwellHist(i,j).xy(dwellHist(i,j).xy(:,2)./max(dwellHist(i,j).xy(:,2))>0.02,:);
            p = polyfit(xy(:,1),log10(xy(:,2)),1);
            numtransitionfit(i,j) = 10^(p(2));
       %     numtransitionfit(i,j) = 10^(polyval(p,bintime));
            rateoffit(i,j) = -p(1);
            f = polyval(p,dwellHist(i,j).xy(:,1));
            hold on;  subplot(numStates,numStates, (i-1)*numStates+j); 
            bar(dwellHist(i,j).xy(:,1),log10(dwellHist(i,j).xy(:,2))); % log10 scale. Show all raw cummulated distribution including the last 2%.
            %axis([0 0.4 0 4]);
            hold on; plot(dwellHist(i,j).xy(:,1),f, 'red');
            title(['slop -', num2str(-p(1)), '; A ', num2str(10^(p(2)))]);
        end
        end
     end
     
    figure;
    for i = 1: numStates %cumulative histgram of each transition
       if ~isempty(dwellHistState(i).xy)
            % remove the last 2% data before fitting
            xy = dwellHistState(i).xy(dwellHistState(i).xy(:,2)./max(dwellHistState(i).xy(:,2))>0.02,:);
            p = polyfit(xy(:,1),log10(xy(:,2)),1);
            ratesum(:,i)= [-p(1); 10^(p(2))];
            f = polyval(p,dwellHistState(i).xy(:,1));
            hold on;  subplot(1,numStates, i); 
            bar(dwellHistState(i).xy(:,1),log10(dwellHistState(i).xy(:,2))); % log10 scale. Show all raw cummulated distribution including the last 2%.
            % axis([0 0.4 0 4]);
            hold on; plot(dwellHistState(i).xy(:,1),f, 'red');
            title(['slop -', num2str(-p(1)), '; A ', num2str(10^(p(2)))]);
        end
    end
% end
%}
%{
 figure; 
 for i = 1: numStates %cumulative histgram of each transition
   for j = 1: numStates
    if ~isempty(dwellHist(i,j).xy)
        %%%remove the last 2% data before fitting
        xy = dwellHist(i,j).xy(dwellHist(i,j).xy(:,2)./max(dwellHist(i,j).xy(:,2))>0.02,:);
        p = polyfit(xy(:,1),log10(xy(:,2)),1);
        numtransitionfit(i,j) = 10^(p(2));
   %     numtransitionfit(i,j) = 10^(polyval(p,bintime));
        rateoffit(i,j) = -p(1);
        f = polyval(p,dwellHist(i,j).xy(:,1));
        hold on;  subplot(numStates,numStates, (i-1)*numStates+j); 
        bar(dwellHist(i,j).xy(:,1),dwellHist(i,j).xy(:,2), 'red'); % log10 scale. Show all raw cummulated distribution including the last 2%.
        %axis([0 0.4 0 4]);
        hold on; plot(dwellHist(i,j).xy(:,1),10.^f, 'black');
        title(['lifetime ', num2str(1/(-2.303*p(1))), ' s; A ', num2str(10^(p(2)))]);
    end
    end
 end
%}

%% make 2D figures of transitions
%{
transitions = [];
radius = 50;
shape = zeros(2*radius+1);
map = zeros(1000);
sigmax = 20; sigmay = 20;
    for i = 1:2*radius+1
        for j = 1:2*radius+1
            shape(i,j)= 1*exp(-(i-radius)^2/(2*sigmax^2)-(j-radius)^2/(2*sigmay^2));
        end
    end

 for f = 1:size(FRETTraj,2)-1
     transitions(f,1) = FRETTraj(f).*1000;
     transitions(f,2) = FRETTraj(f+1).*1000;
     x0 = ceil(transitions(f,1));
     y0 = ceil(transitions(f,2));
     if x0>radius && x0<1000-radius && y0>radius && y0<1000-radius % && x0~=y0  
        map((x0-radius):(x0+radius),(y0-radius):(y0+radius)) =  map((x0-radius):(x0+radius),(y0-radius):(y0+radius)) +shape;
     end
 end   %the diaganol of the map represents the dwell histogram of the sates.

figure; imagesc(map)
%}
return