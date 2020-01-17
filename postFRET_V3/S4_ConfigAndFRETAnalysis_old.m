% 2016 by 
% Dr. Jixin Chen 
% Department of Chemistry and Biochemistry
% Ohio University,
% Athens, Ohio 45701 
% chenj@ohio.edu
% 
%    clear all;
%    close all;
%% manually determine states
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

% config.FRETUL = 1.1; % FRET upper limit
% config.FRETLL = -0.1; % FRET lower limit


%% analyze the raw data
config.RemoveOnePointData = 'true'; % if state A-B-C and B only last one data point, B will be merged to A if 'true'
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

% analyze raw data for the target rate
% data2 = data(:,[1,2,3,5]);
[rateconstantsfit, numtransitionfit, ratesum] =  CountTransitions_difflength(data(:,[1,2,3,5]),datalength, config.FRETStates,config.RemoveOnePointData);
config.rateTarget = rateconstantsfit;


fprintf('\n\n--------------------------------------------------\n');
fprintf('True rates =\n');
disp(num2str(config.rateTrue));

fprintf('\nExperimental rates =\n');
disp(num2str(config.rateTarget));

fprintf('\nwL1_AT score for raw data = ');
wL1 = findWL1(config.rateTarget, config.rateTrue);
fprintf([num2str(wL1),'\n']);
fprintf('\n--------------------------------------------------\n\n');


%--------------------over Writing the single simulation results to average
%results of five simulations:
% config.rateTarget = [0,11.0490484041813,6.99900675780571,4.62117144460293,3.12378208684813,1.77988941058650;18.5350496509884,0,14.4624683203641,9.51391094906703,5.67664047804632,4.14623983875533;14.1908731764572,19.2691468377841,0,17.5685738498993,11.7096255248044,7.99941682980906;9.16397389645158,12.7451014313252,18.4262656751761,0,19.8864730421171,15.3623301340084;4.91208151300297,6.05741911733557,9.27719742196829,14.9463379325173,0,19.4879063883545;1.76809060166991,2.64680105448343,4.00628653127624,6.89471292033787,11.0523423737904,0];
%--------------------





