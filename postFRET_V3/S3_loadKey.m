% load the key

%% Window to select file%%%%
% shown in simulation.cfg

% the true rate is for the training data purpose.
% give a random guess for unknown experimental results:
rate = [0   1; 
        1  0 ];

% rate = [0   0.32   0.2; 
%         0.11  0   0.39;
%         0.06   0.4   0];

% rate = [0   0.32   0.2   0.2; 
%         0.11  0   0.39   0.2;
%         0.06   0.4   0   0.2;
%         0.2    0.2  0.2   0];


% In postFRET the rate matrix is defined
%  [0,       2 to 1    3 to 1 ...
%  1 to 2,    0        3 to 2 ...
%  ...]


config.rateTrue = rate';  % flip the clm and row in the key because it was difined clm to row in postFRET. and change diagonol to 0.