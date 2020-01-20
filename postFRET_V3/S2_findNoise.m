% use total signal to find the noise level of the donor and acceptor
%% use an intensity-dependent Gaussian model to fit the noise of donor and acceptor channels
% then save the noise model in config using the same timestep.

%% normalize the total photocounts to global mean to even experimental variations such as focus drifting.
data = rawdata;
rawmean = mean(data(:,4));

s = 1; e = datalength(1);
tracemean = mean(data(s:e, 4));
data(s:e, 2:4) = data(s:e, 2:4)*rawmean/tracemean;

for i = 2:length(datalength)
s = datalength(i-1)+1; e = datalength(i);
tracemean = mean(data(s:e, 4));
data(s:e, 2:4) = data(s:e, 2:4)*rawmean/tracemean;
end

figure; plot(data(:,2)); hold on;
plot(data(:,3)); plot(data(:,4)); title('normalized A D total');

sumIaa = mean(data(:,4)); % nomalized A+D total;

%% fit the gaussian peaks

[ds, dx] = hist(data(:,2), 500);
[as, ax] = hist(data(:,3), 500);

figure; plot(dx, ds); hold on; plot(ax, as);

%----------------- donor channel
figure;plot(dx, ds); title('donor channel: click the Gaussian peaks from left to right, and hit keybord enter to end');
pD  = ginput; %Gaussian peak positions, [x1, y1; x2, y2;...]
numP = size(pD,1);
% fit donor channel with gaussian
    x = dx; y = ds;
%    X = repmat(x, numP);
%    Y = repmat(y, numP);
    % set fitting options
    option.maxiteration = 10;  % number of iteration fixed, the fitting will stop either this iteration or convergence reached first 
    option.accuracy = 0.0001;  % best searching accuracy, fraction of guessed value
    option.convgtest = 1e-10; % difference between two iterations on the square difference between fitting and data.

    % ----------------Attn: change below for different fitting equations-----------------
    % set the fitting equation to double exponential decay with a base line
   % noisefitMdl = @(para, X) para(1) + para(2)*x;
    % equation grammar: modle name 'mdl' use as function y = mdl(para, x), the rest is the equation.
    % you can also generate a function to replace the above equation: 
    % function newy = mdl(para, x)
    
    % initial guess
    pg = pD;
    pg(:,3) = (pD(end, 1) - pD(1,1))/numP/2;
    pg = pg';
    paraGuess = pg(:)';  
    %------para structure------------------------
    % miu1, A1, sigma1, miu2, A2, sigma2, ...
    %-------
    
    
    % boundarys
    bounds = [paraGuess/10;   % lower boundary
              paraGuess*10]; % upper boundary
    %-------------------------------

    d1 = paraGuess-bounds(1,:);
    d2 = bounds(2,:)-paraGuess;
    if prod(d1)*prod(d2)<=0
        display('WARNING: initial guess out of boundary');
    end
    %--------------END of fitting option setting, equation, initial guess, and 

    %------------------and start fitting:------------------------
    [parafinalD, yfitD, chisqD, rsqD] = jcfit(x, y, paraGuess, bounds, option);
    fprintf(['\n rsq = ', num2str(rsqD), '\n']);
    % parafinal is the fitted results; yfit is the fittered curve; 
    % use residual = y-yfit; to get the residual
    % rsq: root mean sqare value best will be close to 1

%---------------- acceptor channel ------------  
figure;plot(ax, as); title('acceptor channel: click the Gaussian peaks from left to right, and hit keybord enter to end');
pA  = ginput; %Gaussian peak positions, [x1, y1; x2, y2;...]
numP = size(pA,1);
% fit donor channel with gaussian
    x = ax; y = as;
%    X = repmat(x, numP);
%    Y = repmat(y, numP);
    % set fitting options
    option.maxiteration = 10;  % number of iteration fixed, the fitting will stop either this iteration or convergence reached first 
    option.accuracy = 0.0001;  % best searching accuracy, fraction of guessed value
    option.convgtest = 1e-10; % difference between two iterations on the square difference between fitting and data.

    % ----------------Attn: change below for different fitting equations-----------------
    % set the fitting equation to double exponential decay with a base line
   % noisefitMdl = @(para, X) para(1) + para(2)*x;
    % equation grammar: modle name 'mdl' use as function y = mdl(para, x), the rest is the equation.
    % you can also generate a function to replace the above equation: 
    % function newy = mdl(para, x)
    
    % initial guess
    pg = pA;
    pg(:,3) = (pA(end, 1) - pA(1,1))/numP/2;
    pg = pg';
    paraGuess = pg(:)';
    
    % boundarys
    bounds = [paraGuess/10;   % lower boundary
              paraGuess*10]; % upper boundary
    %-------------------------------

    d1 = paraGuess-bounds(1,:);
    d2 = bounds(2,:)-paraGuess;
    if prod(d1)*prod(d2)<=0
        display('WARNING: initial guess out of boundary');
    end
    %--------------END of fitting option setting, equation, initial guess, and 

    %------------------and start fitting:------------------------
    [parafinalA, yfitA, chisqA, rsqA] = jcfit(x, y, paraGuess, bounds, option);
    fprintf(['\n rsq = ', num2str(rsqA), '\n']);
    
    
    
% %% calculate noise
% sumnoise = std(data(:,4)); % two channels sum noise
% config.sumIaa = rawmean;
% % config.rawnoise = sumnoise*sqrt(2); % noise of each channel
% config.rawnoise = sumnoise/sqrt(2);
% 

%% Determine the FRET states
    %------para structure------------------------
    % miu1, A1, sigma1, miu2, A2, sigma2, ...
    %-------
 numstate = numP;
 Donor = reshape(parafinalD, 3, numstate);
 Donor = Donor';
 Accpt = reshape(parafinalA, 3, numstate);
 Accpt = Accpt';
 Dc = Donor(:, 1);
 Ac = Accpt(:, 1);
 
 FRETStates = Ac./(Ac + flip(Dc));
 FRETStates = FRETStates';
 
 figure; hist(data(:,5), 100); title(['FRET states determined: ', num2str(FRETStates), '. Check code if wrong']);
 
%% determine the noise equation
% Use the count-dependent Gaussian sigma model for the noise now. Change
% when necessary.
[noisePD, SD] = polyfit(Donor(:,1), Donor(:,3), 1);
[noisePA, SA] = polyfit(Accpt(:,1), Accpt(:,3), 1);

%----------- manually boost the noise and manually determine the states-
% noisePD(1) = 2*noisePD(1);
% noisePA(1) = 1.5*noisePA(1);
   
%-----------

noiseMdl = @(para, x) para(1)*x + para(2); % linear model,

Ac = FRETStates*sumIaa;
Dc = (1-FRETStates)*sumIaa;

% test run the noise model and compair to the data.

L = length(data(:,1));
f = Accpt(:,2).*Accpt(:,3);
li = round(L*f./sum(f));

% f = Donor(:,2).*Donor(:,3);
% li = round(L*f./sum(f));

trajD = [];
for i = 1:numstate
    % trajD = [trajD; ones(li(i), 1)*Donor(i,1)];
    trajD = [trajD; ones(li(i), 1)*Dc(i)];
end

% f = Accpt(:,2).*Accpt(:,3);
% li = round(L*f./sum(f));
trajA = [];
%for i = numstate:-1:1;
    % trajA = [trajA; ones(li(i), 1)*Accpt(i,1)];
 for i = 1:numstate
    trajA = [trajA; ones(li(i), 1)*Ac(i)];
end

simutimestep = datatimestep;
trajD = noiseAdding(trajD, noiseMdl, noisePD, datatimestep, simutimestep);
trajA = noiseAdding(trajA, noiseMdl, noisePA, datatimestep, simutimestep);
trajFRET = trajA./(trajA+trajD);

figure; subplot(2, 2, 1), hist(data(:,2),500); title('donor exp'); axis([min(Donor(:,1)-3*Donor(:,3)) max(Donor(:,1)+3*Donor(:,3)) 0 inf]); hold on;
subplot(2, 2, 2), hist(data(:,3),500); title('acceptor exp'); axis([min(Accpt(:,1)-3*Accpt(:,3)) max(Accpt(:,1)+3*Accpt(:,3)) 0 inf]);
subplot(2, 2, 3), hist(trajD, 500); title('donor simulated noise'); axis([min(Donor(:,1)-3*Donor(:,3)) max(Donor(:,1)+3*Donor(:,3)) 0 inf]);
subplot(2, 2, 4), hist(trajA, 500); title('acceptor simulated noise'); axis([min(Accpt(:,1)-3*Accpt(:,3)) max(Accpt(:,1)+3*Accpt(:,3)) 0 inf]);


figure; hr = round((max(data(:,5))-min(data(:,5)))/0.01); subplot(2, 1, 1), hist(data(:,5), hr); title('Exprimental FRET'); axis([-0.1 1.1 0 inf]); hold on;
hr = round((max(trajFRET)-min(trajFRET))/0.01); subplot(2, 1, 2), hist(trajFRET, hr); title('Estimated simulated FRET'); axis([-0.1 1.1 0 inf]);

%---------manually determines the states! overwrite the above -----
display(num2str(Ac + flip(Dc)));
display('if the sum of donor and accepter channel are not the same, and the distributions are different, manually determine the FRET state below');
 %-----------force states:
 numstate = 2;
 FRETStates = [0.2, 0.8]; 
% % the above states were determined using STaSI
%-----------

%% summerise results
config.numState = numstate;
config.FRETStates = FRETStates;
config.noiseMdl = noiseMdl;
config.noisePD = noisePD;
config.noisePA = noisePA;
config.sumIaa = sumIaa;

