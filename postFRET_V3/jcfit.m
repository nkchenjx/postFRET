% All rights reserved by Jixin Chen
% First coding: 2016/04/05 by Jixin Chen @ Department of Chemistry and Biochemistry, Ohio University
% 20170110 Jixin Chen modified it to a function
% 20180609 Jixin Chen simplified it into single curve fitting

% Copyright (c) 2018 Jixin Chen @ Ohio University
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

%% an example for double exponential decay fitting of a single curve.
%{
    load('jcfitExampleData.mat');

    x1 = exampledata_x; % a row vector
    y1 = exampledata_y; % a row vector
      
    
    %% -------% the rest can be a double expnential function which only need feed
    % with x and y data.---------------------------
    
    %function [parafinal, rsq] = yournamedoubleExp(x1,y1)
    
    figure; plot(x1,y1); title('raw data');
    % this is a double exponential decay curve that have a dark part
    % find time 0 and remove the data before it starts, and normalize the curve 
    [m,ind] = max(y1);
    l = length(y1);
    x = x1-x1(ind);
    x = x(ind:l);
    y = y1(ind:l);
    y = y/m; 
    figure; plot(x,y); title('treated data'); 

    %------------END of loading data: x and y in row vectors--------
    
    % set fitting options
    option.maxiteration = 10;  % number of iteration fixed, the fitting will stop either this iteration or convergence reached first 
    option.accuracy = 0.0001;  % best searching accuracy, fraction of guessed value
    option.convgtest = 1e-10; % difference between two iterations on the square difference between fitting and data.

    % ----------------Attn: change below for different fitting equations-----------------
    % set the fitting equation to double exponential decay with a base line
    noisefitMdl = @(para, x) para(1)*exp(-(x/para(2))) + para(3)*exp(-(x/para(4))) + para(5);
    % equation grammar: modle name 'noisefitMdl' use as function y = noisefitMdl(para, x), the rest is the equation.
    % you can also generate a function to replace the above equation: 
    % function newy = noisefitMdl(para, x)
    
    % initial guess
    paraGuess = [0.5, 10, 0.5, 100, 0.1];  % A1, tau1,  A2, tau2, baseline
    % boundarys
    bounds = [0, 0.1, 0, 0.1, -0.1;   % lower boundary
              1, 2000, 1, 20000, 0.2]; % upper boundary
    %-------------------------------

    d1 = paraGuess-bounds(1,:);
    d2 = bounds(2,:)-paraGuess;
    if prod(d1)*prod(d2)<=0
        display('WARNING: initial guess out of boundary');
    end
    %--------------END of fitting option setting, equation, initial guess, and 

    %------------------and start fitting:------------------------
    [parafinal, yfit, chisq, rsq] = jcfit(noisefitMdl, x, y, paraGuess, bounds, option);
    fprintf(['\n rsq = ', num2str(rsq), '\n']);
    % parafinal is the fitted results; yfit is the fittered curve; 
    % use residual = y-yfit; to get the residual
    % rsq: root mean sqare value best will be close to 1
%}


%% main function with data error bar not fitted
function [parafinal, yfit, chisq, rsq] = jcfit(x, y, paraGuess, bounds, option)
    % load options
    if isempty(option.maxiteration)
     option.maxiteration = 50;  % number of iteration fixed
    end
    if isempty(option.accuracy)
     option.accuracy = 0.0001;    % best searching accuracy
    end
    if isempty(option.convgtest)
     option.convgtest = 1e-10; % difference between two iterations on the square difference between fitting and data.
    end
    maxiteration = option.maxiteration;
 %   accuracy = option.accuracy;
    convgtest = option.convgtest;
    minerrorlast = 0;
       
    para = paraGuess;
     
     tic
     for iteration = 1:maxiteration % fixed number of iteration.
         fprintf('.');
         if rem(iteration, 100) == 0 % progressing indicator
             fprintf('\n');
         end
         for i = 1:length(paraGuess) % scan each parameter
             %set the scanning scale withing the boundary.
            p = para(i);
            if abs(p) > option.accuracy
                accuracy = abs(p)*option.accuracy;
            else
                accuracy = option.accuracy;
            end        
            lb = bounds(1, i);
            ub = bounds(2, i);
            ll = p-lb;
            nl = floor(log2(ll/accuracy+1)):-0.5:1;
            ul = ub-p;
            nu = 1:0.5:floor(log2(ul/accuracy+1));
            ps = [lb, p-2.^nl*accuracy, p+2.^nu*accuracy, ub];
            error = [];
             % scan the parameter across the scale
             for j = 1: length(ps)
                para(i) = ps(j);
                error(j) = sum((noisefitMdl(para, x)-y).^2); %---the key equation: sum square of residual
             end
             % find the best
             [minerror, ind] = min(error); % find the least square.
             para(i) = ps(ind);
         end
         %test if converged
         if abs(minerror-minerrorlast)< convgtest; %convergence test positive
             fprintf('\n converged at iteration = ', num2str(iteration));
             break
         end
         minerrorlast = minerror;
     end
    fprintf('\n');
    toc
    %options=optimset('Maxiter',100000,'TolFun',1e10);
    % [para,r,J,COVB,mse] = nlinmultifit(x_cell, y_cell, noisefitMdl_cell, para0);

    % Calculate results

    %[ypred1, delta1] = nlpredci(noisefitMdl1, x, parareal, r, 'covar', mse);
    % figure

    parafinal = para;
    yfit = noisefitMdl(para, x);
    
    residual = y - yfit;
    a = (y - yfit).^2./yfit;
    a(isinf(a)) = 0;
    chisq = sum(a);

    sumy = sum(y.^2);
    N = length(x);

    sumr = sum(residual.^2);
    rsq = 1- sumr/sumy;
    
    % plot figures
    figure; plot(x,y,'linewidth',1.5); hold on; plot(x,yfit,'linewidth',1.5); plot(x, residual,'linewidth',1.5);
    title(['rsq = ', num2str(rsq)]);
    ax = gca;
    ax.LineWidth = 1.5;
    ax.Box = 'on';
    ax.TickLength = [0.02, 0.02];
    ax.FontName = 'Arial';
    ax.FontSize = 20;
    ax.FontWeight = 'Bold';

end

% by Jixin Chen