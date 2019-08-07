%% Summary
fprintf('\n\n--------------------------------------------------\n');
fprintf('FRET state selected = \n');
disp(config.FRETStates);

fprintf('True rates =\n');
disp(num2str(config.rateTrue));  %% rate = rate' to convert to kinsoft format.

fprintf('\nExperimental rates =\n');
disp(num2str(config.rateTarget));

fprintf('\nGuess rates =\n');
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


fprintf('\n average rates = ');
mr = mean(rateHistory,3)

fprintf('\n standard deviation of the rates = ');
rs = std(rateHistory, 0, 3)*sqrt(config.repeatTime)


[rateSimu, DataBin, DataBinLength] = MCsimulation(config, mr); % simulate a new trejectory using the final results.


figure; plot(data(:,1), data(:,5)); hold on; plot(DataBin(:,1), DataBin(:,4)+ 1); title('lower: norm. raw data; upper: simulated data');

[Nr,Xr] = hist(data(:,5),100);
[Ns,Xs] = hist(DataBin(:,4),100);
figure;  plot(Xr, Nr/max(Nr)); hold on; 
plot(Xs, ones(length(Xs),1)); plot(Xs, Ns/max(Ns)+1); title('lower raw data, upper simulated data FRET histogram');


