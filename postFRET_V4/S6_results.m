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


fprintf('\n average rates = (cols: start states; rows: end states) ');
mr = mean(rateHistory(:, :, 4:end),3)

fprintf('\n standard deviation of the rates = ');
rs = std(rateHistory(:, :, 4:end), 0, 3)*sqrt(config.repeatTime)


[rateSimu, DataBin, DataBinLength] = MCsimulation(config, mr); % simulate a new trejectory using the final results.


figure; plot(data(:,1), data(:,5)); hold on; plot(DataBin(:,1), DataBin(:,4)+ 1); title('lower: norm. raw data; upper: simulated data');
jcPlotStyle

rhistDA = config.rhistDA;
rhistFRET = config.rhistFRET;

figure; 
hr = round((max(data(:,2))-min(data(:,2)))/rhistDA);
subplot(2, 3, 1), hist(data(:,2),hr); title('donor exp.');  jcPlotStyle
   axis([min(Donor(:,1)-3*Donor(:,3)) max(Donor(:,1)+3*Donor(:,3)) 0 inf]); 
hold on;

hr = round((max(data(:,3))-min(data(:,3)))/rhistDA);
subplot(2, 3, 2), hist(data(:,3),hr); jcPlotStyle
  title('acceptor exp.'); axis([min(Accpt(:,1)-3*Accpt(:,3)) max(Accpt(:,1)+3*Accpt(:,3)) 0 inf]);

hr = round((max(data(:,5))-min(data(:,5)))/rhistFRET); 
  subplot(2, 3, 3), hist(data(:,5), hr); jcPlotStyle
  title('FRET exp.'); axis([-0.1 1.1 0 inf]);

trajD = DataBin(:,2);
hr = round((max(trajD)-min(trajD))/rhistDA);
subplot(2, 3, 4), hist(trajD, hr); jcPlotStyle
   title('donor simu.'); axis([min(Donor(:,1)-3*Donor(:,3)) max(Donor(:,1)+3*Donor(:,3)) 0 inf]);

trajA = DataBin(:,3);
hr = round((max(trajA)-min(trajA))/rhistDA);
subplot(2, 3, 5), hist(trajA, hr); jcPlotStyle
   title('acceptor simu.'); axis([min(Accpt(:,1)-3*Accpt(:,3)) max(Accpt(:,1)+3*Accpt(:,3)) 0 inf]);

trajFRET = DataBin(:,4);
hr = round((max(trajFRET)-min(trajFRET))/rhistFRET); 
subplot(2, 3, 6), hist(trajFRET, hr); jcPlotStyle
   title('FRET simu.'); axis([-0.1 1.1 0 inf]);

