% use total signal to find the noise level of the donor and acceptor


%% normalize the total photocounts to global mean
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

%% calculate noise
sumnoise = std(data(:,4)); % two channels sum noise
config.sumIaa = rawmean;
% config.rawnoise = sumnoise*sqrt(2); % noise of each channel
config.rawnoise = sumnoise/sqrt(2);