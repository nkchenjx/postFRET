%  1.  add Gaussian noise with std propotional to the true signal at ratio
%  of amplifiernoise;
%  2.  add background noise with Poisson distribution at average value
%  backgrandnoise. Then shift down the whole vector the average
%  backgroundnoise.


function [APDData] = addAPDNoise(trueData,amplifiernoise, backgroundnoise)
% amplifier noise
noise = zeros(1,size(trueData,2));
noise = awgn(noise,amplifiernoise);
APDData = trueData+ noise.*trueData; 
APDData(find(APDData<0)) = 0; %amplifier noise added

%background noise
bknoise = ones(1,size(trueData,2)).*backgroundnoise;
bknoise = poissrnd(bknoise);
APDData = APDData + bknoise; % background noise added

APDData = round(APDData)-backgroundnoise;
end