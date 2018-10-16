function [minX, stdX] = findRateMin(valX, valY, valE)

[minY, indminY] = min(valY);

%{
Eqn1 = 'a./x+b*x+c+d*exp(-x*e)';
startPoints = [valY(1)*valX(1)/2 max(valY)/max(valX) minY, max(valY)/max(valX), 1/valX(indminY)*5];
f1 = fit(valX', valY', Eqn1, 'start', startPoints);
Eqn2 = @(x)f1.a./x+f1.b*x+f1.c+f1.d*exp(-x*f1.e);
findx = min(valX):(max(valX)-min(valX))/10000:max(valX);
findy = Eqn2(findx);
[val,ind] = min(findy);
minX = findx(ind);

% find the index of the adjecent valX values
count = 0;
for i = 1:length(valX)
    if minX>valX(i)
        count = count+1;
    end
end
if count < 1
    count =1;
end
if count > length(valX)-1
    count = length(valX)-1;
end
% the minX is between count to count+1 element of valX

% error bar of the valY at this point set to the average of the two
s = (valE(count)+valE(count+1))/2;
eY = f1.a./minX+f1.b*minX+f1.c+f1.d*exp(-minX*f1.e)+s;

% find the intersections of the curve with the line y = y(minX)+s
x = min(valX): (minX-min(valX))/1000 :minX; % screen x at 1000 steps on the left
y = f1.a./x+f1.b*x+f1.c+f1.d*exp(-x*f1.e);
d = abs(y-eY);
[v,i] = min(d);
lx = x(i);
% figure; plot(x, y); hold on; plot(x, eY*ones(1,1001));

x = minX: (max(valX)-minX)/1000 :max(valX); % screen x at 1000 steps on the right
y = f1.a./x+f1.b*x+f1.c+f1.d*exp(-x*f1.e);
d = abs(y-eY);
[v,i] = min(d);
rx = x(i);

% find the error bar of the rate constants
stdX = minX-lx;  % set to be the left distance
if isempty(stdX)
 stdX = 0;
end
%}
%% cancel all the previous calculation, :( just use the minimum.
minX = valX(indminY);
stdX = 0;
% plot the fitting results
% figure; errorbar(valX, valY, valE);
% figure; plot(f1); hold on; errorbar(valX, valY, valE);