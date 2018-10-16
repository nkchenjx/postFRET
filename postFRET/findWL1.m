function [wl1score] = findWL1(Rate, refRate)
[x,y] = size(Rate);

count = 0;
wl1score = 0;

for i = 1:x
    for j= 1:y
        if refRate(i,j)~=0
            count = count+1;
            wl1score = wl1score + abs(Rate(i,j)/refRate(i,j)-1);
        end
    end
end

wl1score = wl1score/count;
end