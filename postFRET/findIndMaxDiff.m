function [mf,mi] = findIndMaxDiff(rate, rateRef)

wl = abs(rate-rateRef);
[r,c] = find(wl==max(wl(:)));
mf = r;
mi = c;
end