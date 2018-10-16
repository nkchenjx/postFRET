function [mi] = findIndMaxDiffState(rate, rateRef)

wl = abs((rate-rateRef)./rateRef);
wl(isnan(wl))=0;
clsum = sum(wl, 1);
mi = find(clsum==max(clsum(:)));
mi = mi(1);
end