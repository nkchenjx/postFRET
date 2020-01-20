function [mf,mi] = findIndMaxErr(rate, rateRef)

wl = abs((rate-rateRef)./rateRef);
wl(isnan(wl)) = 0;
[r,c] = find(wl==max(wl(:)));
mf = r;
mi = c;
end