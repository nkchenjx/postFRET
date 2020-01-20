% fit the histgram of the photocounts with Gaussian functions.

function y = noisefitMdl(para, x)
    i = 1;
    y = para(i*3-1)*exp(-(x-para(i*3-2)).^2/2/para(i*3)^2);
    for i = 2:length(para)/3
        yt = para(i*3-1)*exp(-(x-para(i*3-2)).^2/2/para(i*3)^2);
        y = y+yt;
    end
end