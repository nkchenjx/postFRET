
%*************************************************************************************
%         This function pertains to the addition of AWGN with mean zero and
%         parameter 'sigma' to an input signal.
%
%         AUTHOR: Jixin Chen
% DATE  : 07/15/2015
%
% SYNOPSIS: y = awgn(x,sigma)
%         x   ---> input signal  (row)
%         sigma ---> standard deviation of the Gaussian noise
%         y   ---> y = x + AWGN   
%***********************************************************************************
function y = awgn(x,sigma)
    L = length(x);
    y = sigma*randn(1,L)+x;
return
    
    
    
