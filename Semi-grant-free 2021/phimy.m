function [ y ] = phimy( p,mu,alpha0,alpha1,alpha2,alphas,Ps )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
 
    y = exp(-p*alphas)*gmu(mu,alpha1,alpha2) +...
        exp(p/Ps)*gmu(mu+p/Ps/alpha0,alpha0,alpha1);

end

