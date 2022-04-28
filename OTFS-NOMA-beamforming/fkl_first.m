function [ output_args ] = fkl_first( w,g,eps )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
output_args = -2*(1/sqrt(w'*g*g'*w)*g*g'*w - sqrt(eps)...
    /sqrt(w'*w)*w)/(sqrt(w'*g*g'*w) - sqrt(eps*w'*w))^3;

end

