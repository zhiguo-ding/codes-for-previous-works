function [ z ] = f_n( y, lamb,t,m,p,al,ep1, ximin,snr)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
gy1 = ((1+ep1)/ximin/(snr-ep1*y^al))^(-1/al);
 
z =  exp(-lamb*pi*y^2)*y^(2*(t-m-1)-2*p+1)/(2*m+2*p)*(y^(2*m+2*p) -gy1^(2*m+2*p) );
end

