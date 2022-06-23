clear all%function [p1, p2] = fund_mac(K,ct)
%close all
R=1; Rx=10;al=4; ct=10000;              snrdb = [5: 5 : 30];     M=5;
           I=5;
m=4; n=m+1;
for k =  1:length(snrdb)
 sum1=0; sum2=0; sum3=0; sum4=0; sum5=0; 
 snr = 10^(snrdb(k)/10);     a=1+I; b=I/snr;eps1=(2^R-1)/snr;

 for i = 1 :ct
             
     h = complex(sqrt(0.5)*randn(M,1),sqrt(0.5)*randn(M,1));       
     h = abs(h).^2;
     h=sort(h,'ascend'); h1=h(m); h2=h(n);

     %%%%%%%%%%%%%%%%%%%% random use
     Rdx1(i) = log2(1+snr*h1);
     %%%%%%%%%%%%%%%%%%%%%%% end
    
     a2 = max(0, (h1-I/snr)/h1/(1+I));
  
     Rx2(i) = log2(1+snr*h2*a2);     

     
 end
    c1(k) = mean(Rx2);
    %analytical resultw
    w1=factorial(M)/factorial(m-1)/factorial(n-1-m)/factorial(M-n);
    w2=(1-2*a2)/a2^2;
    w3=factorial(M)/factorial(n-1)/factorial(M-n);
    w4=factorial(M)/factorial(m-1)/factorial(M-m);
    a=1+I; b=I/snr;eps1=(2^R-1)/snr;
%     
%     sum1 = 0;
%     upperx = 10;
%     stepx = (upperx-b)/5000;
%     xi = [b+stepx/1: stepx : upperx];
%     for i = 1 : length(xi)
%         x = xi(i);
%         fx= exp(-x); Fx = 1 - exp(-x);
%         sum1 = sum1 + w1/(M-n+1) * fx*Fx^(m-1)* ( log2(1+(x-b)/a*snr)*(1-Fx)^(M-n+1) + ...
%             1/log(2)*exp(x^2*a/(x-b)/snr)*expint((M-n+1)*x + (M-n+1)*x*a/(x-b)/snr )   )*stepx;
%     end
%     c2(k) = sum1;   
%     
cc=5000;
    sum1 = 0;
    upperx = 10;
    stepx = (upperx-b)/cc;
    xi = [b+stepx: stepx : upperx];
    for i = 1 : length(xi)
        x = xi(i);
        stepy = (upperx-x)/cc;
        yi = [x+stepx: stepy : upperx];
        fx= exp(-x); Fx = 1 - exp(-x);
        sum2=0;
        for j =  1 : length(yi)
            y = yi(j);
            fy= exp(-y); Fy = 1 - exp(-y);
            sum1 = sum1 + w1*log2(1+(x-b)/x/a*snr*y)*fx*Fx^(m-1)*fy*(1-Fy)^(M-n)*stepy*stepx;
        end
    end
    c3(k) = sum1;    
    
end

plot(snrdb,c1,snrdb,c3)%,snrdb,c3)
 