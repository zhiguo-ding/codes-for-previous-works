clear all%function [p1, p2] = fund_mac(K,ct)
%close all
R=0.5; Rx=10;al=4; ct=500000;              snrdb = [5: 5 : 40];     M=5;
           
m=2; n=M;
for k =  1:length(snrdb)
 sum1=0; sum2=0; sum3=0; sum4=0; sum5=0; 
 snr = 10^(snrdb(k)/10);
 for i = 1 :ct
             
     h = complex(sqrt(0.5)*randn(M,1),sqrt(0.5)*randn(M,1));       
     h = abs(h).^2;
     h=sort(h,'ascend'); h1=h(m); h2=h(n);

     %%%%%%%%%%%%%%%%%%%% random use
     Rdx(i) = log2(1+snr*h1)/2+log2(1+snr*h2)/2;
     %%%%%%%%%%%%%%%%%%%%%%% end
    
     a1=4/5;%h2/(h1+h2);
     a2=1-a1;%0.9;%h1/(h1+h2);
  
     Rx1 = log2( 1+snr*h1*a1/(snr*h1*a2+1) );     
     Rx2 = log2(1+snr*h2*a2);
     
     Rxn(i) = Rx1+Rx2;
     
     if Rxn(i)-Rdx(i)<R
         sum1 = sum1+1;
     end
 end
    p1(k) = sum1/ct;
    
    %analytical resultw
    sum2=0;
    for j1= 0: m-1
        sum3=0;
        for j2 = 0:n-m-1
            t1 = j1-j2+n-m;
            t2 = M-n+1+j2;
            sum3 = sum3+ (-1)^(j1+j2)* factorial(m-1)/factorial(j1)/factorial(m-1-j1)...
                *factorial(n-m-1)/factorial(j2)/factorial(n-m-1-j2)/t1...
                *(1/(t2+2^(-2*R)*t1) - 1/(t2+t1) );
        end
        sum2 = sum2 + factorial(M)/factorial(m-1)/factorial(n-m-1)/factorial(M-n)*sum3;
    end
    p2(k) = sum2;
  
end
semilogy(snrdb,p1,'-*',snrdb,p2)

% figure
% plot(snrdb,c1,snrdb,c2)