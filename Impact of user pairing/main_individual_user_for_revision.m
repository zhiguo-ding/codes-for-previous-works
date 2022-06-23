clear all%function [p1, p2] = fund_mac(K,ct)
%close all
R=0; Rx=10;al=4; ct=10000;              snrdb = [5: 5 : 30];     M=5;
           
m=1; n=M-2;
for k =   1:length(snrdb)
 sum1=0; sum2=0; sum3=0; sum4=0; sum5=0; 
 snr = 10^(snrdb(k)/10);
 for i = 1 :ct
             
     h = complex(sqrt(0.5)*randn(M,1),sqrt(0.5)*randn(M,1));       
     h = abs(h).^2;
     h=sort(h,'ascend'); h1=h(m); h2=h(n);

     %%%%%%%%%%%%%%%%%%%% random use
     Rdx1 = log2(1+snr*h1)/2;
     Rdx2 = log2(1+snr*h2)/2;
     %%%%%%%%%%%%%%%%%%%%%%% end
    
     a1=10/11;%h2/(h1+h2);
     a2=1-a1;%0.9;%h1/(h1+h2);
  
     Rx1 = log2( 1+snr*h1*a1/(snr*h1*a2+1) );     
     Rx2 = log2(1+snr*h2*a2);
     
     if Rx1>Rdx1
         sum1 = sum1+1;
     end
     if Rx2<Rdx2
         sum2 = sum2+1;
     end
 end
    p1(k) = sum1/ct;
    p2(k) = sum2/ct;
    
    %analytical resultw
    w1=factorial(M)/factorial(m-1)/factorial(n-1-m)/factorial(M-n);
    w2=(1-2*a2)/a2^2;
    w3=factorial(M)/factorial(n-1)/factorial(M-n);
    w5=factorial(M)/factorial(m-1)/factorial(M-m);
    
    sumx1=0;
    for i =0: m-1
        sumx1 = sumx1 + w5*factorial(m-1)/factorial(i)/factorial(m-1-i)*(-1)^i/(M-m+i+1)*...
            (1-exp(- (1-2*a2)*(M-m+i+1)/snr/a2^2));
    end
    p3(k) = sumx1;
    sumx2=0;
    for j=0: n-1
        sumx2 = sumx2 + w3*factorial(n-1)/factorial(j)/factorial(n-1-j)*(-1)^j/...
            (M-n+j+1)*(1 - exp(-(1-2*a2)*(M-n+j+1)/snr/a2^2));
    end
    p4(k) = sumx2;
  
end
semilogy(snrdb,p3,snrdb,p4)

% figure
% plot(snrdb,c1,snrdb,c2)