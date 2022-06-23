clear all%function [p1, p2] = fund_mac(K,ct)
%close all
R=0; Rx=10;al=4; ct=5000;              snrdb = [5: 5 : 30];     M=5;
           
m=1; n=3;
for k =   1:length(snrdb)
 sum1=0; sum2=0; sum3=0; sum4=0; sum5=0; 
 snr = 10^(snrdb(k)/10);
 for i = 1 :ct
             
     h = complex(sqrt(0.5)*randn(M,1),sqrt(0.5)*randn(M,1));       
     h = abs(h).^2;
     h=sort(h,'ascend'); h1=h(m); h2=h(n);

     %%%%%%%%%%%%%%%%%%%% random use
     Rdx(i) = log2(1+snr*h1)/2+log2(1+snr*h2)/2;
     %%%%%%%%%%%%%%%%%%%%%%% end
    
     a1=6/11;%h2/(h1+h2);
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
    w1=factorial(M)/factorial(m-1)/factorial(n-1-m)/factorial(M-n);
    w2=(1-2*a2)/a2^2;
    w3=factorial(M)/factorial(n-1)/factorial(M-n);
    
    sumx1=0;
    sumx2=0;
    for i =0: n-1-m
        sumx1 = sumx1 + w1*factorial(n-1-m)/factorial(i)/factorial(n-1-m-i)*(-1)^i/(m+i)*...
            intg1(w2,n,m,i,snr,M);
        sumx2 = sumx2 + factorial(n-1-m)/factorial(i)/factorial(n-1-m-i)*(-1)^i/(m+i)*...
            intg2(w2,n,m,i,snr,M);
    end
    for j=0: n-1
        sumx1 = sumx1 + w3/snr*factorial(n-1)/factorial(j)/factorial(n-1-j)*(-1)^j*snr/...
            (M-n+j+1)*exp(-(M-n+j+1)*w2/snr);
    end
    p2(k) = 1-sumx1;
    
    p3(k) = w3*w2^n/n/snr^n - w1/snr^n*sumx2   ;
  
end
semilogy( snrdb,p2 )

% figure
% plot(snrdb,c1,snrdb,c2)