clear all%function [p1, p2] = fund_mac(K,ct)
%close all
R=1; Rx=10;al=4; ct=50000;              snrdb = [5: 5 : 30];     M=5;
           I=1;
m=1; n=M-0;
for k =  1:length(snrdb)
 sum1=0; sum2=0; sum3=0; sum4=0; sum5=0; 
 snr = 10^(snrdb(k)/10);     a=1+I; b=I/snr;eps1=(2^R-1)/snr;

  
    p1(k) = sum1/ct; 
    p2(k) = sum2/ct; 
    px(k) = sum3/ct;
    %analytical resultw
    w1=factorial(M)/factorial(m-1)/factorial(n-1-m)/factorial(M-n);
    %w2=(1-2*a2)/a2^2;
    w3=factorial(M)/factorial(n-1)/factorial(M-n);
    w4=factorial(M)/factorial(m-1)/factorial(M-m);
    a=1+I; b=I/snr;eps1=(2^R-1)/snr;
        
     
    sum1 = 0;
    for i =0 : n-1-m 
        sum1 = sum1 + w1*factorial(n-m-1)/factorial(i)/factorial(n-1-m-i)*(-1)^i...
            *intg31(0,n,m,i,snr,M,b,eps1,a);
    end
    
    for i =0 : n-1-m 
        sum1 = sum1 + w1*factorial(n-m-1)/factorial(i)/factorial(n-1-m-i)*(-1)^i...
            *intg3(0,n,m,i,snr,M,b,eps1,a);
    end
    for i =0 : n-1-m 
        sum1 = sum1 + w1*factorial(n-m-1)/factorial(i)/factorial(n-1-m-i)*(-1)^i...
            *intg4(0,n,m,i,snr,M,b,eps1,a);
    end
    
    for i =0:m-1
        sum1 = sum1 + w4*factorial(m-1)/factorial(i)/factorial(m-1-i)*(-1)^i/(M-m+i+1)...
            *(1 - exp(- (M-m+i+1)*b));
    end

    p3(k) = sum1;
    
    p4(k) = min(1,100/snr^((m)));
end
semilogy( snrdb,p3,'-*' )
