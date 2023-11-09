clear all%function [p1, p2] = fund_mac(K,ct)  
 snrdb = [0: 5 : 20]; ct=500000;
 M=5; m=2; n=5;
 
 NT = 1;
 eps = 2^NT-1; beta=1/10; N=20;
 
 for k = 1: length(snrdb)
    snr = 10^((snrdb(k))/10);
   sum1=0; sum2=0; sum3=0;sum4=0;  
  for ix = 1 :ct       
     h = abs(complex(sqrt(0.5)*randn(M,1),sqrt(0.5)*randn(M,1))).^2;
     h = sort(h, 'ascend');
     
     alphan = max(0,(snr*h(m)-eps)/snr/h(m)/(1+eps));
     % two rates 
     rn1 = log2(1+alphan*snr*h(n)) +log2(1+beta*snr*h(n));
     rn2 = log2(1+snr*h(n));  
    if rn1<rn2 
        sum1 = sum1+1;
    end 
        
  end
    p1(k) = sum1/ct;
    
    %analytical results
    cmn = factorial(M)/factorial(m-1)/factorial(n-1-m)/factorial(M-n);
    suma1 = 0;
    for p=0: n-1-m
        cp = factorial(n-1-m)/factorial(p)/factorial(n-1-m-p)*(-1)^(n-1-m-p);
        suma2 = 0;
        for i=1:N
            thi = cos((2*i-1)/2/N*pi);
            xi = (eps/2/beta/snr+eps/2/snr) +(eps/2/beta/snr-eps/2/snr)*thi;
            
            suma2 = suma2 + pi/N*(eps/2/beta/snr-eps/2/snr)*...
                sqrt(1-thi^2)*exp(-(p+1)*xi)*(1-exp(-xi))^(m-1)/(M-m-p)...
                *(exp(-(M-m-p)*xi) - exp(-(M-m-p)* ( snr*xi*((1-beta)*(1+eps)-1) +eps)...
                /(snr*beta*(snr*xi-eps)))  );
        end
        suma1 = suma1 + cmn*cp*suma2;
    end
    %%%the first term
    suma11 = 0;
    for l=0: m-1
        cl = factorial(m-1)/factorial(l)/factorial(m-1-l)*(-1)^l;
        suma11 = suma11 + factorial(M)/factorial(m-1)/factorial(M-m)*cl...
            *exp( -(M-m+l+1)*eps/snr) /(M-m+l+1);
    end
    pa1(k) = 1-suma11+suma1;
    
    
end
semilogy(snrdb,p1,'-d',snrdb,pa1)