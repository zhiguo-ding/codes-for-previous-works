clear all%function [p1, p2] = fund_mac(K,ct)  
 snrdb = [0: 5 : 30]; ct=500000;
 M=5; m=1; n=5;eta=5; %pn/pm
 beta=1/3;
 
 for k = 1: length(snrdb)
    snr = 10^((snrdb(k))/10);pm=snr; pn=pm*eta;
   sum1=0; sum2=0; sum3=0;sum4=0;  
  for ix = 1 :ct       
     h = abs(complex(sqrt(0.5)*randn(M,1),sqrt(0.5)*randn(M,1))).^2;
     h = sort(h, 'ascend');

     % two rates 
     rn1 = log2(1+pn*h(n));
     rn2 = log2(1+beta*pn*h(n)/(pm*h(m)+1)) +log2(1+beta*pn*h(n));  
    if rn1>rn2
        sum1 = sum1+1;
    end 
    
    kka=(1-2*beta)/( beta^2*pn - (1-beta)*pm ); 
    tth=((1-beta)*(1+pm*h(m))-beta)/beta^2/pn;
    if h(m)<kka&h(n)>tth% &h(m)<kka
        sum2 = sum2 +1;
    end
    
  end
    p1(k) = sum1/ct;
    p2(k) = sum2/ct;
    
    
    %analytical
    if (1-beta)*pm>beta^2*pn
        cmn = factorial(M)/factorial(m-1)/factorial(n-1-m)/factorial(M-n);
        suma1 = 0;
        for p=0: n-1-m
            cp = factorial(n-1-m)/factorial(p)/factorial(n-1-m-p)*(-1)^(n-1-m-p);
            suma2 = 0;
            for l=0:m-1
                cl = factorial(m-1)/factorial(l)/factorial(m-1-l)*(-1)^l;
                a = pm*(1-beta)*(M-m-p)/beta^2/pn+p+l+1; 
                suma2 = suma2 + cl/(a)*exp(-(M-m-p)*(1-2*beta)/beta^2/pn);
            end
            suma1 = suma1 + cmn*cp/(M-m-p)*suma2;
        end
        pa1(k) =  1 - suma1;
    else
        cmn = factorial(M)/factorial(m-1)/factorial(n-1-m)/factorial(M-n);
        suma1 = 0;
        for p=0: n-1-m
            cp = factorial(n-1-m)/factorial(p)/factorial(n-1-m-p)*(-1)^(n-1-m-p);
            suma2 = 0;
            for l=0:m-1
                cl = factorial(m-1)/factorial(l)/factorial(m-1-l)*(-1)^l;
                a = pm*(1-beta)*(M-m-p)/beta^2/pn+p+l+1; 
                ka = (1-2*beta)/( beta^2*pn - (1-beta)*pm );
                suma2 = suma2 + cl/(a)*exp(-(M-m-p)*(1-2*beta)/beta^2/pn)*...
                    (1-exp(-a*ka));% + cl/(M-m+l+1)*exp(-ka*(M-m+l+1));
            end
            suma1 = suma1 + cmn*cp/(M-m-p)*suma2;
        end
        sumat2 = 0;
        for l=0:m-1
            cl = factorial(m-1)/factorial(l)/factorial(m-1-l)*(-1)^l;
            sumat2 = sumat2 + cl*exp( -(M-m+l+1)*ka )/(M-m+l+1);
        end
        pa1(k) =  1 - factorial(M)/factorial(m-1)/factorial(M-m)*sumat2-suma1;
    end
    
    
    %approximation user m
    p2x(k) = factorial(M)/factorial(M-m)*2^(2*m)*(1-beta)^m/factorial(m)/beta^(2*m)/snr^m;
    
        %analytical
    cmn = factorial(M)/factorial(m-1)/factorial(n-1-m)/factorial(M-n);
    suma1 = 0;
    for p=0: n-1-m
        cp = factorial(n-1-m)/factorial(p)/factorial(n-1-m-p)*(-1)^(n-1-m-p);
        suma2 = 0;
        for l=0:m-1
            cl = factorial(m-1)/factorial(l)/factorial(m-1-l)*(-1)^l;
            a = (2-beta)*(M-m-p)/beta; 
            suma2 = suma2 + cl/(a+p+l+1)*1;
        end
        mm=[1:1:m];
        suma2 = factorial(m-1)/(prod(a+p+mm));
        suma1 = suma1 + cmn*cp/(M-m-p)*suma2;
    end
    p1x(k) = 1-suma1;
 end
semilogy(snrdb,p1,'-d',snrdb,pa1,'-*' )
%semilogy(snrdb,p1,'-d',snrdb,p2 ,snrdb,pa1 ,snrdb,pa2 ,snrdb,min(1,p1x),snrdb,min(1,p2x))