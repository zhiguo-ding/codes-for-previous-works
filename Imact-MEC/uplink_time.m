clear all%function [p1, p2] = fund_mac(K,ct)  
 snrdb = [0: 10 : 50]; ct=500000;
 M=5; m=1; n=5; eta=1; %pn/pm
 
 for k =  1:length(snrdb)
    snr = 10^((snrdb(k))/10); pm=snr; pn=eta*pm;
   sum1=0; sum2=0; sum3=0;sum4=0;  
  for ix = 1 :ct       
     h = abs(complex(sqrt(0.5)*randn(M,1),sqrt(0.5)*randn(M,1))).^2;
     h = sort(h, 'ascend');

     % two rates 
     rn = log2(1+pn*h(n)/(1+pm*h(m)));
     rm = log2(1+pm*h(m));  
    if rn>rm
        sum1 = sum1+1;
    end 
%     if h(m)>(pn-pm)/pm^2 & h(n)>pm/pn*h(m)+pm^2/pn*h(m)^2 
%         sum2 = sum2+1;
%     end
    
  end
  p1(k) = sum1/ct;
%   p2(k) = sum2/ct;
  
    %analytical
    cmn = factorial(M)/factorial(m-1)/factorial(n-1-m)/factorial(M-n);
    suma1 = 0;
    for p=0: n-1-m
        cp = factorial(n-1-m)/factorial(p)/factorial(n-1-m-p)*(-1)^(n-1-m-p);
        suma2 = 0;
        for l=0:m-1
            cl = factorial(m-1)/factorial(l)/factorial(m-1-l)*(-1)^l;
            a = pm^2/pn*(M-m-p);
            b =  p+l+1+(M-m-p)*pm/pn;
            sum2x = 0;
            stepx = 0.0001;
            xx = [(max(0,pn-pm))/pm^2:stepx:20];
            for t = 1:length(xx)
                x = xx(t);
                sum2x = sum2x + cl*exp(-a*(x+b/2/a)^2+b^2/4/a)*stepx;
            end
            suma2 = suma2 + sum2x;
            
%            suma2 = suma2 + cl * exp(b^2/4/a)*sqrt(pi)/2/sqrt(a)...
%                *(1 - erf( (pn-pm)/pm^2 + b*sqrt(a)/2/a));
        end
        suma1 = suma1 + cmn*cp/(M-m-p)*suma2;
    end
        sumat2 = 0;
        for l=0:m-1
            cl = factorial(m-1)/factorial(l)/factorial(m-1-l)*(-1)^l;
            sumat2 = sumat2 + cl*exp( -(M-m+l+1)*(max(0,pn-pm))/pm^2 )/(M-m+l+1);
        end

    pa(k) = suma1+1-sumat2*factorial(M)/factorial(m-1)/factorial(M-m);
%     paa(k) = suma1;
    
    %%% approximiation with EVEN NUMBER M
    cmn = factorial(M)/factorial(m-1)/factorial(n-1-m)/factorial(M-n);
    suma1x = 0;
    for p=0: n-1-m
        a = pm^2/pn*(M-m-p);
        cp = factorial(n-1-m)/factorial(p)/factorial(n-1-m-p)*(-1)^(n-1-m-p);
        lam = p+1+(M-m-p)/eta;
        if mod(m,2)==0
            qsum = 0;             
                sumlx=0;
                for lx=0:m-1
                    sumlx = sumlx + factorial(m-1)/factorial(lx)/factorial(m-1-lx)*(-1)^lx*lx^m;
                end                 
 
            q1 = sqrt(pi)/ factorial(m/2)/2^(m+1)/a^((m+1)/2)*...
                (m*(lam)*(-1)^(m-1)*factorial(m-1) +  sumlx  );
        else
            q1 =  sqrt(pi)/ factorial((m-1)/2)/2^(m)/a^((m)/2)*(-1)^(m-1)*factorial(m-1);
        end;
        
        if mod(m,2)==0
            q2 =  1/(prod([m-1:-2:1]))/2^(m/2)/a^(m/2)*(-1)^(m-1)*factorial(m-1);
            
        else
            qsum = 0;             
                sumlx=0;
                for lx=0:m-1
                    sumlx = sumlx + factorial(m-1)/factorial(lx)/factorial(m-1-lx)*(-1)^lx*lx^m;
                end                 
 
            q2 = 1/(prod([m:-2:1]))/2^((m+1)/2)/a^((m+1)/2)*...
                (m*(lam)*(-1)^(m-1)*factorial(m-1) +  sumlx  );
        end;
 
        suma1x = suma1x + cmn*cp/(M-m-p) * ( q1- q2);
    end
    pax(k) = suma1x;% +factorial(M)/factorial(M-m)*(eta-1)^m/factorial(m)/pm^m;
    
    %low SNR approximation
    cmn = factorial(M)/factorial(m-1)/factorial(n-1-m)/factorial(M-n);
    suma3 = 0;
    for p=0: n-1-m
        cp = factorial(n-1-m)/factorial(p)/factorial(n-1-m-p)*(-1)^(n-1-m-p);
        suma4 = 0;
        for l=0:m-1
            cl = factorial(m-1)/factorial(l)/factorial(m-1-l)*(-1)^l;
            a = snr*(M-m-p);
            b = M-m+l+1;
            sumkk = 0;
            for kk = 0: 10
                sumkk = sumkk + (-1)^kk*(prod([2*kk:-2:1]))*2^kk*a^kk/b^(2*kk);
            end
            suma4 = suma4 + cl /b*sumkk;
        end
        suma3 = suma3 + cmn*cp/(M-m-p)*suma4;
    end
    pa2(k) = suma3;
    
end
semilogy(snrdb,p1,snrdb,pa ,snrdb, pax  )