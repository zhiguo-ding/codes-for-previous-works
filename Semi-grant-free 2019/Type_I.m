
%%%%% REMEMBER EPS0 IS EPS0/SNR0 in this code


clear all 
 snrdb = [0: 5 : 30]; ct=500000;
 M=20; 
 %snr0=10^((10)/10);
 
  R0 = 1; %snr0 = (2^R0-1)/h0;
   
   
 for k =  1:length(snrdb)
    snr = 10^((snrdb(1))/10);  
    
    snr0 = 10^((snrdb(k))/10);  
    eps0 = (2^R0-1)/snr0;
    tau = 1/snr0;         
        
   sum1=0; sum2=0; sum3=0;sum4=0;  
  for ix = 1 :ct       
      h0 = abs(complex(sqrt(0.5)*randn(1,1),sqrt(0.5)*randn(1,1))).^2;
      epsx = (h0*snr0/snr/(2^R0-1))-1/snr; 
      
     h = abs(complex(sqrt(0.5)*randn(M,1),sqrt(0.5)*randn(M,1))).^2;
     h = sort(h, 'ascend');     
     
     n = M-sum(sign(h-tau)+1)/2;
     nk(ix)=n;
     xm = sum(h(1:n));
     
     Rx0 = log2(1+h0*snr0/(xm*snr +1 ));
     
     if Rx0< R0
         sum1 = sum1 +1;
     end
     if log2(1+snr0*h0)<R0
         sum2 = sum2 +1;
     end
     
     if log2(1+h0*snr0/(sum(h)*snr +1 ))< R0
         sum3 = sum3 +1;
     end
  end
  p1(k) = sum1/ct;
  nkk(k) = mean(nk);
  p2(k) = sum2/ct;
  p3(k) = sum3/ct;
   
    %analytical
    suma1 = 0;
    for n=1:M
        sumtemp =0;
        for p=0:n
            sumtemp2 = 0;
            for l=0:n-1
                sumtemp2 = sumtemp2 + eps0*snr*exp(-eps0*(1+snr*p*tau))...
                    /(1+eps0*snr)^(l+1);
            end
            sumtemp = sumtemp + factorial(n)/factorial(n-p)/factorial(p)*(-1)^p...
                *exp(-p*tau)/(1-exp(-tau))^n *(sumtemp2 +exp(-eps0) - exp(-eps0-eps0*snr*p*tau)); 
        end
        suma1 = suma1 + factorial(M)/factorial(n)/factorial(M-n)*exp(-(M-n)*tau)*...
            (1-exp(-tau))^n* sumtemp;        
    end
    pa(k) = suma1 +  (1-exp(-eps0));  
     
    %%% fix bar{P}, increase P0 and tau=1/snr0
%     paa(k) = tau*eps0*snr*M*exp(-(M-1)*tau) +eps0;

%     %%% fix bar{P}, fix tau,  increase P0  
    summm1 = 0;
    for n = 1 : M
        summm1 = summm1 + factorial(M)/factorial(n-1)/factorial(M-n)*...
            exp(-(M-n)*tau)*(1-exp(-tau))^(n-1)*...
            ( 1 - exp(-tau) - tau*exp(-tau) );
    end
    paa(k) = eps0*snr*summm1 + eps0;

%     %%% fix bar{P}, fix P0,  reduce snr  
%     summm1 = 0;
%     for n = 1 : M
%         sump2 = 0;
%         for p = 0: n
%             sump2 = sump2 +factorial(n)/factorial(p)/factorial(n-p)...
%                 *(-1)^p*exp(-p*tau)*( exp(-eps0) +p*tau);
%         end
%         summm1 = summm1 + factorial(M)/factorial(n-1)/factorial(M-n)*...
%             exp(-(M-n)*tau)* sump2; 
%     end
%     paa(k) = eps0*snr*summm1 + 1 - exp(-eps0);

 end
semilogy(snrdb,p3,snrdb,p2,snrdb,p1,'-*',snrdb,pa,'-d',snrdb,min(1,paa),'-+')% ,snrdb,  (10.^snrdb/10).^(-1)  )