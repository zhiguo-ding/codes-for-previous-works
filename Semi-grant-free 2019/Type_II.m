clear all 
 snrdb = [0: 5 : 40]; ct=5000;
 M=10; 
 
 
  R0 = 1; %snr0 = (2^R0-1)/h0;
  Ri=0.5;
   
   
 for k =  1:length(snrdb)
    snr = 10^((snrdb(k))/10);  %snr0 = snr;
    snr0=10^((snrdb(5))/10); 
    eps0 = (2^R0-1)/snr0;          %%%%% BE CAREFUL WITH snr0
    tau = 1;%/snr;         
        
   sum1=0; sum2=0; sum3=0;sum4=0;  
  for ix = 1 :ct       
      h0 = abs(complex(sqrt(0.5)*randn(1,1),sqrt(0.5)*randn(1,1))).^2;
      
     h = abs(complex(sqrt(0.5)*randn(M,1),sqrt(0.5)*randn(M,1))).^2;
     h = sort(h, 'ascend');     
     
     n = sum(sign(h-tau)+1)/2;
     nk(ix)=n;
     xm = sum(h(M-n+1:M));
     
     Rs = log2(1+xm*snr/(h0*snr0 +1 ));
     
     if (Rs==0|Rs>n*Ri) & log2(1+h0*snr0)>R0
         sum1 = sum1 +1;
     end
     
     Rall = log2(1+sum(h)*snr/(h0*snr0 +1 ));
     if (Rall==0|Rall>M*Ri) & log2(1+h0*snr0)>R0
         sum2 = sum2 +1;
     end
     
  end
  p1(k) = 1- sum1/ct;
  p2(k) = 1- sum2/ct;
  nkk(k) = mean(nk);
   
    %analytical    
    suma1 = 0;
    for n=1:M
        eps_s = 2^(n*Ri)-1;
        tau_n = (n*tau/eps_s*snr-1)/snr0;
        tau_bar = max(tau_n, eps0);
        Q5 =0;
        for l = 1 : n-1
            Q5 = Q5 + eps_s^l*snr^(-l)*snr0^l/factorial(l)...
                *exp(-tau_n)*gammainc((tau_bar-tau_n)*(1+eps_s*snr^(-1)*snr0), l+1, 'upper')*factorial(l)...
                /(1+eps_s/snr*snr0)^(l+1);                
        end
        Q5 = Q5 +exp(-eps0)-exp(-tau_bar) + exp(-tau_n - (tau_bar-tau_n)*(1+eps_s*snr^(-1)*snr0))/(1+eps_s*snr^(-1)*snr0);
        suma1 = suma1 + factorial(M)/factorial(n)/factorial(M-n)*exp(-n*tau)*...
            (1-exp(-tau))^(M-n)* Q5;        
    end
    pa(k) = 1 - suma1 - (1-exp(-tau))^M*exp(-eps0); 
     
    %%approximation
    eps_s1 = 2^(n*Ri)-1;
    paa(k) = 1 - exp(-eps0) +  exp(1/snr0)*snr0*M*exp(-tau)*...
        (1-exp(-tau))^(M-1)*eps_s1/snr/exp(tau/eps_s1*snr/snr0);
    
    nkk(k) = mean(nk);
end
semilogy(snrdb,p2,snrdb, 1-exp(-eps0)*ones(1,length(snrdb)),snrdb,p1,snrdb,pa,snrdb,min(1,paa))% ,snrdb,  (10.^snrdb/10).^(-1)  )