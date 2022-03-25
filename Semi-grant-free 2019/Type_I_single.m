clear all 
 snrdb = [-10: 10 : 30]; ct=100000;
 M=10; 
 snr0=10^((30)/10);
 
  R0 = 2; %snr0 = (2^R0-1)/h0;

   
 for k =  1:length(snrdb)
    snr = 10^((snrdb(k))/10);  %snr0 = snr;
    tau = 100/snr;         
        
   sum1=0; sum2=0; sum3=0;sum4=0;  
  for ix = 1 :ct       
      h0 = abs(complex(sqrt(0.5)*randn(1,1),sqrt(0.5)*randn(1,1))).^2;
      epsx = (h0*snr0/snr/(2^R0-1))-1/snr; 
      
     h = abs(complex(sqrt(0.5)*randn(M,1),sqrt(0.5)*randn(M,1))).^2;
     h = sort(h, 'ascend');     
     
     n = M-sum(sign(h-tau)+1)/2;
     nk(ix)=n;
     %xm = sum(h(1:n));
     if n == 0
         xm = 0;%sum(h(min(1,n)));
     else
         xm = sum(h(1));
     end
     
     Rx0 = log2(1+h0*snr0/(xm*snr +1 ));
     
     if Rx0< R0
         sum1 = sum1 +1;
     end
     
  end
  p1(k) = sum1/ct;
  nkk(k) = mean(nk);
    
    %analytical
    suma1 = 0;
    eps0 = 2^(R0)-1; 
    for n=1:M
        Q6 =0;
        for p=0:n
             Q6 = Q6 + factorial(n)/factorial(p)/factorial(n-p)*(-1)^(n-p)...
                 *exp(p/snr - (n-p)*tau )/(1-exp(-tau))^n...
                 *( exp(-( eps0/snr0+p/snr )) - exp( -(1+snr*tau)*(eps0/snr0+p/snr)))...
                 /(1+p/snr/eps0*snr0);
        end
        suma1 = suma1 + factorial(M)/factorial(n)/factorial(M-n)*exp(-(M-n)*tau)*...
            (1-exp(-tau))^n* (1-exp(-eps0/snr0)+Q6);        
    end
    pa(k) = suma1 +  exp(-M*tau)*(1-exp(-eps0/snr0));  

 end
semilogy(snrdb,p1,snrdb,pa )% ,snrdb,  (10.^snrdb/10).^(-1)  )