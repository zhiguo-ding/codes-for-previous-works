clear all 
 snrdb = [0 :10 : 30]; ct=5000;
 M=10; 
 snr0=10^((20)/10);
 
  R0 = 1; %snr0 = (2^R0-1)/h0;
  Ri=0.2;
   eps0 = (2^R0-1)/snr0;
   
 for k =  1:length(snrdb)
    snr = 10^((snrdb(k))/10);  %snr0 = snr;
    tau = 1;%/snr;         
        
   sum1=0; sum2=0; sum3=0;sum4=0;  
  for ix = 1 :ct       
      h0 = abs(complex(sqrt(0.5)*randn(1,1),sqrt(0.5)*randn(1,1))).^2;
      
     h = abs(complex(sqrt(0.5)*randn(M,1),sqrt(0.5)*randn(M,1))).^2;
     h = sort(h, 'ascend');     
     
     n = sum(sign(h-tau)+1)/2;
     nk(ix)=n;
     if n==0
         xm = 0;
     else
         xm = sum(h(M));%sum(h(M-n+1:M));
     end
     
     Rs = log2(1+xm*snr/(h0*snr0 +1 ));
     
     %use the best
     if (Rs==0|Rs>Ri) & log2(1+h0*snr0)>R0
         sum1 = sum1 +1;
     end
 
     
  
 end
  p1(k) = 1- sum1/ct;
  nkk(k) = mean(nk);
  
    %analytical
    suma1 = 0;
    eps0 = 2^(R0)-1; eps1 = 2^(Ri)-1; 
    for n=1:M
        Q7x =0;
        thetah = max(eps0/snr0, snr/snr0/eps1*tau-1/snr0);
        for p=0:n
             Q7x = Q7x + factorial(n)/factorial(p)/factorial(n-p)*(-1)^(p)*exp(p*tau)...
                 *exp(-p/snr*eps1) *exp(-(1+p/snr*eps1*snr0)*thetah)/(1+p/snr*eps1*snr0);
        end
        if snr/snr0/eps1*tau -1/snr0>eps0/snr0
            Q7 = exp(-thetah) - Q7x + exp(-eps0/snr0) ...
                - exp(- (snr/snr0/eps1*tau - 1/snr0));
        else
            Q7 = exp(-thetah) -  Q7x ;
        end
        suma1 = suma1 + factorial(M)/factorial(n)/factorial(M-n)*exp(-n*tau)*...
            (1-exp(-tau))^(M-n)* Q7;        
    end
    pa(k) = 1-suma1 - (1- exp(-tau))^M*exp(-eps0/snr0);  
    
end
semilogy(snrdb,p1,snrdb,pa)% ,snrdb,  (10.^snrdb/10).^(-1)  )