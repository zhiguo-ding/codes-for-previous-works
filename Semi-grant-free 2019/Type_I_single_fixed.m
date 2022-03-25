clear all 
 snrdb = [0: 5 : 40]; ct=50000;
 M=20; N=20;
 
 Ri=0.1;
  R0 = 0.5; %snr0 = (2^R0-1)/h0;
   
   
 for k =  1:length(snrdb)
    snr = 10^((snrdb(k))/10);  %snr0 = snr;
    snr0=  10^((snrdb(k))/10); 
    eps0 = (2^R0-1)/snr0;
    tau = 0.5;%1/snr;         
        
   sum1=0; sum2=0; sum3=0;sum4=0;  
  for ix = 1 :ct       
      h0 = abs(complex(sqrt(0.5)*randn(1,1),sqrt(0.5)*randn(1,1))).^2;
      epsx = (h0*snr0/snr/(2^R0-1))-1/snr; 
      
     nxx = 0;
     while nxx<N 
         hx = abs(complex(sqrt(0.5)*randn(1,1),sqrt(0.5)*randn(1,1))).^2;
         if hx<tau
             nxx = nxx+1;
             h(nxx) = hx;
         end
     end
     h = sort(h, 'ascend');     
     
     if log2(1+h0*snr0/(h(1)*snr +1 ))>R0 & log2(1+snr*h(1))>Ri
         sum1 = sum1+1;
     end
      
     
  end
  p1(k) = 1-sum1/ct;
  
  %analytical results
  sumt =0;
  eps0 = 2^(R0)-1; eps1 = 2^(Ri)-1;
  if eps1/snr<tau      
      for p = 0:N
          sumt = sumt + factorial(N)/factorial(p)/factorial(N-p)*(-1)^(N-p)/(1-exp(-tau))^N...
              *exp(p/snr - tau*(N-p))/(1+p/snr/eps0*snr0)...
              *(exp(-(eps0/snr0+p/snr)) - exp(-(eps0/snr0+p/snr)*(1+snr*tau)) );
      end
      pa(k) = 1 -  (exp(-eps1/snr)-exp(-tau))^N*exp(-eps0/snr0)/(1-exp(-tau))^N +sumt;  
  else
      pa(k) = 1;
  end
    
end
semilogy(snrdb,p1,snrdb,pa)% ,snrdb,  (10.^snrdb/10).^(-1)  )