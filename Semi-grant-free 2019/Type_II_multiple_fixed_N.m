clear all 
 snrdb = [0 :5 : 30]; ct=5000;
 M=10;
   N=5; i =M;
 snr0=10^((20)/10);
 
  R0 = 1; %snr0 = (2^R0-1)/h0;
  Ri=0.5;
   eps0 = (2^R0-1)/snr0;
   
 for k = 1: length(snrdb)
    snr = 10^((snrdb(k))/10);  %snr0 = snr;
    tau = 1;%/snr;         
        
   sum1=0; sum2=0; sum3=0;sum4=0;  
  for ix = 1 :ct       
      h0 = abs(complex(sqrt(0.5)*randn(1,1),sqrt(0.5)*randn(1,1))).^2;
      
     nxx = 0;
     while nxx<N 
         hx = abs(complex(sqrt(0.5)*randn(1,1),sqrt(0.5)*randn(1,1))).^2;
         if hx>tau
             nxx = nxx+1;
             h(nxx) = hx;
         end
     end
     h = sort(h, 'ascend');    
     
     %the best user's rate - single user
     if log2(1+h(i-(M-N+1)+1)*snr/(snr*sum(h(1:i-(M-N+1)  ))+h0*snr0 +1 ))<Ri
         sum1 = sum1+1;
     end
     
  end
  p1(k) = sum1/ct;
   
 end
semilogy(snrdb,p1  )% ,snrdb,  (10.^snrdb/10).^(-1)  )