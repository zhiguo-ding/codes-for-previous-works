clear all 
 snrdb = [0: 5 : 30]; ct=50000;
 %M=20; 
 N=5;
 %snr0=10^((20)/10);
 Ri=0.6;
   i=N;
   eps1 = 2^(Ri)-1;
   
 for k = 1:length(snrdb)
    snr = 10^((snrdb(k))/10);  %snr0 = snr;
    tau = 0.5;%1/snr;         
        
   sum1=0; sum2=0; sum3=0;sum4=0;  
  for ix = 1 :ct       
      h0 = abs(complex(sqrt(0.5)*randn(1,1),sqrt(0.5)*randn(1,1))).^2;
      
     nxx = 0;
     while nxx<N 
         hx = abs(complex(sqrt(0.5)*randn(1,1),sqrt(0.5)*randn(1,1))).^2;
         if hx<tau
             nxx = nxx+1;
             h(nxx) = hx;
         end
     end
     h = sort(h, 'ascend');     
     
     if log2(1+h(i)*snr/(sum(h(1:i-1))*snr +1 ))>Ri  
         sum1 = sum1+1;
     end
      
     
  end
  p1(k) = 1-sum1/ct;
   
  %analytical results
  n100 = 500;
  if tau<eps1/snr
      stepx = tau/n100;
      xall = [0:stepx : tau];
      sumx1 = 0;
      for ix = 1 : length(xall)
          x = xall(ix);
          fx = factorial(N)/factorial(i-1)/factorial(N-i)...
              *exp(-x)*(1-exp(-x))^(i-1)*(exp(-x)-exp(-tau))^(N-i)...
              /(1-exp(-tau))^N;
          sumx1 = sumx1 + fx*stepx;
      end
      pa(k) = sumx1;
  else
      stepx = eps1/snr/n100;
      xall = [0:stepx : eps1/snr];
      sumx1 = 0;
      for ix = 1 : length(xall)
          x = xall(ix);
          fx = factorial(N)/factorial(i-1)/factorial(N-i)...
              *exp(-x)*(1-exp(-x))^(i-1)*(exp(-x)-exp(-tau))^(N-i)...
              /(1-exp(-tau))^N;
          sumx1 = sumx1 + fx*stepx;
      end
      sumx2 = 0;
      stepx = (tau-eps1/snr)/n100;
      xall = [eps1/snr:stepx : tau];
      for ix = 1 : length(xall)
          x = xall(ix);
          fx = factorial(N)/factorial(i-1)/factorial(N-i)...
              *exp(-x)*(1-exp(-x))^(i-1)*(exp(-x)-exp(-tau))^(N-i)...
              /(1-exp(-tau))^N;
          
          stepy = (x*(i-1) - (x/eps1-1/snr))/n100;
          yall = [(x/eps1-1/snr):stepy : x*(i-1)];
          sumx3 = 0;
          for iy = 1: length(yall)
              y = yall(iy);
              if y>tau*(i-1)
                  cfdfd=0;
              end
              fy_temp = 0;
              for p = 0:i-1
                  fy_temp = fy_temp + factorial(i-1)/factorial(p)/factorial(i-1-p)*(-1)^p...
                      *exp(-p*x)/(1-exp(-x))^(i-1)/factorial(i-2)...
                      *(y-p*x)^(i-2)*exp(-(y-p*x))*max(0,sign(y-p*x)); %%%%% THIS IS WRONG
              end
              fy_tempx(iy) = fy_temp;
              sumx3 = sumx3 +  fy_temp *stepy;
          end
          
          sumx2 = sumx2 + fx*sumx3*stepx;
      end
      
      pa(k) = sumx1 + sumx2;
  end
    
end
semilogy(snrdb,p1,snrdb,pa )% ,snrdb,  (10.^snrdb/10).^(-1)  )