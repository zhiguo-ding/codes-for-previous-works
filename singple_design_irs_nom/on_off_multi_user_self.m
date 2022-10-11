clear all 
 snrdb = [0: 5 : 30]; ct=2000;
 M=4; N=40; K=3;
 al1 = 0.8; al2=0.2;
 R = 1;
 P=N; Q=N/P;
 V=kron(eye(P), ones(Q,1))/sqrt(Q);
 
for ki =  1:length(snrdb)
     snr = 10^((snrdb(ki))/10);
        sum0=0; sum1=0; sum2=0; sum3=0;sum4=0;   sum5=0;
  for ix = 1 :ct            
     %%%%%%%%%%%%%%%%%%%%%%%%% High mobility user
     h1 = complex(sqrt(0.5)*randn(N,1),sqrt(0.5)*randn(N,1));
     h2 = complex(sqrt(0.5)*randn(N,K-1),sqrt(0.5)*randn(N,K-1));
      
     D1 = diag( complex(sqrt(0.5)*randn(N,1),sqrt(0.5)*randn(N,1)) );  
     
      for fi =1 : P
          fvector = V(:,fi);
          SINR(fi) = abs(fvector'*D1*h1)^2*al1/(abs(fvector'*D1*h1)^2*al2 ...
              + sum(abs(fvector'*D1*h2).^2)+ 1/snr);
      end
      if max(SINR)<2^R -1
          sum1 = sum1 +1;
      end
       
  end
  
  p(ki) = sum1/ct;
   
  %analytical result
  epidot =2^R -1;
  tau = al1 - al2*epidot;  
  tz = 2*sqrt(epidot/snr/tau);
  pa(ki) =  (1 -   tz * besselk(1,tz)...
      /(1+ epidot/tau)^(K-1))^N;
 
  %approximated results
  paa(ki) = (1 - 1 /(1+ epidot/tau)^(K-1))^N;
  
end
semilogy(snrdb, p ,snrdb,pa,snrdb,paa )
%semilogy(snrdb, p  ,snrdb,po,snrdb,pdft)
 
  