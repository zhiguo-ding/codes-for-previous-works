clear all 
 snrdb = [0: 5 : 30]; ct=5000;
 M=4; N=2;
 al1 = 0.8; al2=0.2;
 R = 2;
 P=N; Q=N/P;
 V=kron(eye(P), ones(Q,1))/sqrt(Q);
 
for ki =  1:length(snrdb)
     snr = 10^((snrdb(ki))/10);
        sum0=0; sum1=0; sum2=0; sum3=0;sum4=0;   sum5=0;
  for ix = 1 :ct            
     %%%%%%%%%%%%%%%%%%%%%%%%% High mobility user
     h1 = complex(sqrt(0.5)*randn(N,1),sqrt(0.5)*randn(N,1));
     h2 = complex(sqrt(0.5)*randn(N,1),sqrt(0.5)*randn(N,1));
      
     D1 = diag( complex(sqrt(0.5)*randn(N,1),sqrt(0.5)*randn(N,1)) );
     
      for fi =1 : P
          fvector = V(:,fi);
          SINR(fi) = abs(fvector'*D1*h1)^2*al1/(abs(fvector'*D1*h1)^2*al2  + 1/snr);
      end
      if max(SINR)<2^R -1
          sum1 = sum1 +1;
      end
      %optimal
      fo = D1*h1/sqrt(h1'*D1'*D1*h1);
      SINRo = abs(fo'*D1*h1)^2*al1/(abs(fo'*D1*h1)^2*al2  + 1/snr);
      if  (SINRo)<2^R -1
          sum2 = sum2 +1;
      end
     
      %DFT
      d= dftmtx(N)/sqrt(N);
      for fi =1 : P
          fvectorx = d(:,fi);
          SINRx(fi) = abs(fvectorx'*D1*h1)^2*al1/(abs(fvectorx'*D1*h1)^2*al2  + 1/snr);
      end
      if max(SINRx)<2^R -1
          sum3 = sum3 +1;
      end
  end
  
  p(ki) = sum1/ct;
  po(ki) = sum2/ct;
  pdft(ki) = sum3/ct;
  
  %analytical result
  epidot =2^R -1;
  xi = Q*epidot/snr/(al1 - al2*epidot);
  pa(ki) = 1/(factorial(Q-1))^P * xi^(P*(Q+1)/2) *(xi^(- (Q+1)/2)*factorial(Q-1)...
      -2*xi^(-1/2)*besselk(Q, 2*xi^(1/2)))^P;
  
  %approximated results
  if Q>1
      paa(ki) = xi^P/(Q-1)^P;
  elseif Q==1
      paa(ki) =  xi^N*(-log(xi))^N;
  end
  
end
paa(1:3) = 0;
semilogy(snrdb, p,'-d',snrdb,pa,'-*',snrdb,min(1,paa),'-x',snrdb,po,'-s',snrdb,pdft,'-o')
%semilogy(snrdb, p,snrdb,po,snrdb,pdft) 
  