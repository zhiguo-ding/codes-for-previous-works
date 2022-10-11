clear all 
 snrdb = [0: 5 : 30]; ct=2000;
 M=2; N=20;
 al1 = 0.8; al2=0.2;
 R = 1;
 
for ki =  1:length(snrdb)
     snr = 10^((snrdb(ki))/10);
        sum0=0; sum1=0; sum2=0; sum3=0;sum4=0;   sum5=0;
  for ix = 1 :ct            
     %%%%%%%%%%%%%%%%%%%%%%%%% High mobility user
     h1 = complex(sqrt(0.5)*randn(N,1),sqrt(0.5)*randn(N,1));
     h2 = complex(sqrt(0.5)*randn(N,1),sqrt(0.5)*randn(N,1));
      
     D1 = diag( complex(sqrt(0.5)*randn(N,1),sqrt(0.5)*randn(N,1)) );

     d= dftmtx(N)/sqrt(N);
     
      for fi =1 : N
          fvector = d(:,fi);
          SINR(fi) = abs(fvector'*D1*h1)^2*al1/(abs(fvector'*D1*h1)^2*al2 ...
              +abs(fvector'*D1*h2)^2 + 1/snr);
      end
      if max(SINR)<2^R -1
          sum1 = sum1 +1;
      end
     
  end
  
  p(ki) = sum1/ct;
  
end
semilogy(snrdb, p)
 
  