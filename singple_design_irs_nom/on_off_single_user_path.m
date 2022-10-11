clear all 
 snrdb = [0: 2 : 30]; ct=5000;
  N=64;
 al1 = 0.8; al2=0.2;
 R = 2;
 P=N; Q=N/P;
 V=kron(eye(P), ones(Q,1))/sqrt(Q);
 dsr=20;drd=10;
 al=3.5;
 
 sigma2_dbm= -70;%+10*log10(BW)+Nf; %Thermal noise in dBm
sigma_square=10^((sigma2_dbm-30)/10);
 
for ki =  1:length(snrdb)
     snr = 10^((snrdb(ki)-30)/10)/sigma_square;
        sum0=0; sum1=0; sum2=0; sum3=0;sum4=0;   sum5=0;
  for ix = 1 :ct            
     %%%%%%%%%%%%%%%%%%%%%%%%% High mobility user
     h1 = complex(sqrt(0.5)*randn(N,1),sqrt(0.5)*randn(N,1));
     h2 = complex(sqrt(0.5)*randn(N,1),sqrt(0.5)*randn(N,1));
      
     D1 = diag( complex(sqrt(0.5)*randn(N,1),sqrt(0.5)*randn(N,1)) );
     
      for fi =1 : P
          fvector = V(:,fi);
          SINR(fi) = abs(fvector'*D1*h1)^2*al1/dsr^al/drd^al/(abs(fvector'*D1*h1)^2*al2/dsr^al/drd^al  + 1/snr);
      end
      if max(SINR)<2^R -1
          sum1 = sum1 +1;
      end
      %optimal
      fo = D1*h1/sqrt(h1'*D1'*D1*h1);
      SINRo = abs(fo'*D1*h1)^2*al1/dsr^al/drd^al/(abs(fo'*D1*h1)^2*al2/dsr^al/drd^al  + 1/snr);
      if  (SINRo)<2^R -1
          sum2 = sum2 +1;
      end
      
      %relaying
      if log2(1+abs(h1(1))^2/dsr^al*snr)<4*R
          sum3 = sum3 +1;
      elseif log2(1+abs(h2(1))^2/drd^al*snr)<4*R
          sum3 = sum3 +1;
      end
      
  end
  
  p(ki) = sum1/ct;
  po(ki) = sum2/ct;
  pc(ki) = sum3/ct;
   
   
end
paa(1:3) = 0;
semilogy(snrdb,pc,snrdb, p,'-d',snrdb,po,'-*')
%semilogy(snrdb, p,snrdb,po,snrdb,pdft) 
  