clear all 
 snrdb = [0: 5 : 30]; ct=50;
 M=4; N=12;
 al1 = 0.8; al2=0.2;
 R = 2;
 P=12; Q=N/P;
 V=kron(eye(P), ones(Q,1))/sqrt(Q);
 
for ki =  1:length(snrdb)
     snr = 10^((snrdb(ki))/10);
        sum0=0; sum1=0; sum2=0; sum3=0;sum4=0;   sum5=0; 
  %analytical result
  epidot =2^R -1;
  xi = Q*epidot/snr/(al1 - al2*epidot);
  pa(ki) = 1/(factorial(Q-1))^P * xi^(P*(Q+1)/2) *(xi^(- (Q+1)/2)*factorial(Q-1)...
      -2*xi^(-1/2)*besselk(Q, 2*xi^(1/2)))^P;
   
  
end
paa(1:3) = 0;
semilogy( snrdb,pa,'-*' )
%semilogy(snrdb, p,snrdb,po,snrdb,pdft) 
  