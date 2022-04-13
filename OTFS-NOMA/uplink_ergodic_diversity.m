clear all 
 snrdb = [0: 5 : 30]; ct=10000;
 P0 = 4; Pi = P0; K=8;
 %M=8;  N=8; 
 %DD = [  0 1 ; 0 3;  1 5 ;  1 7 ];  %for M=8
 
 M=16;  N=16;   
 DD = [  0 2 ; 0 6;  2 10 ;  2 14 ]; %for M=16

  R0 = 1; %snr0 = (2^R0-1)/h0;
  Ri= 0.8;
  
  gamma0=3/4; gamma1=1-gamma0;
      
 for ki =  1:length(snrdb)
    snr = 10^((snrdb(ki))/10);  %snr0 = snr;
    eps0 = (2^R0-1);  
    epsi = (2^Ri-1);  
        
   sum0=0; sum1=0; sum2=0; sum3=0;sum4=0;   sum5=0;
  for ix = 1 :ct            

     %%%%%%%%%%%%%%%%%%%%%%%%% High mobility user channel
     h0 = complex(sqrt(0.5)*randn(P0,1),sqrt(0.5)*randn(P0,1))/sqrt(P0);
     
     H0 = zeros(N*M,N*M);
     
     %for [nu tau] - last line of A0: zeros(1, M-1-tau) 1 zeros(1,tau)
     % so the matrix will be 
     for i =1 : 1 : N
         A(:,:,i)  = zeros(M,M);
     end

     for i =1 : 1 : P0
         nu = DD(i,1);
         tau = DD(i,2);
         % nu+1 is the index for A from the right-hand            
         % N- (nu+1) +1 is the index for A from the left-hand
         
         %for [nu tau] - last line of A0: zeros(1, M-1-tau) 1 zeros(1,tau)
         % so the matrix will be 
         A(:,:,  (nu+1)   )  = A(:,:,  (nu+1)   )+...
                          h0(i)*[ zeros(tau,M-tau)     eye(tau); ...
                                  eye(M-tau)           zeros(M-tau,tau)];
     end
        
     H = zeros(N*M,N*M);
     for i = 1 :N
         for m = 1 : N
             H((i-1)*M+1:i*M,(m-1)*M+1:m*M ) = A(:,:,mod(m-i,N)+1);
         end
     end     
        
     A0 = vec2mat(H(:,1), M); A0 = A0.';
     d=sqrt(N*M) *dftmtx(M)/sqrt(M) * A0 * dftmtx(N)/sqrt(N);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOMA users    channel
    
     hi = complex(sqrt(0.5)*randn(Pi,K),sqrt(0.5)*randn(Pi,K))/sqrt(Pi);

     DFT = dftmtx(M);
     dix=  DFT(:,1:Pi)*hi;   
     dix_ab= abs(dix).^2;
     if K==1
         di = dix_ab;
     else
         di = max((dix_ab)'); %the M best users on M subcarriers
     end
        

    %%%%%%%%%%%%%%NOMA FIRST STAGE
    SINR_NOMA= snr*di(1)/(snr*abs(d(1,1))^2 +1);
    ratei(ix) = log2(1+ SINR_NOMA);
    
      
  end
   ri(ki) = mean(ratei); % THE FIRST stage  
 end
 
 plot(snrdb,ri)
   