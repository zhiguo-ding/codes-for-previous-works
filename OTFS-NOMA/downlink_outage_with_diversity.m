clear all 
 snrdb = [0: 5 : 30]; ct=50000;
 P0 = 4; Pi = P0;
 M=8;  N=8; 
 DD = [  0 1 ; 0 3;  1 5 ;  1 7 ];  %for M=8
 
 %M=16;  N=16;   
 %DD = [  0 2 ; 0 6;  2 10 ;  2 14 ]; %for M=16

  R0 = 1.9; %snr0 = (2^R0-1)/h0;
  Ri= 0;
  
  gamma0=3/4; gamma1=1-gamma0;
      
 for ki =  1:length(snrdb)
    snr = 10^((snrdb(ki))/10);  %snr0 = snr;
    eps0 = (2^R0-1);  
    epsi = (2^Ri-1);  
        
   sum0=0; sum1=0; sum2=0; sum3=0;sum4=0;   sum5=0; sum6=0; sum50=0;
  for ix = 1 :ct            
     %%%%%%%%%%%%%%%%%%%%%%%%% High mobility user
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
         
          hx=d(:);%sqrt(N*M)*kron(dftmtx(N)/sqrt(N), dftmtx(M)/sqrt(M))*H(:,1);
     
     %OMA -ZF
    SINR0 = snr/(  sum(sum(1./abs(d).^2))/N/M);
    if P0*SINR0 <eps0
        sum0 = sum0 +1;
    end
    

    %%%%%% FD-ZF
    SINR1 = snr*gamma0/(snr*gamma1 + sum(sum(1./abs(d).^2))/N/M);
    
    if P0*SINR1 <eps0
        sum1 = sum1 +1;
    end
    
    %%%%%% OMA FD-DFE random
    SINR20 = snr*sum(abs(hx(1)).^2);
    
    if SINR20 <eps0
        sum5 = sum5 +1;
    end      
    %%%%%% OMA FD-DFE
    SINR20 = snr*sum(abs(hx(1:P0)).^2);
    
    if SINR20 <eps0
        sum50 = sum50 +1;
    end   
    %%%%%% FD-DFE no grouping
    SINR2 = snr*gamma0/(snr*gamma1 +1/abs(d(1,1))^2);
    
    if SINR2 <eps0
        sum2 = sum2 +1;
    end      
    %%%%%% FD-DFE grouping
    SINR2 = snr*gamma0/(snr*gamma1 +1/sum(abs(hx(1:P0)).^2));
    
    if SINR2 <eps0
        sum6 = sum6 +1;
    end
    
    
  end
  p0(ki) = sum0/ct; %oma zf
  p01(ki) = sum5/ct; % oma DFE
    p010(ki) = sum50/ct; % oma DFE with fading

  p1(ki) = sum1/ct; %ZF
    p2(ki) = sum2/ct;%DFE wot grouping
    p3(ki) = sum6/ct; %DFE with grouping
    
    p4(ki) = 1-sum3/ct; % NOMA user ZF
    

  
 end
 semilogy(snrdb, p0, snrdb,p01,snrdb,p010,snrdb,p1,snrdb,p2,snrdb,p3)