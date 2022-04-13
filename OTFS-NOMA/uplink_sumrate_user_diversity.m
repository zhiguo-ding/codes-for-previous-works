clear all 
 snrdb = [0: 5 : 30]; ct=5000;
 P0 = 4; Pi = P0; K=16;
 %M=8;  N=8; 
 %DD = [  0 1 ; 0 3;  1 5 ;  1 7 ];  %for M=8
 
 M=16;  N=16;   
 DD = [  0 1 ; 0 3;  1 5 ;  1 7 ]; %for M=16

  R0 = 0.5; %snr0 = (2^R0-1)/h0;
  R0t = R0*ones(1,N*M);%[1:1:N*M]*N*M*R0/(N*M/2*(1+N*M));
  Ri= 1;
  
      
 for ki =  1:length(snrdb)
    snr = 10^((snrdb(ki))/10);  %snr0 = snr;
    eps0 = (2^R0-1);  
    epsi = (2^Ri-1);  
        
   sum0=0; sum1=0; sum2=0; sum3=zeros(N*M,1); sum4=0;   sum5=zeros(N*M,1);
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
    
     %OMA -ZF
    SINR0 = snr/(  sum(sum(1./abs(d).^2))/N/M);
    if SINR0 <eps0
        sum0 = sum0 +1;
    end
    %%%%%% OMA FD-DFE
    [L,p] = chol(H'*H);
    dcl = diag(L);

    dvec = abs(dcl(:)).^2;
    dvec = sort(dvec,'ascend');    
    SINR20 = snr*dvec;
    for ni =1 : N*M    
        if SINR20(ni) <2^R0t(ni)-1
            sum5(ni) = sum5(ni) +1;
        end     
    end
        
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
    if SINR_NOMA<epsi
        sum1 = sum1 +1;
    end
    
    %%%%%%%%%%%%%%NOMA SECOND STAGE
    for m =1:M %M NOMA users
        SINR_X(m)= snr*di(m)/(snr*max(abs(d(m,:)).^2) +1);
    end
    
        
    %%%%%% FD-ZF
    SINR0 = snr/(  sum(sum(1./abs(d).^2))/N/M);
    
    if min(SINR_X)>epsi & SINR0 >eps0 %the weakest NOMA can do it and
        sum2 = sum2 +1;
    end    

    %%%%%% FD-DFE
     SINR20 = snr*dvec;
    for ni =1 : N*M    
        if min(SINR_X)>epsi & SINR20(ni) >2^R0t(ni)-1
            sum3(ni) = sum3(ni) +1;
        end     
    end
    
%     SINR20 = snr/( 1/abs(d(1,1))^2);
%     
%     if min(SINR_X)>epsi & SINR20 >eps0
%         sum3 = sum3 +1;
%     end      
    
      
  end
  p0(ki) = sum0/ct; %oma zf
  p20(1:N*M,ki) = sum5/ct;%oma DFE
  %p01(ki) = sum5/ct; % oma DFE
  p1(ki) = sum1/ct; % THE FIRST stage
  
  
    p2(ki) = 1- sum2/ct;%ZF
    %p3(ki) = 1- sum3/ct;%dfe
    p3(1:N*M,ki) = 1 - sum3/ct;%DFE

 
  
 end
  oma_zf = (1-p0)*R0;
 oma_dfe = ((1-p20)'*R0t')/N/M;
 noma_zf=(1-p1)*Ri +(1-p2)*R0 ;%(1-p1)*R0+(1-p3)*Ri;
 noma_dfe=(1-p1')*Ri +    ((1-p3)'*R0t')/N/M ;

 %plot(snrdb, oma_zf, snrdb,oma_dfe, snrdb,noma_zf, snrdb, noma_dfe)
 
