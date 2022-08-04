clear all%function [p1, p2] = fund_mac(K,ct)
            
N=6;  K=4;   snrdb=[10 : 5 : 30]; Pu=2;  

load coorelation_data.mat;
%%%%%%%%%%% get the correlation matrix
R = Rr(:,:,1);
r=rank(R); Mt = M-r*(K-1);
[a1,b1 ] = eig(R); a11=a1'; U = a11(M-r+1:M,:); Lam = (b1(M-r+1:M,M-r+1:M))^(1/2) ;
% AX =   complex(sqrt(0.5)*randn(M,Mt),sqrt(0.5)*randn(M,Mt)); [ax,bx,cx]=svd(AX);
% P = ax(:,1:Mt);
A = Lam*U*P;
Ay = inv(P'*R*P); ai = Ay(1,1); ai=1/ai;
%%%%%%%%%%%%%%%%%%%%%% start the iteration
a1 = 3/4; a2 = 1/4;  R1=0.5; R2=5;   ct=5000;  RR =[R1 R2];
for k =  1:length(snrdb)
 sum1=0; sum2=0; sum3=0; sum4=0; sum5=0; sum6=0;
 snr = 10^(snrdb(k)/10);   eps1=(2^R1-1);eps2=(2^R2-1); 
 xix = [max(1e-20,a1 - eps1*(a2))/eps1  (a2)/eps2 ]*snr; taux = 0.5;%1/snr;

sum1 = 0;
 for i = 1 :ct       
     for j = 1 : Pu
         G = complex(sqrt(0.5)*randn(N,r),sqrt(0.5)*randn(N,r));   
         W1 = A'*G'*G*A;
         W2 = inv(W1);
         chy(j) = real(W2(1,1));
     end
     chy = sort(chy, 'descend');
     chx = chy(1,1);
     set2 = sum( (sign(1./chy-taux)+1)/2 ); %%% the size of SET2 
     set1 = Pu-set2;
     
    %%%%%%%%%%%%%%%%%%% with one bit
     if 1/chx > taux % user in SET 2
         ui = randsrc(1,1, [1 : Pu]);
         xi_ui = min( xix(1:ui));
         if 1/chx > 1./xi_ui
             sum1 = sum1+RR(ui);
         end
     else %%%% user in SET 1
         ui = randsrc(1,1, [1: set1]);
         xi_ui = min( xix(1:ui));
         if 1/chx > 1./xi_ui
             sum1 = sum1+RR(ui);
         end
     end

      %%%% with perfect CSI
         xi_ui2 = xix(1);
         if 1/chx > 1./xi_ui2
             sum2 = sum2+R1;
         end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% user 2     
    %%%%%%%%%%%%%%%%%%% with one bit
    chx = chy(2);
     if 1/chx > taux % user in SET 2
         ui = randsrc(1,1, [Pu-set2+1 : Pu]);
         xi_ui = min( xix(1:ui));
         if 1/chx > 1./xi_ui
             sum3 = sum3+RR(ui);              
         end
     else %%%% user in SET 1
         ui = randsrc(1,1, [1: set1]);
         xi_ui = min( xix(1:ui));
         if 1/chx > 1./xi_ui
             sum3 = sum3+RR(ui);
         end
        
     end

      %%%% with perfect CSI
         xi_ui2 = min(xix(1:2));
         if 1/chx > 1./xi_ui2
             sum4 = sum4+R2;
         end
    
         %%%%%%%%%%%%%%%%%% OMA
         
         if snr/chy(1)>(2^(Pu*R1)-1)
             sum5 = sum5 +R1;
         end
         if snr/chy(2)>(2^(Pu*R2)-1)
             sum6 = sum6 +R2;
         end
 end 
 p1(k) = sum1/ct;   p2(k) = sum2/ct;      p3(k) = sum3/ct;   p4(k) = sum4/ct; p5(k) = sum5/ct;   p6(k) = sum6/ct; 
 
 xi1 = min((a1 - eps1*(a2))/eps1);
 xi2 = min((a1 - eps1*(a2))/eps1,a2/eps2);
 %%% analytical
% 
%   pa2(k) = 2*max(0, Fp(1/taux,ai, N, Mt) - Fp(snr*xi2,ai, N, Mt) )*(1 - Fp(1/taux,ai, N, Mt)) +  Fp(1/taux,ai, N, Mt) *max(0, Fp(1/taux,ai, N, Mt) - Fp(snr*xi1,ai, N, Mt) ) + Fp(1/taux,ai, N, Mt) *max(0, Fp(1/taux,ai, N, Mt) - Fp(snr*xi2,ai, N, Mt) ) ...
%      -   0.5*max(0, (Fp(1/taux,ai, N, Mt))^2 - (Fp(snr*xi1,ai, N, Mt))^2 ) -   0.5*max(0, (Fp(1/taux,ai, N, Mt))^2 - (Fp(snr*xi2,ai, N, Mt))^2 ) ...
%      +0.5*(1-Fp(max(1/taux, snr*xi1),ai, N, Mt))^2 +0.5*(1-Fp(max(1/taux, snr*xi2),ai, N, Mt))^2 
%  pa1(k) = (1-2*(1 - Fp(max(1/taux, snr*xi1),ai, N, Mt))*Fp(1/taux,ai, N, Mt))*R1...
%      +(1-0.5*(1 - (Fp(max(1/taux, snr*xi1),ai, N, Mt))^2)  -  (1 - 0.5*Fp(max(1/taux, snr*xi1),ai, N, Mt))*Fp(1/taux,ai, N, Mt))*R1...
%      + (1-0.5*(1 - (Fp(max(1/taux, snr*xi2),ai, N, Mt))^2)- (1 - Fp(max(1/taux, snr*xi2),ai, N, Mt))*Fp(1/taux,ai, N, Mt))*R2 ...
%      +   (1- max(0,0.5*( Fp(1/taux,ai, N, Mt))^2 -0.5*( Fp(snr*xi1,ai, N, Mt))^2))*R1   + (1-max(0,0.5* (Fp(1/taux,ai, N, Mt))^2 - 0.5*(Fp(snr*xi2,ai, N, Mt))^2))*R2;
  
end
 plot(snrdb, p1+p3, snrdb, p2+p4,snrdb,p5+p6)
  