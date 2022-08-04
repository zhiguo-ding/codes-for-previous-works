clear all%function [p1, p2] = fund_mac(K,ct)
 
            
N=6;  K=4;   snrdb=[10 : 5 : 30]; Pu=3;  

load coorelation_data.mat;
%%%%%%%%%%% get the correlation matrix
% M=50; Delta=15/360*2*pi;lamd=1/(2*10^10); thetag = -pi+Delta; 
% D=0.5/sqrt( (1 - cos(2*pi/M))^2 +(sin(2*pi/M))^2);
% stepd = Delta/100;
% dd = [-Delta : stepd : Delta];
% for m =1 : M
%     B(1:2,m) = [lamd*D*cos((m-1)*2*pi/M) lamd*D*sin((m-1)*2*pi/M)];
% end
% for m = 1: M
%     for n = 1 :M
%         sum1 = 0;
%         for p = 1 : length(dd)
%             al = dd(p);
%             kx = -2*pi/lamd*[cos(al+thetag) ; sin(al+thetag)];
%             sum1 = sum1 + 1/2/Delta* exp(sqrt(-1)*kx'*(B(:,m) -B(:,n))  )*stepd;
%         end
%         R(m,n)=sum1;
%     end
% end
R = Rr(:,:,1);
r=rank(R); Mt = M-r*(K-1);
[a1,b1 ] = eig(R); a11=a1'; U = a11(M-r+1:M,:); Lam = (b1(M-r+1:M,M-r+1:M))^(1/2) ;
AX =   complex(sqrt(0.5)*randn(M,Mt),sqrt(0.5)*randn(M,Mt)); [ax,bx,cx]=svd(AX);
P = ax(:,1:Mt);
A = Lam*U*P;
Ay = inv(P'*R*P); ai = Ay(1,1); ai=1/ai;
%%%%%%%%%%%%%%%%%%%%%% start the iteration
a1 = 5/8; a2 = 2/8; a3 = 1-a1-a2; R1=1.4; R2=1.4; R3 =2; ct=50000; 
for k =  1:length(snrdb)
 sum1=0; sum2=0; sum3=0; sum4=0; sum5=0; 
 snr = 10^(snrdb(k)/10);   eps1=(2^R1-1);eps2=(2^R2-1);eps3=(2^R3-1);

sum1 = 0;
 for i = 1 :ct       
     for j = 1 : Pu
         G = complex(sqrt(0.5)*randn(N,r),sqrt(0.5)*randn(N,r));   
         W1 = A'*G'*G*A;
         W2 = inv(W1);
         chx(j) = W2(1,1);
     end
     chx = sort(chx,'descend');
     
     %%%%% user 1    
     sinr1 = snr*a1/chx(1)/ (1 + snr*a2/chx(1)+snr*a3/chx(1));
     
     if real(sinr1)<eps1
         sum1 = sum1+1;
     end    
     %%%%% user 2    
     sinr1 = snr*a1/chx(2)/ (1 + snr*a2/chx(2)+snr*a3/chx(2));
     sinr2 = snr*a2/chx(2)/ (1 + snr*a3/chx(2));
     
     if real(sinr1)<eps1 | real(sinr2)<eps2
         sum2 = sum2+1;
     end    
     %%%%% user 3    
     sinr1 = snr*a1/chx(3)/ (1 + snr*a2/chx(3)+snr*a3/chx(3));
     sinr2 = snr*a2/chx(3)/ (1 + snr*a3/chx(3));
     sinr3 = snr*a3/chx(3);
     
     if real(sinr1)<eps1 | real(sinr2)<eps2 | real(sinr3)<eps3
         sum3 = sum3+1;
     end    
%      if snr/chx(3) <max([(a1 - eps1*(a2+a3))/eps1 (a2 - eps1*(a3))/eps2 (a3)/eps3 ] )
%          sum3 = sum3+1;
%      end
     
 end 
 p1(k) = sum1/ct;
  p2(k) = sum2/ct;
   p3(k) = sum3/ct;
% % pa1(k) = gammainc( real(1/snr/ai/xix), N-Mt+1);
 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%% analytical 
 
     pa1(k) =0;
     p=1;xix = (a1 - eps1*(a2+a3))/eps1;
     for l = 1: p-1+1
        i=l-1;
        pip = factorial(Pu)/factorial(Pu-p)/factorial(p-1);
        Fp = 1 - gammainc( real(1/snr/ai/xix), N-Mt+1);
        pa1(k) = pa1(k) + factorial(p-1)/factorial(i)/factorial(p-1-i)*(-1)^i*pip/(Pu-p+i+1)*(1- (Fp)^(Pu-p+i+1) );
     end
     %%%%%%%%%%%%%%%%%%%user 2
   pa2(k) =0;
    p=2;xix = min((a1 - eps1*(a2+a3))/eps1,(a2 - eps1*(a3))/eps2);
     for l = 1: p-1+1
        i=l-1;
        pip = factorial(Pu)/factorial(Pu-p)/factorial(p-1);
        Fp = 1 - gammainc( real(1/snr/ai/xix), N-Mt+1);
        pa2(k) = pa2(k) + factorial(p-1)/factorial(i)/factorial(p-1-i)*(-1)^i*pip/(Pu-p+i+1)*(1- (Fp)^(Pu-p+i+1) );
     end
    %%%%%%%%%%%%%%%%%%%user 3
   pa3(k) =0;
    p=3;xix = min([(a1 - eps1*(a2+a3))/eps1 (a2 - eps1*(a3))/eps2 (a3)/eps3 ] );
     for l = 1: p-1+1
        i=l-1;
        pip = factorial(Pu)/factorial(Pu-p)/factorial(p-1);
        Fp = 1 - gammainc( real(1/snr/ai/xix), N-Mt+1);
        pa3(k) = pa3(k) + factorial(p-1)/factorial(i)/factorial(p-1-i)*(-1)^i*pip/(Pu-p+i+1)*(1- (Fp)^(Pu-p+i+1) );
     end
end
semilogy(snrdb,p1,snrdb,p2,snrdb,p3 )
 
semilogy(snrdb,pa1,'-d',snrdb,pa2,'-o',snrdb,pa3,'-x' )
 