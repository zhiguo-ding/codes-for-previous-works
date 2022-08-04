clear all%function [p1, p2] = fund_mac(K,ct)
%close all
            
N=5;  K=4;   snrdb=[10 : 5 : 30];

%%%%%%%%%%% get the correlation matrix
M=50; Delta=15/360*2*pi;lamd=1/(2*10^10); thetag = -pi+Delta; 
D=0.5/sqrt( (1 - cos(2*pi/M))^2 +(sin(2*pi/M))^2);
stepd = Delta/100;
dd = [-Delta : stepd : Delta];
for m =1 : M
    B(1:2,m) = [lamd*D*cos((m-1)*2*pi/M) lamd*D*sin((m-1)*2*pi/M)];
end
for m = 1: M
    for n = 1 :M
        sum1 = 0;
        for p = 1 : length(dd)
            al = dd(p);
            kx = -2*pi/lamd*[cos(al+thetag) ; sin(al+thetag)];
            sum1 = sum1 + 1/2/Delta* exp(sqrt(-1)*kx'*(B(:,m) -B(:,n))  )*stepd;
        end
        R(m,n)=sum1;
    end
end
r=rank(R); Mt = M-r*(K-1);
[a1,b1 ] = eig(R); a11=a1'; U = a11(M-r+1:M,:); Lam = (b1(M-r+1:M,M-r+1:M))^(1/2) ;
AX =   complex(sqrt(0.5)*randn(M,Mt),sqrt(0.5)*randn(M,Mt)); [ax,bx,cx]=svd(AX);
P = ax(:,1:Mt);
A = Lam*U*P;
Ay = inv(P'*R*P); ai = Ay(1,1); ai=1/ai;
%%%%%%%%%%%%%%%%%%%%%% start the iteration
a1 = 3/4; a2 = 1-a1; R1=1.9; R2=1; ct=50000; 
for k =  1:length(snrdb)
 sum1=0; sum2=0; sum3=0; sum4=0; sum5=0; 
 snr = 10^(snrdb(k)/10);   eps1=(2^R1-1);

sum1 = 0;
 for i = 1 :ct             
     G = complex(sqrt(0.5)*randn(N,r),sqrt(0.5)*randn(N,r));   
     W1 = A'*G'*G*A;
     W2 = inv(W1);
     
     x = W2(1,1);
     
     sinr1 = snr*a1/x/ (1 + snr*a2/x);
     
     if real(sinr1)<eps1
         sum1 = sum1+1;
     end    
     
 end 
 p1(k) = sum1/ct;
 
 xix = (a1 - eps1*a2)/eps1;
 pa1(k) = gammainc( real(1/snr/ai/xix), N-Mt+1);
end
semilogy(snrdb,p1,snrdb,pa1)

 