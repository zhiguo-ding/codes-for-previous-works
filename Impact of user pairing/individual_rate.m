clear all%function [p1, p2] = fund_mac(K,ct)
%close all
R=0; Rx=10;al=4; ct=5000;              snrdb = [5: 5 : 40];     M=5;
           
m=2; n=M;
for k =  1:length(snrdb)
 sum1=0; sum2=0; sum3=0; sum4=0; sum5=0; 
 snr = 10^(snrdb(k)/10);
 for i = 1 :ct
             
     h = complex(sqrt(0.5)*randn(M,1),sqrt(0.5)*randn(M,1));       
     h = abs(h).^2;
     h=sort(h,'ascend'); h1=h(m); h2=h(n);

     %%%%%%%%%%%%%%%%%%%% random use
     Rdx1(i) = log2(1+snr*h1)/2;
     Rdx2(i) = log2(1+snr*h2)/2;
     %%%%%%%%%%%%%%%%%%%%%%% end
    
     a1=4/5;%h2/(h1+h2);
     a2=1-a1;%0.9;%h1/(h1+h2);
  
     Rx1(i) = log2( 1+snr*h1*a1/(snr*h1*a2+1) );     
     Rx2(i) = log2(1+snr*h2*a2);
     
     if Rx1(i)-Rdx1(i)<R
         sum1 = sum1+1;
     end
     
     if Rx2(i)-Rdx2(i)<R
         sum2 = sum2+1;
     end
     
     if h1>(1-2*a2)/snr/a2^2
         sum3 = sum3 +1;
     end
 end
    p1(k) = sum1/ct;
    p2(k) = sum2/ct;
    p3(k) = sum3/ct;
    
    %analytical resultw
    w5=factorial(M)/factorial(m-1)/factorial(M-m);
    sum1 = 0;
    for i = 0 : m-1
        sum1 = sum1+ factorial(m-1)/factorial(m-1-i)/factorial(i) *(-1)^i*w5/(M-m+i+1)...
            *( 1 - exp(-(1-2*a2)*(M-m+i+1)/snr/a2^2 )  );
    end
    p11(k) = 1-sum1;
    p11a(k) = 1- w5*(1-2*a2)^m/m/snr^m/a2^(2*m);

    w3=factorial(M)/factorial(n-1)/factorial(M-n);
    sum1 = 0;
    for i = 0 : n-1
        sum1 = sum1+ factorial(n-1)/factorial(n-1-i)/factorial(i) *(-1)^i*w3/(M-n+i+1)...
            *( 1 - exp(-(1-2*a2)*(M-n+i+1)/snr/a2^2 )  );
    end
    p22(k) = sum1;
    p22a(k) =  w3*(1-2*a2)^n/n/snr^n/a2^(2*n);

end
%p11a(1:3) =0;
%semilogy(snrdb,p1,'-*',snrdb,p2,snrdb,p11,snrdb,p11a,'-s',snrdb,p22,snrdb,p22a,'-s')
semilogy(snrdb,p1,'-*',snrdb,p2,'-o',snrdb,p22,snrdb,p22a,'-s')

% figure
% plot(snrdb,c1,snrdb,c2)