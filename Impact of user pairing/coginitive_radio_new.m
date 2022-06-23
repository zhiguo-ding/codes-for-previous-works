clear all%function [p1, p2] = fund_mac(K,ct)
%close all
R=1; Rx=10;al=4; ct=1000000;              snrdb = [5: 5 : 30];     M=5;
           I=5;
m=2; n=M;
for k =  1:length(snrdb)
 sum1=0; sum2=0; sum3=0; sum4=0; sum5=0; 
 snr = 10^(snrdb(k)/10);     a=1+I; b=I/snr;eps1=(2^R-1)/snr;

 for i = 1 :ct
             
     h = complex(sqrt(0.5)*randn(M,1),sqrt(0.5)*randn(M,1));       
     h = abs(h).^2;
     h=sort(h,'ascend'); h1=h(m); h2=h(n);

     %%%%%%%%%%%%%%%%%%%% random use
     Rdx1(i) = log2(1+snr*h1);
     %Rdx2(i) = log2(1+snr*h2)/2;
     %%%%%%%%%%%%%%%%%%%%%%% end
    
     a2 = max(0, (h1-I/snr)/h1/(1+I));
  
     Rx2(i) = log2(1+snr*h2*a2);
     
     if Rx2(i)<R
         sum1 = sum1+1;
     end
     
%      if h1<b
%          sum3 = sum3+1;
%      else
%          if h2<a*eps1
%              if h2>b & h2<a*eps1 & h1>b
%                  sum3 = sum3+1;
%              end
%          else
%              if h2 > b+eps1 
%                  if h1>b & h1<b/(1 - a*eps1/h2)
%                     sum3 = sum3+1;
%                  end
%              else
%                  if h2>a*eps1 & h2<b+eps1 & h1>b
%                      sum3 = sum3+1;
%                  end       
%              end
%          end
%      end
     
%      b=I/snr; a=1+I;eps1=(2^R-1)/snr;
%      if h1>b & h1<b/(1-a*eps1/h2) & h2>b+eps1 %h1>b& h1<min( h2, b/(1- a*eps1/h2))
%          sum2 = sum2+1;
%      end
     
 end
    p1(k) = sum1/ct; 
    p2(k) = sum2/ct; 
    px(k) = sum3/ct;
    %analytical resultw
    w1=factorial(M)/factorial(m-1)/factorial(n-1-m)/factorial(M-n);
    w2=(1-2*a2)/a2^2;
    w3=factorial(M)/factorial(n-1)/factorial(M-n);
    w4=factorial(M)/factorial(m-1)/factorial(M-m);
    a=1+I; b=I/snr;eps1=(2^R-1)/snr;
        
     
    sum1 = 0;
    for i =0 : n-1-m 
        sum1 = sum1 + w1*factorial(n-m-1)/factorial(i)/factorial(n-1-m-i)*(-1)^i...
            *intg31(w2,n,m,i,snr,M,b,eps1,a);
    end
    
    for i =0 : n-1-m 
        sum1 = sum1 + w1*factorial(n-m-1)/factorial(i)/factorial(n-1-m-i)*(-1)^i...
            *intg3(w2,n,m,i,snr,M,b,eps1,a);
    end
    for i =0 : n-1-m 
        sum1 = sum1 + w1*factorial(n-m-1)/factorial(i)/factorial(n-1-m-i)*(-1)^i...
            *intg4(w2,n,m,i,snr,M,b,eps1,a);
    end
    
    for i =0:m-1
        sum1 = sum1 + w4*factorial(m-1)/factorial(i)/factorial(m-1-i)*(-1)^i/(M-m+i+1)...
            *(1 - exp(- (M-m+i+1)*b));
    end

    p3(k) = sum1;
    
    p4(k) = min(1,100/snr^((m)));
end
semilogy(snrdb,p1,'-*',snrdb,p3,'-s',snrdb,p4)
