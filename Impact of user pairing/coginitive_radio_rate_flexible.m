clear all%function [p1, p2] = fund_mac(K,ct)
%close all
R=0.9; Rx=10;al=4; ct=50000;              snrdb = [5: 5 : 30];     M=5;
           I=5;
m=4; n=M;
for k =  1:length(snrdb)
 sum1=0; sum2=0; sum3=0; sum4=0; sum5=0; 
 snr = 10^(snrdb(k)/10);     a=1+I; b=I/snr;eps1=(2^R-1)/snr;

 for i = 1 :ct
             
     h = complex(sqrt(0.5)*randn(M,1),sqrt(0.5)*randn(M,1));       
     h = abs(h).^2;
     h=sort(h,'ascend'); h1=h(m); h2=h(n);

     %%%%%%%%%%%%%%%%%%%% random use
     Rdx1(i) = log2(1+snr*h1);
     %%%%%%%%%%%%%%%%%%%%%%% end
    
     a2 = max(0, (h1-I/snr)/h1/(1+I));
  
     Rx2(i) = log2(1+snr*h2*a2);     

     
 end
    c1(k) = mean(Rx2);
 
end

plot(snrdb,c1 )%,snrdb,c3)
 