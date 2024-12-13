clear all%function [p1, p2] = fund_mac(K,ct)
%close all
R=0.01; Rx=10;al=4; ct=5000;              snrdb = [5: 5 : 20];     N=20;
           
for k =  1:length(snrdb)
 sum1=0; sum2=0; sum3=0; sum4=0; sum5=0;ct2=0;
 snr = 10^(snrdb(k)/10);eps = (2^(2*R)-1)/snr;ep=(2^(2*R)-1); ep0 = (eps+sqrt(eps^2+eps))/2;
 for i = 1 :ct
             
     h = complex(sqrt(0.5)*randn(N,1),sqrt(0.5)*randn(N,1));       
     h = abs(h).^2;
     h=sort(h,'descend'); h1=h(N); h2=h(1);

     %%%%%%%%%%%%%%%%%%%% random use
     Rdx(i) = log2(1+snr*h(1))/2+log2(1+snr*h(N))/2;
     %%%%%%%%%%%%%%%%%%%%%%% end

     
     if h1>h2
         xx1 = h2;
         h2 = h1;
         h1 = xx1;
     end
     
     a1=4/5;%h2/(h1+h2);
     a2=1-a1;%0.9;%h1/(h1+h2);
  
     Rx1 = log2( 1+snr*h1*a1/(snr*h1*a2+1) );     
     Rx2 = log2(1+snr*h2*a2);
     
     Rxn(i) = Rx1+Rx2;
     
     x1(i) = log2(1+snr*h1)/2;
     x2(i) = log2(1+snr*h2)/2;
     
     y1(i) = Rx1;
     y2(i) = Rx2;
 end
     c1(k) = mean(Rdx);
     c2(k) = mean(Rxn);
     d1(k) = mean(x1);
     d2(k) = mean(x2);
     e1(k) = mean(y1);
     e2(k) = mean(y2);
     
end
 plot(snrdb,c1,snrdb,c2,'-o',snrdb,d1,snrdb,d2,snrdb,e1,'-x',snrdb,e2,'-s')
