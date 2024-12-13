clear all%function [p1, p2] = fund_mac(K,ct)
%close all
R=0.01; Rx=3;al=2; ct=5000;         M=2;     snrdb = [0 : 5 : 50];  N=10;      ax=2*Rx; bx=ax; K=10;
           
nx = [1: N];
w = pi/N; thetan = cos((2*nx-1)*pi/2/N); b=-w*sqrt(1-thetan.^2).*(Rx*thetan/2+Rx/2); c=1+(Rx*thetan/2+Rx/2).^al; c0 = 0; b0=-sum(b);
A = matrix_kx(N);

for k =  1:length(snrdb)
 sum1=0; sum2=0; sum3=0; sum4=0; sum5=0;ct2=0;
 snr = 10^(snrdb(k)/10);eps = (2^(2*R)-1)/snr;ep=(2^(2*R)-1); ep0 = (eps+sqrt(eps^2+eps))/2;
 for i = 1 :ct
     h1 = complex(sqrt(0.5)*randn(1,1),sqrt(0.5)*randn(1,1));  
     h2 = complex(sqrt(0.5)*randn(1,1),sqrt(0.5)*randn(1,1));  
     h1 = abs(h1).^2;
     h2 = abs(h2).^2;
     h =  complex(sqrt(0.5)*randn(1,K),sqrt(0.5)*randn(1,K)); 
     h = abs(h).^2; h =sort(h,'ascend');
     h1 = h(1);
     h2 = h(K);
     
     %%%%%%%%%%%%%%%%%%%% random use
     Rdx = (log2(1+snr*h1)+log2(1+snr*h2))/2;
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
     Rxn = Rx1+Rx2;     
     dx(i) = Rxn - Rdx;     
     
     %%%%%%% cooperative
     g = complex(sqrt(0.5)*randn(1,1),sqrt(0.5)*randn(1,1));  g=abs(g)^2;
     Rx12 = min(log2( 1+snr*h1*a1/(snr*h1*a2+1) +snr*g ),  log2( 1+snr*h2*a1/(snr*h2*a2+1) )) ;     
     Rxn2 = Rx12+Rx2;     
     dx2(i) = Rxn2 - Rdx;     

 end
     c1(k) = mean(dx);
     c2(k) = mean(dx2);
     
end
 plot(snrdb,c1,'-d',snrdb,c2,'-o')
 %legend('Conventional, worst user', 'NOMA, worst user',  'Conventional1, sum rate','NOMA, sum rate','Conventoional2, sum rate')