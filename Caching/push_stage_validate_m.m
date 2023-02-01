clear all%function [p1, p2] = fund_mac(K,ct)
fc = 2e9; %Carrier frequency (Hz)
c= 3e8; %speed of light (m/s)
wavelength = c/fc; % in m
eta= (wavelength/(4*pi))^2;
BW=10*10^6; %10 MHz 
Nf=10;%dB
sigma2_dbm= -180+10*log10(BW)+Nf; %Thermal noise in dBm
sigma_square=10^((sigma2_dbm-30)/10);

 snrdb = [10: 5 : 40]; ct=50000;
 K=2; Rc = 50; al=3; R1 = 1; R2=1; R3=1;  rmax = 5500;    lamb=0.01/pi/Rc^2; 
 bsi = 1; %which base staiton
 t = 5; %threshold
 F = 3; %number of files
 beta1 = 3/4; beta2 = 1/4;
 N=50; %Gassus variables
 m=1;
  
 for k = 1:length(snrdb)
    snr = 10^((snrdb(k)-30)/10)*eta/sigma_square*10;  
    ep1 = 2^R1-1; ep2=2^R2-1; ep3=2^R3-1;  
     
    xi1 = beta1 - ep2*beta2; %file 2, so it is related to ep2
    xi2=beta2; %for file 3
    
  sum1=zeros(t,F); sum2=0; sum3=0; sum4=0; sum5=0; sum6=0;   
  
    %theorethic  

%      %%%%%%%%%  % user m to decode f3 m=2
    xc1 = snr/ep1; %the first part excluding Q1
    sumx = 0;
    for kl = 0: t-1
        sumx = sumx +  xc1^(2*kl/al)/factorial(kl)*(lamb*pi)^(kl);
    end
       ximin = min(xi1/ep2,xi2/ep3);%min(xi1,xi2); %f3
       tau2=(ep1/snr)^(-1/al);
       tau1 =   ( snr*ximin/(1+ep1+ep1*ximin) )^(1/al) ;
      Q1 = 0;
          for p = 0 : t-m-1
              sumx1 = 0;
              for l = 1: N
                  wl = cos( (2*l-1)/2/N*pi );
                  sumx1 = sumx1 + pi* (tau2-tau1)  /2/N*...
                      f_n(0.5* ( (tau2-tau1) *wl +  (tau2+tau1) ), lamb,t,m,p,al,ep1, ximin,snr)...
                      *sqrt(1-wl^2);
              end  
              Q1 = Q1 + (-1)^p*factorial(t-m-1)/factorial(p)/factorial(t-m-1-p)*sumx1;
          end
       Q1 = Q1*4*(lamb*pi)^(t)/factorial(t-m-1)/factorial(m-1);
      pam3(k) = Q1 + exp(-lamb*pi*xc1^(2/al))*sumx;
   end
semilogy( snrdb,pam3  )