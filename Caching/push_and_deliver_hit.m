clear all%function [p1, p2] = fund_mac(K,ct)
fc = 2e9; %Carrier frequency (Hz)
c= 3e8; %speed of light (m/s)
wavelength = c/fc; % in m
eta= (wavelength/(4*pi))^2;
BW=10*10^6; %10 MHz 
Nf=10;%dB
sigma2_dbm= -180+10*log10(BW)+Nf; %Thermal noise in dBm
sigma_square=10^((sigma2_dbm-30)/10);

fl =  [3:-1:1];
gammaf = 1.5;
pf = 1./fl.^gammaf./sum(1./fl.^gammaf);

%load Data_A_11_2.mat;%A=polynomial(11,1);
 snrdb = [10: 10 : 40]; ct=5000;
 K=2; Rc = 50; al=3; R1=1; R2=8; rmax = 1000;   lamb=0.00001; thd = 100;
 R0=0.5; R1 = 3; R2=3.5;R3=8;
 beta0 = 4/8; beta1=3/8; beta2=2/8; beta3=1/8;  
 
 t = 5; %threshold
 F = 3; %number of files
 beta1 = 3/4; beta2 = 1/4;
 N=20; %Gassus variables
 
  
 for k =  1: length(snrdb)
    snr = 10^((snrdb(k)-30)/10)*eta/sigma_square*10;     
    ep0=(2^(R0/4)-1); zeta0 =    beta0 - ep0* ( beta1+beta2+beta3); 
    ep1=(2^(R1/4)-1); zeta1 =    beta1 - ep1* ( beta2+beta3); 
    ep2=(2^(R2/4)-1); zeta2 =    beta2 - ep2* ( beta3); 
    ep3=(2^(R3/4)-1); zeta3 =    beta3; 
    tau0 = (1/min( [snr*zeta0/ep0   ]  ))^(1/al);
    tau1 = (1/min( [snr*zeta0/ep0 snr*zeta1/ep1   ]  ))^(1/al);
    tau2 = (1/min( [snr*zeta0/ep0 snr*zeta1/ep1 snr*zeta2/ep2    ]  ))^(1/al);
    tau3 = (1/min( [snr*zeta0/ep0 snr*zeta1/ep1 snr*zeta2/ep2 snr*zeta3/ep3  ]  ))^(1/al);
   sum1=0; sum21=0; sum41=0;sum22=0; sum42=0;sum23=0; sum43=0;sum20=0; sum40=0;
   sum3=0;  sum5=0; sum6=0;  
  for ix = 1 :ct       
      %get the locations for the base stations
     bsn  = poissrnd(pi*rmax^2*lamb);
     pppind = 1;pppind2 = 1; pppind3 = 1;
     cx1=zeros(bsn,1); cy1=zeros(bsn,1); dx=zeros(bsn,1); cx2=0; cy2=0;cx3=0; cy3=0;
     while pppind<=bsn                 
         cx1(pppind,1) = sign(randn(1,1))*rand(1,1)*rmax;
         cy1(pppind,1) = sign(randn(1,1))*rand(1,1)*rmax;
         dx(pppind,1) = (cx1(pppind))^2+(cy1(pppind))^2;
         if dx(pppind,1)<rmax^2
             pppind = pppind+1;
             hx = abs(complex(sqrt(0.5)*randn(1,1),sqrt(0.5)*randn(1,1)))^2;
             zx = hx/(sqrt(dx(pppind-1,1))^al);
             if zx>tau3^al %only the good ones will be included
                 pppind2 = pppind2 +1;
                 cx2(pppind2-1,1) = cx1(pppind-1,1);
                 cy2(pppind2-1,1) = cy1(pppind-1,1);
             end
             if zx>(2^R3-1)/snr %only the good ones will be included
                 pppind3 = pppind3 +1;
                 cx3(pppind3-1,1) = cx1(pppind-1,1);
                 cy3(pppind3-1,1) = cy1(pppind-1,1);
             end
         end
     end

      %the location of the newcomer 
      cx=rmax/2; cy=rmax/2;
%      pppind4=1;
%      while pppind4<=1                 
%          cx(pppind4,1) = sign(randn(1,1))*rand(1,1)*rmax;
%          cy(pppind4,1) = sign(randn(1,1))*rand(1,1)*rmax;
%          dx3(pppind4,1) = (cx(pppind4))^2+(cy(pppind4))^2;
%          if dx3<rmax^2
%              pppind4 = pppind4+1;
%          end
%      end
     
     if pppind2==1 %no sucessful nodes for NOMA
         sum1 = sum1+1;  
     else
         distancex = sqrt( (cx2-cx).^2 + (cy2-cy).^2);
         if min(distancex)>thd
             sum1 = sum1+1;
         end         
     end
     if pppind3==1 %no sucessful nodes for OMA
         sum20 = sum20+1;  
     else
         distancex2 = sqrt( (cx3-cx).^2 + (cy3-cy).^2);
         if min(distancex2)>thd
             sum20 = sum20+1;
         end         
     end

  end
  
    p1(k) = sum1/ct; %user 2 to decode f3
    p2(k) = sum20/ct; %user 2 to decode f3
    
    %theorethic 
    nx=[1:1:N];
    wl = cos((2*nx-1)/2/N*pi);
    rm = sqrt((cx)^2+(cy)^2);
    ba = 2*lamb*thd * sum( pi/N*sqrt(1-wl.^2).*...
        (exp(-tau3^al*(rm+thd*wl).^al)) .* (rm+thd*wl) .* ...
        acos( ((rm+thd*wl).^2 +rm^2 - thd^2)/2/rm./(rm+thd*wl) )    ); 
    pa(k) =   exp(-ba);
    
  end
semilogy(snrdb,p2,snrdb,p1,snrdb,pa)