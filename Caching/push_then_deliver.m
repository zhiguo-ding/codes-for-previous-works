clear all%function [p1, p2] = fund_mac(K,ct)
fc = 2e9; %Carrier frequency (Hz)
c= 3e8; %speed of light (m/s)
wavelength = c/fc; % in m
eta= (wavelength/(4*pi))^2;
BW=10*10^6; %10 MHz 
Nf=10;%dB
sigma2_dbm= -180+10*log10(BW)+Nf; %Thermal noise in dBm
sigma_square=10^((sigma2_dbm-30)/10);

 snrdb = [-10: 10 : 20]; ct=10000;
 K=2; Rc = 100; Rs =50; al=4; R1 = 1; R2=6; rmax = 40000;    lamb=0.01/pi/Rc^2; 
 N=20; %Gassus variables
 alpha1 = 3/4; alpha2 = 1/4;
 for k = 1:length(snrdb)
    snr = 10^((snrdb(k)-30)/10)*eta/sigma_square*10;  ep1 = 2^R1-1;   ep2 = 2^R2-1;
  sum1=zeros(K,1); sum2=zeros(K,1); sum3=zeros(K,1); sum4=0; sum5=0; sum6=0;  
  for ix = 1 :ct       
      %get the locations for the base stations
     bsn  = poissrnd(pi*rmax^2*lamb);
     pppind = 1;
     cx1=zeros(bsn,1); cy1=zeros(bsn,1); dx=zeros(bsn,1);
     while pppind<=bsn                 
         cx1(pppind,1) = sign(randn(1,1))*rand(1,1)*rmax;
         cy1(pppind,1) = sign(randn(1,1))*rand(1,1)*rmax;
         dx(pppind,1) = (cx1(pppind))^2+(cy1(pppind))^2;
         if dx(pppind,1)<rmax^2
             pppind = pppind+1;
         end
     end
%      %order these base stations to the comp user
%      [t1,t2] = sort(dx,'ascend');
%      cx1 = cx1(t2); cy1 = cy1(t2); %now all the base statins are ordered. 

     bsi = 1; %%% only focus on this particular cell

     %get the locations for those users      
     cx2=zeros(bsn,K); cy2=zeros(bsn,K); %each rwo contains K users  relative locations
 
     %the far user
     pppind = 1;         
         while pppind<=1                 
             cx2(1,pppind) = sign(randn(1,1))*rand(1,1)*Rc;
             cy2(1,pppind) = sign(randn(1,1))*rand(1,1)*Rc;
             dx2 = (cx2(1,pppind))^2+(cy2(1,pppind))^2;
             if dx2<Rc^2 & dx2>Rs^2
                 pppind = pppind+1;
             end
         end 
             cx3(bsi,1) = cx2(1)+cx1(bsi); cy3(bsi,1) = cy2(1)+cy1(bsi); %actual user locations

     %the near user
         pppind = 1;         
         while pppind<=1                 
             cx2(2,pppind) = sign(randn(1,1))*rand(1,1)*Rc;
             cy2(2,pppind) = sign(randn(1,1))*rand(1,1)*Rc;
             dx2 = (cx2(2,pppind))^2+(cy2(2,pppind))^2;
             if dx2<Rs^2
                 pppind = pppind+1;
             end
         end
         cx3(bsi,2) = cx2(2)+cx1(bsi); cy3(bsi,2) = cy2(2)+cy1(bsi); %actual user locations
         
         %the performance at BS_i  
         
     for kx = 1: K   
         Iinter = 0;
         for ii = 1 : bsn
             if ii ~= bsi
                 hij = abs(complex(sqrt(0.5)*randn(1,1),sqrt(0.5)*randn(1,1)))^2;
                 dij = sqrt( (cx3(bsi,kx)- cx1(ii)).^2 +(cy3(bsi,kx)- cy1(ii)).^2 ); % from BS_ii to user 1 in bsi
                 Iinter = Iinter + hij/dij^al;
             end
         end
         g = abs(complex(sqrt(0.5)*randn(1,1),sqrt(0.5)*randn(1,1)))^2;
         gd = sqrt( (cx3(bsi,kx)- cx1(bsi)).^2 +(cy3(bsi,kx)- cy1(bsi)).^2 );
         SINRi = g/gd^al * alpha1/( g/gd^al * alpha2 + Iinter +1/snr);   
         if SINRi<ep1
             sum1(kx) = sum1(kx) +1;
             sum2(kx) = sum2(kx) +1;
         else
             SINRi2 = g/gd^al * alpha2/(   Iinter +1/snr);   
             if SINRi2<ep2
                 sum2(kx) = sum2(kx) +1;
             end
         end 
         Rxx=[R1 R2];
         if g/gd^al  /(   Iinter +1/snr)<2^(2*Rxx(kx))-1
             sum3(kx) = sum3(kx) +1;
         end

     end
         %Iinter=0;

    end
    p1(k) = sum1(1)/ct;
    p2(k) = sum2(2)/ct;
    p3(k) = sum3(1)/ct;
    p4(k) = sum3(2)/ct;    
    
    %theorethic 
    nx=[1:1:N];
    thetanx = cos((2*nx-1)/2/N*pi);
    wnx = pi/2/N*sqrt(1-thetanx.^2).*(thetanx+1);  
    cnc = (Rc/2*thetanx+Rc/2).^al;
    cns = (Rs/2*thetanx+Rs/2).^al;

    muc = cnc*ep1/(alpha1 - ep1*alpha2);
    mus = cns*ep1/(alpha1 - ep1*alpha2);
    Liinterc = exp( -2*pi*lamb* (muc).^(2/al) /al * beta(2/al, (al-2)/al)  );
    Liinters = exp( -2*pi*lamb* (mus).^(2/al) /al * beta(2/al, (al-2)/al)  );
    pa(k) = Rc^2/(Rc^2-Rs^2)*sum(wnx) -Rs^2/(Rc^2-Rs^2)*sum(wnx) ...
        +Rs^2/(Rc^2-Rs^2)*  sum(wnx.*exp( -cns*ep1/snr/(alpha1-ep1*alpha2)  ).*Liinters) -...
         Rc^2/(Rc^2-Rs^2)* sum(wnx.*exp( -cnc*ep1/snr/(alpha1-ep1*alpha2)  ).*Liinterc);
    
    % second message at the near user
    
        mu = cns/(min(alpha2/ep2,(alpha1 - ep1*alpha2)/ep1));
    Liinter = exp( -2*pi*lamb* (mu).^(2/al) /al * beta(2/al, (al-2)/al)  );
    pa2(k) = sum(wnx) - sum(wnx.*exp( -cns/snr/min(alpha2/ep2,(alpha1-ep1*alpha2)/ep1)  ).*Liinter);
end
semilogy(snrdb, p3, snrdb,p1,snrdb,pa,snrdb,p4,snrdb,p2,snrdb,pa2)