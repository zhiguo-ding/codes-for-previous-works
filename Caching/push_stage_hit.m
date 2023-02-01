clear all%function [p1, p2] = fund_mac(K,ct)
fc = 2e9; %Carrier frequency (Hz)
c= 3e8; %speed of light (m/s)
wavelength = c/fc; % in m
eta= (wavelength/(4*pi))^2;
BW=10*10^6; %10 MHz 
Nf=10;%dB
sigma2_dbm= -180+10*log10(BW)+Nf; %Thermal noise in dBm
sigma_square=10^((sigma2_dbm-30)/10);

fl = [1:10];
gammaf = 1.5;
pf = 1./fl.^gammaf./sum(1./fl.^gammaf);

 snrdb = [10: 5 : 40]; ct=50;
 K=2; Rc = 50; al=3; R1 = 1; R2=1; R3=1;  rmax = 5500;    lamb=0.01/pi/Rc^2; 
 bsi = 1; %which base staiton
 t = 5; %threshold
 F = 3; %number of files
 beta1 = 3/4; beta2 = 1/4;
 N=20; %Gassus variables
 m=1;
  
 for k = 1:length(snrdb)
    snr = 10^((snrdb(k)-30)/10)*eta/sigma_square*10;  
    ep1 = 2^R1-1; ep2=2^R2-1; ep3=2^R3-1;  
     
    xi1 = beta1 - ep2*beta2; %file 2, so it is related to ep2
    xi2=beta2; %for file 3
    
  sum1=zeros(t,F); sum2=0; sum3=0; sum4=0; sum5=0; sum6=0;  
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
     %order these base stations to the comp user
     [t1,t2] = sort(dx,'ascend');
     %bsn_small = abs(sum((sign(t1-rsmall^2)-1)/2)); %the number of bs in the small circle.
     cx1 = cx1(t2); cy1 = cy1(t2); %now all the base statins are ordered. 
     dx_order = sqrt(t1);
     
     alpha1 = min(1, ep1*(snr/dx_order(t)^al +1)/snr/(1+ep1)*dx_order(t)^al );
     Pr = 1-alpha1;
     
     for tx = 1: t %the t nearest base stations
         rate1 = log2(1+snr*alpha1/dx_order(tx)^al / (snr*(1-alpha1)/dx_order(tx)^al +1) );
         rate2 = log2(1+snr*beta1*Pr/dx_order(tx)^al / (snr*beta2*Pr/dx_order(tx)^al +1) );
         rate3 = log2(1+snr*beta2*Pr/dx_order(tx)^al / ( 1) );
         if rate1<R1-0.0000001
             sum1(tx,1) = sum1(tx,1)+1; sum1(tx,2) = sum1(tx,2)+1; sum1(tx,3) = sum1(tx,3)+1;
         elseif rate2<R2
             sum1(tx,2) = sum1(tx,2)+1; sum1(tx,3) = sum1(tx,3)+1;
         elseif rate3<R3
             sum1(tx,3) = sum1(tx,3)+1;
         end 

     end
     %test
  end
  
    pt1(k) = sum1(t,1)/ct; %user t to decode f1 
    pt2(k) = sum1(t,2)/ct; %user t to decode f2 
    pt3(k) = sum1(t,3)/ct; %user t to decode fi
    pm1(k) = sum1(m,1)/ct; %user 2 to decode f3
    pm2(k) = sum1(m,2)/ct; %user 2 to decode f2
    pm3(k) = sum1(m,3)/ct; %user 2 to decode f3
    %theorethic 
    %user t to decode f1
    xc1 = snr/ep1;
    sumx = 0;
    for kl = 0: t-1
        sumx = sumx +  xc1^(2*kl/al)/factorial(kl)*(lamb*pi)^(kl);
    end
    pat1(k) =   exp(-lamb*pi*xc1^(2/al))*sumx;

        %user t to decode f2
    xc1 = ep1/snr+(1+ep1)/snr/min(xi1/ep2); %f2
    sumx = 0;
    for kl = 0: t-1
        sumx = sumx +  xc1^(-2*kl/al)/factorial(kl)*(lamb*pi)^(kl);
    end
    pat2(k) = exp(-lamb*pi*xc1^(-2/al))*sumx;

    
    %user t to decode f3
    xc1 = ep1/snr+(1+ep1)/snr/min(xi1/ep2,xi2/ep3); %f3
    sumx = 0;
    for kl = 0: t-1
        sumx = sumx +  xc1^(-2*kl/al)/factorial(kl)*(lamb*pi)^(kl);
    end
    pat3(k) = exp(-lamb*pi*xc1^(-2/al))*sumx;

    %user m to decode f1
    xc1 = ep1/snr;%+epi*(1+epi)/snr/min(xi1,xi2); %f3
     
    sumx = 0;
    for kl = 0: m-1
        sumx = sumx +  xc1^(-2*kl/al)/factorial(kl)*(lamb*pi)^(kl);
    end
    pam1(k) = exp(-lamb*pi*xc1^(-2/al))*sumx;

    %      %%%%%%%%%  % user m to decode f3 m=2
    xc1 = snr/ep1; %the first part excluding Q1
    sumx = 0;
    for kl = 0: t-1
        sumx = sumx +  xc1^(2*kl/al)/factorial(kl)*(lamb*pi)^(kl);
    end
       ximin = min(xi1/ep2);%min(xi1,xi2); %f2
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
      pam2(k) = Q1 + exp(-lamb*pi*xc1^(2/al))*sumx;
      
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
      
      ht(k) = (1-pt1(k))*pf(1) + (1-pt2(k))*pf(2) + (1-pt3(k))*pf(3);
      hat(k) = (1-pat1(k))*pf(1) + (1-pat2(k))*pf(2) + (1-pat3(k))*pf(3);
      hm(k) = (1-pm1(k))*pf(1) + (1-pm2(k))*pf(2) + (1-pm3(k))*pf(3);
      ham(k) = (1-pam1(k))*pf(1) + (1-pam2(k))*pf(2) + (1-pam3(k))*pf(3);
      hto(k) = (1-pat1(k))*pf(1) ;
      hmo(k) = (1-pam1(k))*pf(1) ;      
   end
semilogy(snrdb,ht, snrdb,hm,snrdb, hat, snrdb, ham )