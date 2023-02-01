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
gammaf = 0.5;
pf = 1./fl.^gammaf./sum(1./fl.^gammaf);

 snrdb = [10: 10 : 40]; ct=5000;
 K=2; Rc = 50;  al=3;   rmax = 3000;    lamb=0.01/pi/Rc^2; deltax=1.1;
 R1 = 1; R2=6;
 beta2=1/4; beta1=1/4;
 N=20; %Gassus variables
 for k = 1: length(snrdb)
    snr = 10^((snrdb(k)-30)/10)*eta/sigma_square*10; 
    taux =   (2^(R1/2)-1)/snr/(1-beta2*2^R1); taux2 = (2^(R2/2)-1)/snr/beta2;
  sum1=zeros(K); sum2=zeros(K); sum3=0; sum4=0; sum5=0; sum6=0;  
  for ix = 1 :ct       
      %get the locations for the base stations
     bsn  = poissrnd(pi*rmax^2*lamb);
     pppind = 1;
     cx1=zeros(bsn,1); cy1=zeros(bsn,1); dx=zeros(bsn,1);
     while pppind<=bsn                 
         cx1(pppind,1) = sign(randn(1,1))*rand(1,1)*rmax;
         cy1(pppind,1) = sign(randn(1,1))*rand(1,1)*rmax;
         dx(pppind,1) = (cx1(pppind))^2+(cy1(pppind))^2;
         if dx(pppind,1)<rmax^2 &dx(pppind,1)>deltax^2*Rc^2
             pppind = pppind+1;
         end
     end
     %order these base stations to the comp user
     [t1,t2] = sort(dx,'ascend');
     if min(sqrt(dx))<Rc
         dd=0;
     end
     cx1 = cx1(t2); cy1 = cy1(t2); %now all the base statins are ordered. 

     bsi = 5; %%% only focus on this particular cell

     %get the relative locations for the users      
     cx2=zeros(bsn,K); cy2=zeros(bsn,K); %each rwo contains K users  relative locations
 
     %the  user in this cell
     pppind = 1;         
         while pppind<=1                 
             cx2(1,pppind) = sign(randn(1,1))*rand(1,1)*Rc;
             cy2(1,pppind) = sign(randn(1,1))*rand(1,1)*Rc;
             dx2 = (cx2(1,pppind))^2+(cy2(1,pppind))^2;
             if dx2<Rc^2  
                 pppind = pppind+1;
             end
         end 
             cx3(bsi,1) = cx2(1)+cx1(bsi); cy3(bsi,1) = cy2(1)+cy1(bsi); %actual user locations
 
         
    %the performance at user  
    dd = sqrt(cx3(bsi,1)^2 + cy3(bsi,1)^2 );
    zmk = abs(complex(sqrt(0.5)*randn(1,1),sqrt(0.5)*randn(1,1)))^2/ dd^al;
    if zmk< taux
        sum1 = sum1+1;
    end
    if zmk<(2^(R1)-1)/snr
        sum3 = sum3 +1; %oma
    end
    
    %the performance at base station  
    dbs = sqrt(cx1(bsi,1)^2 + cy1(bsi,1)^2 );    
    if 1/dbs^al<taux | 1/dbs^al<taux2
        sum2 = sum2+1;
    end
    if 1/dbs^al<(2^(R2)-1)/snr  
        sum4 = sum4+1;
    end
    
  end
    p1(k) = sum1(1)/ct;     
    p2(k) = sum2(1)/ct;     
    o1(k) = sum3/ct;     
    o2(k) = sum4/ct;     
    %theorethic  
    ep1 = 2^(R1/2)-1;
    sumx1 = 0;
    ctt=5;  %%%%%%%%%%%%%%%%%% need to change 
    stepz =0.5; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% need to change 
    zx = [deltax*Rc : stepz  :deltax*Rc+stepz*ctt*20]; %%%%%%%%%%%%%% need to change XXX
    for i =1: length(zx)
        z = zx(i);
        stepr = Rc/ctt;
        rx = [z-Rc+0.000000000000000000001: stepr: z+Rc-0.0000000000000000000001];
        for j =1: length(rx)
            r = rx(j);
            sumx1 = sumx1 + exp(-taux*r^al)*  2*r*acos((z^2+r^2-Rc^2)/2/z/r)/pi/Rc^2  *stepr...
                *2*lamb^bsi*pi^bsi*z*(z^2-deltax^2*Rc^2)^(bsi-1)/factorial(bsi-1)...
                *exp(-lamb*pi*(z^2-deltax^2*Rc^2)) *stepz;
        end
    end
    pa(k) = 1- sumx1 ;
    
    sum2x = 0;
    taum = (min( 1/taux, 1/taux2))^(1/al); 
    temp1 = lamb*pi*(taum^2 - deltax^2*Rc^2);
    for l =0: bsi-1
        sum2x = sum2x + temp1^l/factorial(l)*exp(-temp1);
    end
    pa2(k) = sum2x;
   
end
semilogy(snrdb,o1,snrdb,p1,snrdb,pa,snrdb,o2,snrdb,p2,snrdb,pa2)