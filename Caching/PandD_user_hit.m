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

 snrdb = [10: 10 : 40]; ct=50000;
 K=2; Rc = 50;  al=3;   rmax = 3000;    lamb=0.01/pi/Rc^2; deltax=1.1;
 R0=0.5; R1 = 3; R2=3.5;R3=11;
 beta0 = 4/8; beta1=3/8; beta2=2/8; beta3=1/8; 
 N=20; %Gassus variables
 for k = 1: length(snrdb)
    snr = 10^((snrdb(k)-30)/10)*eta/sigma_square*10; 
    zeta0 =    beta0 - (2^(R0/4)-1)* ( beta1+beta2+beta3); 
    zeta1 =    beta1 - (2^(R1/4)-1)* ( beta2+beta3); 
    zeta2 =    beta2 - (2^(R2/4)-1)* ( beta3); 
    zeta3 =    beta3; 
   sum1=zeros(K); sum21=0; sum41=0;sum22=0; sum42=0;sum23=0; sum43=0;sum20=0; sum40=0;
   sum3=0;  sum5=0; sum6=0;  
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
      cx1 = cx1(t2); cy1 = cy1(t2); %now all the base statins are ordered. 

     bsi = 3; %%% only focus on this particular cell

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
    if zmk< (2^(R0/4)-1)/zeta0/snr     
        sum1 = sum1+1;
    end
    if zmk<(2^(R0)-1)/snr
        sum3 = sum3 +1; %oma
    end
    
    %the performance at base station  
    %%%%% file 0
    dbs = sqrt(cx1(bsi,1)^2 + cy1(bsi,1)^2 );    
    if  dbs^al>min([zeta0*snr/(2^(R0/4)-1)])
        sum20 = sum20+1;
    end
    if 1/dbs^al<(2^(R0)-1)/snr  
        sum40 = sum40+1;
    end
    %%%%% file 1
    if  dbs^al>min([zeta0*snr/(2^(R0/4)-1) zeta1*snr/(2^(R1/4)-1) ])
        sum21 = sum21+1;
    end
    if 1/dbs^al<(2^(R1)-1)/snr  
        sum41 = sum41+1;
    end
     %%%%% file 2
   if  dbs^al>min([zeta0*snr/(2^(R0/4)-1) zeta1*snr/(2^(R1/4)-1) zeta2*snr/(2^(R2/4)-1)  ])
        sum22 = sum22+1;
    end
    if 1/dbs^al<(2^(R2)-1)/snr  
        sum42 = sum42+1;
    end
    %%%%% file 0
    if  dbs^al>min([zeta0*snr/(2^(R0/4)-1) zeta1*snr/(2^(R1/4)-1) zeta2*snr/(2^(R2/4)-1) zeta3*snr/(2^(R3/4)-1)])
        sum23 = sum23+1;
    end
    if 1/dbs^al<(2^(R3)-1)/snr  
        sum43 = sum43+1;
    end
    
    
    
    
    
  end
    p1(k) = sum1(1)/ct;     
    o1(k) = sum3/ct;  
    
   p20(k) = sum20/ct;     
    o20(k) = sum40/ct;     
   p21(k) = sum21/ct;     
    o21(k) = sum41/ct;     
   p22(k) = sum22/ct;     
    o22(k) = sum42/ct;  
       p23(k) = sum23/ct;     
    o23(k) = sum43/ct;     

    %theorethic  
    ep1 = 2^(R1/2)-1;
    sumx1 = 0;
    ctt=10;  %%%%%%%%%%%%%%%%%% need to change 
    stepz =0.1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% need to change 
    zx = [deltax*Rc : stepz  :deltax*Rc+stepz*ctt*20]; %%%%%%%%%%%%%% need to change XXX
    for i =1: length(zx)
        z = zx(i);
        stepr = Rc/ctt;
        rx = [z-Rc+0.000000000000000000001: stepr: z+Rc-0.0000000000000000000001];
        for j =1: length(rx)
            r = rx(j);
            sumx1 = sumx1 + exp(-(2^(R0/4)-1)*r^al/snr/zeta0)*  2*r*acos((z^2+r^2-Rc^2)/2/z/r)/pi/Rc^2  *stepr...
                *2*lamb^bsi*pi^bsi*z*(z^2-deltax^2*Rc^2)^(bsi-1)/factorial(bsi-1)...
                *exp(-lamb*pi*(z^2-deltax^2*Rc^2)) *stepz;
        end
    end
    pa(k) = 1- sumx1 ;
    
    %file 3
    sum2x = 0;
    taum = 1/(min([zeta0*snr/(2^(R0/4)-1) zeta1*snr/(2^(R1/4)-1) zeta2*snr/(2^(R2/4)-1) zeta3*snr/(2^(R3/4)-1)]))^(1/al); 
    temp1 = lamb*pi*(1/taum^2 - deltax^2*Rc^2);
    for l =0: bsi-1
        sum2x = sum2x + temp1^l/factorial(l)*exp(-temp1);
    end
    pa23(k) = sum2x;
    %file 2
    sum2x = 0;
    taum = 1/(min([zeta0*snr/(2^(R0/4)-1) zeta1*snr/(2^(R1/4)-1) zeta2*snr/(2^(R2/4)-1) ]))^(1/al); 
    temp1 = lamb*pi*(1/taum^2 - deltax^2*Rc^2);
    for l =0: bsi-1
        sum2x = sum2x + temp1^l/factorial(l)*exp(-temp1);
    end
    %file 1
    pa22(k) = sum2x;
    sum2x = 0;
    taum = 1/(min([zeta0*snr/(2^(R0/4)-1) zeta1*snr/(2^(R1/4)-1)  ]))^(1/al); 
    temp1 = lamb*pi*(1/taum^2 - deltax^2*Rc^2);
    for l =0: bsi-1
        sum2x = sum2x + temp1^l/factorial(l)*exp(-temp1);
    end
    pa21(k) = sum2x;
    %file 0
    sum2x = 0;
    taum = 1/(min([zeta0*snr/(2^(R0/4)-1)]))^(1/al); 
    temp1 = lamb*pi*(1/taum^2 - deltax^2*Rc^2);
    for l =0: bsi-1
        sum2x = sum2x + temp1^l/factorial(l)*exp(-temp1);
    end
    pa20(k) = sum2x;
    
       hm(k) = (1-p21(k))*pf(1) + (1-p22(k))*pf(2) + (1-p23(k))*pf(3);
       ham(k) = (1-pa21(k))*pf(1) + (1-pa22(k))*pf(2) + (1-pa23(k))*pf(3);
       ho(k) = (1-o21(k))*pf(1) + (1-o22(k))*pf(2) + (1-o23(k))*pf(3);
    
end
semilogy( snrdb,p1,snrdb,pa )%semilogy(snrdb,o1,snrdb,p1,snrdb,pa,snrdb,ho,snrdb,hm,snrdb,ham)