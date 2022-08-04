clear all;
 
 
 
sigma2_dbm= -70;%+10*log10(BW)+Nf; %Thermal noise in dBm
sigma_square=10^((sigma2_dbm-30)/10);

ct=200000;
N=16;
snrdbm = [0: 2: 30];
c1=0.8; c2 =1-c1;
d2=10; dr = d2; dr1 =10; d12 = dr1;  
al=4;
R1 = 1;
Q=64;

betadr_db = -35.1-36.7*log10(dr) +20;
betadr = 10^(betadr_db/10);
betadr1_db = -35.1-36.7*log10(dr1) +20;
betadr1 = 10^(betadr1_db/10);


for k = 1 : length(snrdbm)
    snr = 10^((snrdbm(k)-30)/10)/sigma_square;
    P1=snr; P2=snr; Pr=snr;
    P0 = (P1+Pr)/2;
    Ps= (P1+P2+Pr)/3;
    sum1 = 0; sum2=0; sum3=0; sum4=0;
    for i =1 : ct
        g1=complex(randn(N,1)*sqrt(0.5),randn(N,1)*sqrt(0.5)); 
        g0=complex(randn(N,1)*sqrt(0.5),randn(N,1)*sqrt(0.5));     
         
             x = sum(abs(g1.*g0))^2*betadr*betadr1;%/dr^al/dr1^al;
            %noma
            nomaratex = log2(1 + Ps*c1*x/(Ps*c2*x +1));
            
            %oma-irs
            irsomaratex = 1/2*log2(1 + Ps*x );
         if  (nomaratex)< R1
            sum1 = sum1 +1;
        end

        if  (irsomaratex)< R1
            sum2 = sum2 +1;
        end

        %conventional relaying
        h12=complex(randn(1,1)*sqrt(0.5),randn(1,1)*sqrt(0.5)); 
        h2=complex(randn(1,1)*sqrt(0.5),randn(1,1)*sqrt(0.5)); 
        if 1/4*log2(1+P1*abs(h2)^2*betadr)<R1
            sum3= sum3+1;
        elseif 1/4*log2( 1+Pr*abs(h12)^2*betadr1)<R1
            sum3= sum3 +1;
        end
                
        %NOMA relaying
        if 1/2*log2(1+Ps*c1*abs(h2)^2*betadr/(1+Ps*c2*abs(h2)^2*betadr))<R1
            sum4= sum4+1;
        elseif 1/2*log2( 1+Pr*abs(h12)^2*betadr1)<R1
            sum4= sum4 +1;
        end

    end
    pnoma(k)= sum1/ct;
    pioma(k)= sum2/ct;
    poma(k)= sum3/ct;
    pnomarelay(k)= sum4/ct;
    
    %analysis
    %analysis
    cc1 = 2^(4*R1)-1;
    cc2 = (2^(2*R1)-1)/betadr1/betadr/P0;
    mu1 = pi/4;
    sig1 = sqrt(1-pi^2/16);
    psi1 = sqrt(N)*(sqrt(cc2)/N-mu1);
    
    cc3 = 1/betadr/betadr1*(2^(R1)-1)/Ps/(c1-c2*(2^(R1)-1));
    psi2 = sqrt(N)*(sqrt(cc3)/N-mu1);
 
    pomaa(k) = 1 -exp(-1/betadr*cc1/P1) + exp(-1/betadr*cc1/P1) *...
        ( 1 -exp(-1/betadr1*cc1/Pr)   );
    piomaa(k) =  0.5+0.5*erf(psi1/sqrt(2)/sig1);%normcdf(psi1,0,sig1) ;
    pinomaa(k) =  normcdf(psi2,0,sig1) ;
    
    Nb=N/2;
    pou(k) = 2^N*pi^(N/2)*(gamma(3/2))^N/factorial(3*Nb-1)...
        *2^(-3*Nb)*gammainc(2*sqrt(cc2),3*Nb)*gamma(3*Nb);
    pnu(k) = 2^N*pi^(N/2)*(gamma(3/2))^N/factorial(3*Nb-1)...
        *2^(-3*Nb)*gammainc(2*sqrt(cc3),3*Nb)*gamma(3*Nb);
      
end
%  semilogy(  snrdbm,pioma,snrdbm,pnoma,'-*', snrdbm,piomaa ...
%     ,snrdbm,pinomaa,'-o',snrdbm,min(1,pou),snrdbm,min(1,pnu),'-d' )
% semilogy(snrdbm,poma, snrdbm,pioma,snrdbm,pnoma,'-*',snrdbm,pomaa,snrdbm,piomaa ...
%    ,snrdbm,pinomaa,'-o',snrdbm,min(1,pou),snrdbm,min(1,pnu),'-d' ) %semilogy(snrdbm,pioma,snrdbm,pnoma, snrdbm,piomaa ...
 %    ,snrdbm,pinomaa,snrdbm,min(1,pou),snrdbm,min(1,pnu) ) 
 %semilogy(snrdbm,poma,snrdbm,pnomarelay,snrdbm,pioma,snrdbm,pnoma,'-*' )
  semilogy( snrdbm,pioma,snrdbm,pnoma,'-*' )
