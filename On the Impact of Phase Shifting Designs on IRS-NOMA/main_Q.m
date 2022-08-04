clear all;
 

 
sigma2_dbm= -94;%+10*log10(BW)+Nf; %Thermal noise in dBm
sigma_square=10^((sigma2_dbm-30)/10);

ct=20000000;
N=64;
snrdbm = [0: 2: 30];
c1=0.8; c2 =1-c1;
d2=10; dr = d2; dr1 =15; d12 = dr1;  
al=4;
R1 = 1.5;
Q=8;

betadr_db = -35.1-36.7*log10(dr) +20;
betadr = 10^(betadr_db/10);
betadr1_db = -35.1-36.7*log10(dr1) +20;
betadr1 = 10^(betadr1_db/10);

for k = 1: length(snrdbm)
    snr = 10^((snrdbm(k)-30)/10)/sigma_square;
    P1=snr; P2=snr; Pr=snr;
    P0 = (P1+Pr)/2;
    Ps= (P1+P2+Pr)/3;
    sum1 = 0; sum2=0; sum3=0; sum4=0;
    for i =1 : ct
        g1=complex(randn(N,1)*sqrt(0.5),randn(N,1)*sqrt(0.5)); 
        g0=complex(randn(N,1)*sqrt(0.5),randn(N,1)*sqrt(0.5));     
%        had = hadamard(N);
        th1 = randsrc(N,N,[1 -1]);%had(:,1);
        %thx = rand(N,N);
        %th1 = exp(-complex(0,1)*2*pi*thx);
        
%         for n = 1 : Q
%             x = abs(g1'*diag(th1(:,n))*g0)^2*betadr*betadr1;
%             %noma
%             nomaratex(n) = log2(1 + Ps*c1*x/(Ps*c2*x +1));
%             
%             %oma-irs
%             irsomaratex(n) = 1/2*log2(1 + Ps*x );
%         end
        x = max(abs((g1.*g0).'*th1(:,1:Q)).^2*betadr*betadr1);
        
            znoma = log2(1 + Ps*c1*x/(Ps*c2*x +1));
            
            %oma-irs
            zoma = 1/2*log2(1 + Ps*x );
        if znoma< R1
            sum1 = sum1 +1;
        end

        if zoma< R1
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
    end
    pnoma(k)= sum1/ct;
    pioma(k)= sum2/ct;
    poma(k)= sum3/ct;
    
    %analysis
    cc1 = 2^(4*R1)-1;
    cc2 = 2^(2*R1)-1;
    cc3 = 2^(R1)-1;
    pomaa(k) = 1 -exp(-1/betadr*cc1/P1) + exp(-1/betadr*cc1/P1) *...
        ( 1 -exp(-1/betadr1*cc1/Pr)   );
    piomaa(k) = (1 - exp(-1/betadr/betadr1*cc2/P0/N));
    pnomaa(k) = (1 - exp(-1/betadr/betadr1*cc3/N/Ps/(c1-c2*cc3)));
end
%semilogy(snrdbm,poma,'-d',snrdbm,pioma,snrdbm,pnoma,snrdbm,pomaa,snrdbm,piomaa,snrdbm,pnomaa)
%semilogy( snrdbm,pioma,snrdbm,pnoma, snrdbm,piomaa,snrdbm,pnomaa)
semilogy(snrdbm,poma,snrdbm,pioma,snrdbm,pnoma)
%semilogy( snrdbm,pioma,snrdbm,pnoma)



%  fc = 2e9; %Carrier frequency (Hz)
% c= 3e8; %speed of light (m/s)
% wavelength = c/fc; % in m
% eta= (wavelength/(4*pi))^2;
% BW=10*10^6; %10 MHz 
% Nf=10;%dB
% sigma2_dbm= -170+10*log10(BW)+Nf; %Thermal noise in dBm
% sigma_square=10^((sigma2_dbm-30)/10);
 