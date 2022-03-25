clear all;

ct=300;
snrdb = [0: 5: 30];
R0 = 1; M=5;
Rs = 0.9;
eps = 2^Rs-1; ep0= 2^R0 -1;
sum2=zeros(length(snrdb),1);
sum3=zeros(length(snrdb),ct);

for k = 1 : length(snrdb)
    snr = 10^((snrdb(k))/10);
    P0=snr; Ps=snr; alpha0=ep0/P0; alphas=eps/Ps;
    alpha1 =(1+eps)*alpha0;
    alpha2 = ep0*(eps+1)/P0/(1-ep0*eps);
    sum1 = 0;   sum4=0; sum5 = 0; sum6=0; sum7=0;
    rhybx=0; riix=0; rix=0;
    for i =1 : ct
        g=complex(randn(1,1)*sqrt(0.5),randn(1,1)*sqrt(0.5)); 
        g=abs(g)^2;
        h=complex(randn(M,1)*sqrt(0.5),randn(M,1)*sqrt(0.5));     
        tau = max(0, P0*g/(2^R0-1)-1);  
        
        htt = (abs(h)).^2;
        htoder = sort(htt,'ascend');
        if htoder(1)>tau/Ps % scheme II
            Rsgft = log2(1+Ps*htoder(M)/(1+P0*g));
            if Rsgft<Rs
                sum1 = sum1 +1;
            end
            
        else % scheme I
            tt1 = (sign(htoder*Ps-tau)+1)/2;
            tt2 = sum(tt1); %number above tau
            hsel = htoder(M-tt2);
            Rsgft = max(log2(1+Ps*hsel),log2(1+Ps*htoder(M)/(1+P0*g))) ;
            if Rsgft<Rs
                sum1 = sum1 +1;
            end
            
        end
              
    end
    ps(k)= sum1/ct; 
    
    %analytical
    if M>=2
    % Qm first, 1\leq m \leq M-2
    for m =1 : M-2
        etam = factorial(M)/factorial(m-1)/factorial(M-m-2);
        sum11 = 0; sum12=0;
        for l = 0 : M-m
            etatempl = factorial(M-m)/factorial(l)/factorial(M-m-l);
            for p = 0 : m
                etatempp = factorial(m)/factorial(p)/factorial(m-p);
                
                  mu2 = (l*alphas+(M-m-l)/Ps/ep0)*P0;
                 mu4 = exp(-l*alphas+(M-m-l)/Ps); 
                
                sum11 = sum11 + etam/m/(M-m-1)/(M-m) * etatempl * (-1)^l * etatempp*(-1)^p...
                    *(mu4* phimy(p,mu2,alpha0,alpha1,alpha2,alphas,Ps) );
            end
        end
        qm(m) = sum11;
    end
    
    %QM-1
    sum12=0;
    etaM1 = factorial(M)/factorial(M-2);
    for i = 0 : M-1
        mu7 = P0/Ps/ep0; mu8 = alphas*P0;
        sum12 = sum12 + factorial(M-1)/factorial(i)/factorial(M-1-i)...
            *(-1)^i*etaM1/(M-1)*...
            (exp(1/Ps)*phimy(i,mu7,alpha0,alpha1,alpha2,alphas,Ps)...
            -exp(-alphas)*phimy(i,mu8,alpha0,alpha1,alpha2,alphas,Ps));
    end        
    qm(M-1) = sum12;
    
    %QM
    sum13 = 0;
    for i = 0 : M
        sum13 = sum13 + factorial(M)/factorial(i)/factorial(M-i)...
            *(-1)^i*exp(i/Ps)*gmu(i/alpha0/Ps , alpha0, alpha0*(1+eps));
    end
    qm(M) = sum13 + (1-exp(-alphas))^M*exp(-alpha0*(1+eps));
    
    %Q0
        sum14 = 0;
    for l = 0 : M             
            mu12 = l*alphas*P0+(M-l)/Ps/alpha0;  

        sum14 = sum14 + factorial(M)/factorial(M-2)/M/(M-1)...
            *factorial(M)/factorial(l)/factorial(M-l)*(-1)^l...
            *exp(-l*alphas)*exp((M-l)/Ps)*gmu(mu12, alpha0,alpha2);
    end
    qm(M+1) = sum14; %actually it is Qm(0)

 
    
    %the last one
    sum15 = 0;
    for i = 0 : M
        sum15 = sum15 + factorial(M)/factorial(M-i)/factorial(i)*(-1)^i...
            *exp(-i*alphas)*(1-exp(-(1+i*alphas*P0)*alpha0))/(1+i*alphas*P0);
    end
    Qx = sum15;
    pa(k) = sum(qm) +Qx;
    else
        Q0 = exp(1/Ps)*gmu(1/Ps/alpha0, alpha0,alpha2)...
                - exp(-alphas)*gmu(alphas*P0, alpha0,alpha2) ; 
            %QM
        sum13 = 0;
        for i = 0 : M
            sum13 = sum13 + factorial(M)/factorial(i)/factorial(M-i)...
                *(-1)^i*exp(i/Ps)*gmu(i/alpha0/Ps , alpha0, alpha0*(1+eps));
        end
        qm(M) = sum13 + (1-exp(-alphas))^M*exp(-alpha0*(1+eps));
        %the last one
        sum15 = 0;
        for i = 0 : M
            sum15 = sum15 + factorial(M)/factorial(M-i)/factorial(i)*(-1)^i...
                *exp(-i*alphas)*(1-exp(-(1+i*alphas*P0)*alpha0))/(1+i*alphas*P0);
        end
        Qx = sum15;
        pa(k) = Q0 + qm(M) +Qx;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %approximation
    if M>=2
    % Q0
    suma1=0;
    alx2=(eps+1)/(1-ep0*eps);
    for i =0 : M     
        suma1 = suma1 + factorial(M)/factorial(M-2)/Ps^(M+1)/M/(M-1)...
            *exp(M/Ps)*factorial(M)/factorial(M-i)/factorial(i)...
            *(eps+1)^(M-i)*(eps-1/ep0)^i*ep0^(i+1)*(alx2^(i+1)-1)/(i+1);
    end
    qma(M+1) = suma1;
    
    %Qm 1\leq m\leq M-2    
    for m = 1 : M-2
        suma2=0; suma3=0;
        etabar = factorial(M)/factorial(m-1)/factorial(M-m-2)...
            /m/(M-m-1)/(M-m);
        for i = 0 : M-m
            suma2 = suma2 + etabar/Ps^(M+1)*eps^m*...
                factorial(M-m)/factorial(i)/factorial(M-m-i)...
                *(eps+1)^(M-m-i)*(eps-1/ep0)^i*ep0^(i+1)...
                *(alx2^(i+1) - (1+eps)^(i+1))/(i+1);
            suma3 = suma3 + etabar/Ps^(M+1)/ep0^m*...
                factorial(M-m)/factorial(i)/factorial(M-m-i)...
                *(eps+ep0*eps)^(M-m-i)*(eps-1/ep0)^i...
                *ep0^(m+i+1)*eps^(m+i+1)/(m+i+1);
        end
        qma(m) = suma2+suma3;
    end
    
    %Q_M-1
    qma(M-1) = factorial(M)/factorial(M-2)/(M-1)/Ps^(M+1)...
        *(1+eps)*eps^(M-1)*ep0*(alpha2/alpha0-1-eps+eps/M);
    
    %QM
    qma(M) = 1/(M+1)/ep0^M*alpha0^(M+1)*eps^(M+1)+alphas^M;
    
    %Qx
    Qxa = alphas^M * ((1+ep0)^(M+1)-1)/(M+1)/P0;
    
    paa(k) =  min(1,sum(qma) +Qxa);
    paax(k) =  min(1,alphas^M );
    else
        alx2=(eps+1)/(1-ep0*eps);
        Q0x = 1/Ps^2*(1+eps)*eps*(alx2-1);
        %QM
        qM = 1/(M+1)/ep0^M*alpha0^(M+1)*eps^(M+1)+alphas^M;
        Qxa = alphas^M * ((1+ep0)^(M+1)-1)/(M+1)/P0;
        paa(k) = min(1, Q0x+qM  +Qxa);
        paax(k) =  min(1,alphas^M );
    end
        
end
%semilogy( snrdb,pa,snrdb,paa,snrdb,paax)
%semilogy( snrdb,pa)
semilogy(snrdb,ps,snrdb,pa)
%semilogy(snrdb,poma,'-d',snrdb,p0,'-x',snrdb,pi,snrdb,pii,snrdb,ps)
 
 