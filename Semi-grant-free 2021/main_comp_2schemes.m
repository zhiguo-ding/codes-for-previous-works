clear all;

ct=10000;
snrdb = [0: 5: 30];
R0 = 1; M=1;
Rs = 1.5;
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
        rhybx(i) = Rsgft;
        
         % always scheme II
        rii = log2(1+Ps*htoder(M)/(1+P0*g));
        if rii<Rs
            sum6 = sum6 +1;
        end  
        riix(i) = rii;
        
        %always scheme I
        taux=1;%snr;%XXXX HOW TO SET TAU? need to be large
%         if htoder(1)<taux %if the weakest grant free user below the threshold
            if log2(1 + g*P0/(1+Ps*htoder(1))) > R0 %the first step of SIC is ok
                rix(i) = log2(1+Ps*htoder(1));
                if rix(i)<Rs
                    sum7 = sum7 +1;
                end
            else
                sum7 = sum7 +1;
            end
%         else
%             sum7 = sum7+1;
%         end
        
        %always the best gf user
        rbest(i) = log2(1+Ps*htoder(M));  
              
    end
    ps(k)= sum1/ct; 
    pii(k) = sum6/ct; 
    pi(k) = sum7/ct;
       
    %rate
    Rh(k) = mean(rhybx);
    Rii(k) = mean(riix);
    Ri(k) = mean(rix);
        
end
%semilogy( snrdb,pa,snrdb,paa)
%semilogy(snrdb,pi,'-d',snrdb,pii,snrdb,ps )
plot(snrdb,Ri,'-d',snrdb,Rii,snrdb,Rh)
%semilogy(snrdb,poma,'-d',snrdb,p0,'-x',snrdb,pi,snrdb,pii,snrdb,ps)
 
 