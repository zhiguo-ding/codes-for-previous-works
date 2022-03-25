clear all 
 snrdb = [0: 5 : 40]; ct=10000;
 M=10; 
 
 
  R0 = 1; %snr0 = (2^R0-1)/h0;
  Ri=0.5;
   
   
 for k =  1:length(snrdb)
    snr = 10^((snrdb(k))/10);  %snr0 = snr;
    snr0=10^((snrdb(5))/10); 
    eps0 = (2^R0-1)/snr0;          %%%%% BE CAREFUL WITH snr0
    tau =  1.2;%/snr;         
        
   sum1=0; sum2=0; sum3=0;sum4=0;  
  for ix = 1 :ct       
      h0 = abs(complex(sqrt(0.5)*randn(1,1),sqrt(0.5)*randn(1,1))).^2;
      
     h = abs(complex(sqrt(0.5)*randn(M,1),sqrt(0.5)*randn(M,1))).^2;
     h = sort(h, 'ascend');     
     
     n = sum(sign(h-tau)+1)/2;
     nk(ix)=n;
     xm = sum(h(M-n+1:M));
     
     Rs = log2(1+xm*snr/(h0*snr0 +1 ));
     
     if (Rs==0|Rs>n*Ri) & log2(1+h0*snr0)>R0
         sum1 = sum1 +1;
     end
     
     % SINGLE
     if n==0
         xm = 0;
     else
         xm = sum(h(M));%sum(h(M-n+1:M));
     end     
     Rs = log2(1+xm*snr/(h0*snr0 +1 ));
     
     %use the best
     if (Rs==0|Rs>Ri) & log2(1+h0*snr0)>R0
         sum2 = sum2 +1;
     end

          Rall = log2(1+sum(h)*snr/(h0*snr0 +1 ));
     if (Rall==0|Rall>M*Ri) & log2(1+h0*snr0)>R0
         sum3 = sum3 +1;
     end
  end
  p1(k) = 1- sum1/ct;
  p2(k) = 1- sum2/ct;
  p3(k) = 1- sum3/ct;
  nkk(k) = mean(nk);
   
    %analytical    
    suma1 = 0;
    for n=1:M
        eps_s = 2^(n*Ri)-1;
        tau_n = (n*tau/eps_s*snr-1)/snr0;
        tau_bar = max(tau_n, eps0);
        Q5 =0;
        for l = 1 : n-1
            Q5 = Q5 + eps_s^l*snr^(-l)*snr0^l/factorial(l)...
                *exp(-tau_n)*gammainc((tau_bar-tau_n)*(1+eps_s*snr^(-1)*snr0), l+1, 'upper')*factorial(l)...
                /(1+eps_s/snr*snr0)^(l+1);                
        end
        Q5 = Q5 +exp(-eps0)-exp(-tau_bar) + exp(-tau_n - (tau_bar-tau_n)*(1+eps_s*snr^(-1)*snr0))/(1+eps_s*snr^(-1)*snr0);
        suma1 = suma1 + factorial(M)/factorial(n)/factorial(M-n)*exp(-n*tau)*...
            (1-exp(-tau))^(M-n)* Q5;        
    end
    pa(k) = 1 - suma1 - (1-exp(-tau))^M*exp(-eps0); 
      
   
    %analytical DCC
%     suma1 = 0;
%       eps1 = 2^(Ri)-1; 
%     for n=1:M
%         Q7x =0;
%         thetah = max(eps0 , snr/snr0/eps1*tau-1/snr0);
%         for p=0:n
%              Q7x = Q7x + factorial(n)/factorial(p)/factorial(n-p)*(-1)^(p)*exp(p*tau)...
%                  *exp(-p/snr*eps1) *exp(-(1+p/snr*eps1*snr0)*thetah)/(1+p/snr*eps1*snr0);
%         end
%         if snr/snr0/eps1*tau -1/snr0>eps0 
%             Q7 = exp(-thetah) - Q7x + exp(-eps0 ) ...
%                 - exp(- (snr/snr0/eps1*tau - 1/snr0));
%         else
%             Q7 = exp(-thetah) -  Q7x ;
%         end
%         suma1 = suma1 + factorial(M)/factorial(n)/factorial(M-n)*exp(-n*tau)*...
%             (1-exp(-tau))^(M-n)* Q7;        
%     end
%     pa2(k) = 1-suma1 - (1- exp(-tau))^M*exp(-eps0);  
    suma1 = 0;
      eps1 = 2^(Ri)-1; 
        Q4x =0;
        thetah = max(eps0 , snr/snr0/eps1*tau-1/snr0);
        for p=0:M
             Q4x = Q4x + factorial(M)/factorial(p)/factorial(M-p)*(-1)^(p)...
                 *exp(-p/snr*eps1) *exp(-(1+p/snr*eps1*snr0)*thetah)/(1+p/snr*eps1*snr0);
        end
        if snr/snr0/eps1*tau -1/snr0>eps0 
            Q4 = exp(-thetah) - Q4x + (1-(1-exp(-tau))^M)*(exp(-eps0 ) ...
                - exp(- (snr/snr0/eps1*tau - 1/snr0)));
        else
            Q4 = exp(-thetah) -  Q4x ;
        end     

        pa2(k) = 1-Q4 - (1- exp(-tau))^M*exp(-eps0);  

 end
semilogy(snrdb,p3,snrdb, 1-exp(-eps0)*ones(1,length(snrdb)),snrdb,p1,'-s',snrdb,pa,'-*',snrdb,p2,'-d',snrdb,pa2,'-x')% ,snrdb,  (10.^snrdb/10).^(-1)  )