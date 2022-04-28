function [random_result,SCO,SDR]=SRD_Gaus(V,ct,R0,ixxx)
 %load('data_channel_2D.mat') 
 N=8;
 M=8;
 
 eps = 0;
 snrdb = [10: 5 : 25];
 
for knr = 1 : length(snrdb) 
    snr=10^(snrdb(knr)/10);
    eta=(1/(2^R0-1)-1)*snr*N*M;
 for ict = 1 : ct
     % generate channel
      for k = 1: N
         for l = 1 : M
            h(:,k,l) = complex(sqrt(0.5)*randn(V,1),sqrt(0.5)*randn(V,1));
            g(:,k,l) = complex(sqrt(0.5)*randn(V,1),sqrt(0.5)*randn(V,1));
         end
      end

      %%%% we remove the initalization and use random one. Please go to 

 
     sumx1=0;
     w0=zeros(V,1);
     w0_index = randsrc(1,1,[1:1:V]);
     w0(w0_index)=1;

     %ixxx=10;
     for i = 1 : ixxx

         cvx_begin quiet
           variables w(V)  t
           %dual variables  z1 z2 z3 %y corresponds to nu and z correspons to lambda
           %calculate the second constraint
           sum11 = 0; sum12=0; sum13=0;
           for k = 1: N
             for l = 1 : M
                 sum11 = sum11 + fkl_first(w0,g(:,k,l),eps);%*(w-w0)...
                 sum13 = sum13 + 1/(sqrt(w0'*g(:,k,l)*g(:,k,l)'*w0) - sqrt(eps*w0'*w0))^2;
             end
           end
           %calculate the second constraint
           sum21 = zeros(V,M); sum22=0; sum23=zeros(M,1);
           for k = 1: M %%%%%%%%% PLEASE note that the dimension of H is MXM, NOT NXM
             for l = 1 : M
                 sum21(:,k) = sum21(:,k) + fkl_first(w0,h(:,k,l),0);%*(w-w0)...
                 sum23(k) = sum23(k) + 1/(w0'*h(:,k,l)*h(:,k,l)'*w0);
             end
           end

           %start the cvx
           minimize( t)
           subject to
                for m0 = 1 : M
                    real(1/(w0'*h(:,m0,m0)*h(:,m0,m0)'*w0)) +real(fkl_first(w0,h(:,m0,m0),0)'*(w-w0)) <=t;
                end
                real(sum13) + real(sum11'*(w-w0))  <=eta; 
                for mi = 1 : M
                    real(sum23(mi) + sum21(:,mi)'*(w-w0))  <=eta/N;
                end
                w'*w<=1;
        cvx_end
        if max(isnan(w))==1
            w=w0;
            break            
        else
            w0 = w;
        end
     end

     %SDR
      cvx_begin quiet
           variables X(V,V)  xh(M,M) yg(N,M) t
           dual variables tx
           %start the cvx
           minimize( t )
           subject to
               t>=0;
               tx: X <In> semidefinite(V);  
               for mi = 1 : M
                   for mj = 1 : M
                       real(trace(X*h(:,mi,mj)*h(:,mi,mj)')) == xh(mi,mj);
                       xh(mi,mj) >=0;
                   end
               end
               for kn = 1 : N
                   for lm = 1 : M
                       real(trace(X*g(:,kn,lm)*g(:,kn,lm)')) == yg(kn,lm);
                       yg(kn,lm)>=0;
                   end
               end
               sum(sum(inv_pos(yg))) <= eta;
               for mi = 1 : M
                   inv_pos(xh(mi,mi)) <=t;
                   sum(inv_pos( xh(:,mi) )) <=eta/N;
               end               
               trace(X)<=1;
        cvx_end
        if max(max(isnan(X)))==1
            wo = zeros(V,1);
        else
            [e1,e2]= eig(X);
            for li = 1 : ixxx
                xl(:,li) = randn(V,1);
                wlx(:,li) = e1*sqrt(e2)*xl(:,li);
                   for mj = 1 : M
                       test_k1(mj)=log2(1+snr*(wlx(:,li)'*h(:,mj,mj)*h(:,mj,mj)'*wlx(:,li)));%real(trace(wlx(:,li)*wlx(:,li)'*h(:,mj,mj)*h(:,mj,mj)'));
                   end
                   test_k2(li) = min(test_k1);
            normx(li) = wlx(:,li)'*wlx(:,li);
            end
            [tere, indexw] = max(test_k2);
            wo = wlx(:,indexw); 
         end
        
      %wo = e1(:,V);
      
     w0=zeros(V,1);
     w0(w0_index)=1;
    for mc =1 : M
        t1(mc) = log2(1+snr*(w'*h(:,mc,mc)*h(:,mc,mc)'*w));
        t2(mc) = log2(1+ snr*(w0'*h(:,mc,mc)*h(:,mc,mc)'*w0));
        t3(mc) = log2(1+snr*(wo'*h(:,mc,mc)*h(:,mc,mc)'*wo)); 
    end
    SCOx(ict) = min(t1);
    random_resultx(ict) = min(t2);
     SDRx(ict) =    min(t3);
      
   
 end
 SCO(knr) = mean(SCOx);
 random_result(knr) = mean(random_resultx);
 SDR(knr) = mean(SDRx);
  
 
 %%%%%%%%%%%%%%%%% WRONG the minimal 
end

%plot(snrdb,random_result,snrdb,SCO,snrdb,SDR)
 
 