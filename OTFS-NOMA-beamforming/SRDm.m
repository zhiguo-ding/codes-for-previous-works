function [random_result,SCO,SDR] = SRDm(V,ct,R0) 
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

     for i = 1 : 10

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
                for m0 = 1 : 1
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
               t>0;
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
                   inv_pos(xh(1,1)) <=t;
                   sum(inv_pos( xh(:,mi) )) <=eta/N;
               end               
               trace(X)<=1;
        cvx_end
      [e1,e2]= eig(X);
      wo = e1(:,V);
      if abs(e2(V-1,V-1))>0.1
          dfere=0;
      end
%       
%                  for mo = 1 : M
%                obj0 = obj0 + log(1+real(snr*trace(X*h(:,mo,mo)*h(:,mo,mo)')));
%            end

     w0=zeros(V,1);
     w0(w0_index)=1;
    for mc =1 : 1
        t1(mc) = log2(1+snr*(w'*h(:,mc,mc)*h(:,mc,mc)'*w));
        t2(mc) = log2(1+ snr*(w0'*h(:,mc,mc)*h(:,mc,mc)'*w0));
        t3(mc) = log2(1+snr*(wo'*h(:,mc,mc)*h(:,mc,mc)'*wo)); 
    end
    SCOx(ict) = min(t1);
    random_resultx(ict) = min(t2);
     SDRx(ict) =    min(t3);
     
     if abs(SCOx(ict) - SDRx(ict))>0.001
         dfdfe=1;
     end
  
 end
 SCO(knr) = mean(SCOx);
 random_result(knr) = mean(random_resultx);
 SDR(knr) = mean(SDRx);
 
 
 %%%%%%%%%%%%%%%%% WRONG the minimal 
end

  
 