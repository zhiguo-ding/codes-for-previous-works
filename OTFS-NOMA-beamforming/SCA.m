function [random_result, my_result]=SCA(eps,ct,R0,V)
 %load('data_channel_2D.mat') 
 N=8;
 M=8;
 %V=8;
 
 %eps = 0.0075;
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

   %CVX part
     sumx1=0;
     w0=zeros(V,1);
     w0_index = randsrc(1,1,[1:1:V]);
     w0(w0_index)=1;     
     ixxx=10;
     for i = 1 : ixxx
         cvx_begin quiet
           variables t w(V) 
           %dual variables z0(M) z1 z2 z3 %y corresponds to nu and z correspons to lambda
           %calculate the channel estimation error constraint
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
              for m0 = 1 : M
                    xdd(:,m0)= fkl_first(w0,h(:,m0,m0),0);
              end
     if max(max(isnan(abs(xdd))))==0
       minimize( t )
           subject to
              for m0 = 1 : M
                    real(1/(w0'*h(:,m0,m0)*h(:,m0,m0)'*w0) + fkl_first(w0,h(:,m0,m0),0)'*(w-w0)) <=t;
              end
                real(sum13 + sum11'*(w-w0))  <=eta; 
              for mi = 1 : M
                  real(sum23(mi) + sum21(:,mi)'*(w-w0))  <=eta/N;
              end
                w'*w<=1;
        cvx_end
     else
         w=zeros(V,1);
            i=ixxx;
     end
        if max(isnan(w))==1
            w=zeros(V,1);
            i=ixxx;
        else
            test_o(i) = 1/(w'*h(:,1,1)*h(:,1,1)'*w);
            w0 = w;
        end
     end
     wtr=w;
     
     if w'*w>0.01
         %%% we check whether the resulting w causes a zero, STRONG ERROR
         w0=zeros(V,1);
         w0(w0_index)=1;
         error_check = zeros(N,M);
               for k = 1: N
                 for l = 1 : M
                     cvx_begin quiet
                         variables e(V,1) 
                         minimize( real(e'*e) )
                         subject to
                            real(w*w'*e) ==-real(w*w'*g(:,k,l));
                     cvx_end
                     if max(isnan(e))~=1 %there is a feasible e
                         if e'*e<eps^2 %this feasible e also satisifies the eps constraint 
                            error_check(k,l) = 1;
                         end
                     end
                 end
               end
               if max(max(error_check))>0
                   wtr=zeros(V,1); %there is strong error
                   w0=zeros(V,1);
               end
     else
                   wtr=zeros(V,1); %there is strong error
                   w0=zeros(V,1);
     end
     
     for mc =1 : M
        t1(mc) = log2(1+snr*(w0'*h(:,mc,mc)*h(:,mc,mc)'*w0));
        t2(mc) = log2(1+ snr*(wtr'*h(:,mc,mc)*h(:,mc,mc)'*wtr));
    end
        
    my_resultx(ict) = min(t2);
    random_resultx(ict) = min(t1);
  
 end
 my_result(knr) = mean(my_resultx);
 random_result(knr) = mean(random_resultx);
end

%plot(snrdb,random_result, snrdb,my_result)
 
   %for i = 1 : M
    
  
 