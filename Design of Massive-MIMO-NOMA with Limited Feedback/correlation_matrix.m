clear all;
M=50; Delta=15/360*2*pi;lamd=1/(2*10^10); 
D=0.5/sqrt( (1 - cos(2*pi/M))^2 +(sin(2*pi/M))^2);
stepd = Delta/100;
dd = [-Delta : stepd : Delta];

K=4;

for k = 1 : K
    thetag = -pi+Delta +(k-1)*2*pi/K; 
    for m =1 : M
        B(1:2,m) = [lamd*D*cos((m-1)*2*pi/M) lamd*D*sin((m-1)*2*pi/M)];
    end
    for m = 1: M
        for n = 1 :M
            sum1 = 0;
            for p = 1 : length(dd)
                al = dd(p);
                kx = -2*pi/lamd*[cos(al+thetag) ; sin(al+thetag)];
                sum1 = sum1 + 1/2/Delta* exp(sqrt(-1)*kx'*(B(:,m) -B(:,n))  )*stepd;
            end
            Rr(m,n,k)=sum1;
        end
    end

    r=rank(Rr(:,:,k)); Mt = M-r*(K-1);
    [a1,b1 ] = eig(Rr(:,:,k)); a11=a1'; U(:,:,k) = a11(M-r+1:M,:); Lam(:,:,k) = (b1(M-r+1:M,M-r+1:M))^(1/2) ;
rx(k)=rank(Rr(:,:,k)); 
end

Ux = zeros(r*(K-1),M);
for k = 2 :K
    Ux((k-2)*r+1:(k-1)*r,:) = U(:,:,k);
end
P = null(Ux); 