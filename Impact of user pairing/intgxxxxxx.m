function sumxx = intgxxxxxx(w2,n,m,i,snr,M,b,eps1,a)

sumxx=0;
dx=10;
stepy=(dx - max(b+a*eps1,a*eps1))/500000;
yx = [max(b+a*eps1,a*eps1): stepy: dx];

for j = 1 : length(yx)
    y = yx(j);
    fy = exp(-y);
    Fy = 1 - exp(-y);
    Fyy = 1 - exp(-y);
    Fb = 1 - exp(-b);
    sumxx = sumxx + fy*Fy^(n-1-m-i)*(1-Fy)^(M-n) * ( Fyy^(m+i) - Fb^(m+i) )/(m+i)*stepy;
end
