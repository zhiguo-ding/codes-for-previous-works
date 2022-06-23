function sumxx = intg3(w2,n,m,i,snr,M,b,eps1,a)

sumxx=0;
stepy=(b+a*eps1-a*eps1)/1000;
yx = [a*eps1: stepy: b+a*eps1];

for j = 1 : length(yx)
    y = yx(j);
    gy = exp(-y);
    Gy = 1 - exp(-y);
    Gb = 1 - exp(-b);
    sumxx = sumxx + gy*Gy^(n-1-m-i)*(1-Gy)^(M-n) * ( Gy^(m+i) - Gb^(m+i) )/(m+i)*stepy;
end
