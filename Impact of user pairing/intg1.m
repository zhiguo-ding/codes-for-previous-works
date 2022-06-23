function sumxx = intg1(w2,n,m,i,snr,M)

sumxx=0;
stepy=(w2-sqrt(1+w2)+1)/1000;
yx = [sqrt(1+w2)-1:stepy :w2];

for j = 1 : length(yx)
    y = yx(j);
    fy = exp(-y/snr)/snr;
    Fy = 1 - exp(-y/snr);
    tx=(w2-y)/(1+y);
    sumxx = sumxx + fy*Fy^(n-1-m-i)*(1-Fy)^(M-n) * ( Fy^(m+i) - (1 - exp(-tx/snr))^(m+i) )*stepy;
end
