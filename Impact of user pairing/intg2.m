function sumxx = intg2(w2,n,m,i,snr,M)

sumxx=0;
stepy=(w2-sqrt(1+w2)+1)/1000;
yx = [sqrt(1+w2)-1:stepy :w2];

for j = 1 : length(yx)
    y = yx(j);
    tx=(w2-y)/(1+y);
    sumxx = sumxx + y^(n-1-m-i) * ( y^(m+i) - tx^(m+i) )*stepy;
end
